#include <igl/eigs.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/read_triangle_mesh.h>
#include <Eigen/Sparse>
#include <iostream>
#include <vector>
#include <set>

//Standard function to compute the cotan from two vectors.
double cotan(Eigen::Vector3d a, Eigen::Vector3d b) {
    return a.dot(b) / (a.cross(b)).norm();
}

//The name is pretty explicit
void computeCotangentLaplacian(Eigen::MatrixXd V, Eigen::MatrixXi F, Eigen::SparseMatrix<double>& cotL) {
    cotL.resize(V.rows(), V.rows());
    cotL.data().squeeze();
    
    //We travel each face and we compute the cotan between each edge starting from the vertex not included inside the triangle.
    for (int fi = 0; fi < F.rows(); fi++) {
        Eigen::MatrixXi faces_vertices = F.row(fi);
        for (int ei = 0; ei < faces_vertices.cols(); ei++) {
            int i, j, k;
            i = faces_vertices(ei);
            j = faces_vertices((ei + 1) % 3);
            k = faces_vertices((ei + 2) % 3);
            double w = cotan(Eigen::Vector3d(V.row(k) - V.row(i)), Eigen::Vector3d(V.row(k) - V.row(j)));

            //The cotan angle will only be add two times because the edge i,j is only common to two faces.
            cotL.coeffRef(i, j) -= w;
            cotL.coeffRef(j, i) -= w;

            //For each edge of this triangle, the vertice index will be present two times.
            //The angle for the edge i-j will be the same value as alpha for i and as beta for j.
            //It is easier to imagine with a little sketch.
            cotL.coeffRef(i, i) += w;
            cotL.coeffRef(j, j) += w;
        }
    }
    cotL.makeCompressed();
    cotL = cotL / 2.;
}

//The name is pretty explicit
void computeMassMatrix(Eigen::MatrixXd V, Eigen::MatrixXi F, Eigen::SparseMatrix<double>& massMatrix) {
    massMatrix.resize(V.rows(), V.rows());
    massMatrix.reserve(Eigen::VectorXi::Constant(V.rows(), 10));
    //We travel each face. The following is similar as computeCotLaplacian.
    for (int fi = 0; fi < F.rows(); fi++) {
        Eigen::MatrixXi faces_vertices = F.row(fi);
        //We travel each edge of the face.
        for (int ei = 0; ei < faces_vertices.cols(); ei++) {
            int i, j, k;
            i = faces_vertices(ei);
            j = faces_vertices((ei + 1) % 3);
            k = faces_vertices((ei + 2) % 3);

            //We compote the area of the triangle with vector.
            double area = 1 / .2 * (Eigen::Vector3d(V.row(k) - V.row(i)).cross(Eigen::Vector3d(V.row(k) - V.row(j)))).norm();

            //The edge ij is only present two times because there is only two triangles adjacent to this edge. Then, at the end, the sum will be the adjacent area.
            massMatrix.coeffRef(i, j) += area / 12.;
            massMatrix.coeffRef(j, i) += area / 12.;

            //For each vertex of the face, we add the area of the triangle one time. At the end, we will have the one-ring area because each triangle area with vertex i will be add up.
            massMatrix.coeffRef(i, i) += area / 6.;
        }
    }
    massMatrix.makeCompressed();
}

void computeLumpedMassMatrix(Eigen::MatrixXd V, Eigen::MatrixXi F, Eigen::SparseMatrix<double>& lumpedMassMatrix) {
    computeMassMatrix(V, F, lumpedMassMatrix);

    //The lumped mass matrix is equal to the diagonal of the massmatrix times 2 because it correspond to the third of the one-ring area for each vertex.
    lumpedMassMatrix = Eigen::SparseMatrix<double>((lumpedMassMatrix.diagonal() * 2.).asDiagonal());
    lumpedMassMatrix.makeCompressed();
}

int main(int argc, char * argv[])
{
    std::string model_path = argc > 1 ? argv[1] : "../../input/camel-1.obj";
    std::string method_name = argc > 2 ? argv[2] : "default";
    int neigens = argc > 3 ? std::stoi(argv[3]) : 3;

    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    Eigen::MatrixXd eigenvectors;
    igl::read_triangle_mesh(model_path, V, F);

    std::cout << V.rows() << std::endl;
    std::cout << F.rows() << std::endl;

    // Compute Cotangant Laplacian matrix and Mass matrix
    Eigen::SparseMatrix<double> cotL,massMatrix;
    computeCotangentLaplacian(V, F,cotL);
    if (method_name == "default") {
        computeMassMatrix(V, F, massMatrix);
    }
    else if  (method_name == "lumped") {
        computeLumpedMassMatrix(V, F, massMatrix);
    }

    //Compute the eigenvalues and eigenvectors.
    Eigen::MatrixXd eVec;
    Eigen::VectorXd eVal;
    igl::eigs(cotL, massMatrix, neigens, igl::EIGS_TYPE_SM, eVec, eVal);
    std::cout << "Eigenvalues: " << eVal << std::endl;

    //Verification of the properties for the cotangiant laplacian and mass matrix.
    Eigen::SparseMatrix<double> iglcotL;
    igl::cotmatrix(V, F, iglcotL);
    std::cout << "Cotangiant Laplacian -> Sparse (NonZeros/(Cols * Rows)): " << cotL.nonZeros() << "/" << cotL.cols() * cotL.rows() << std::endl;
    std::cout << "Cotangiant Laplacian -> Symmetric: " << (cotL.transpose().isApprox(cotL)) << std::endl;
    std::cout << "Cotangiant Laplacian -> Non-negative eigenvalues (Semite Definite Positive): " << (eVal.array() >= -0.01).all() << std::endl;
    std::cout << "Cotangiant Laplacian -> Sum of row to 0: " << (Eigen::MatrixXd(cotL * Eigen::VectorXd::Ones(cotL.cols())).cwiseAbs().array() <= 0.1).all() << std::endl;
    std::cout << "Cotangiant Laplacian -> Difference with igl [Note: The cotmatrix of igl is the same as our by a factor of -1. They adopt a different convention.]: " << (Eigen::MatrixXd(iglcotL + cotL).cwiseAbs().array() < 0.001).all() << std::endl;

    if (method_name == "default") {
        std::cout << method_name << " Mass Matrix -> Sum of row must be equal to one-ring area/3: " << ((massMatrix * Eigen::VectorXd::Ones(massMatrix.cols()) - massMatrix.diagonal() * 2.).cwiseAbs().array() < 0.001).all() << std::endl;
    }
    else if (method_name == "lumped") {
        Eigen::SparseMatrix<double> iglmassMatrix;
        igl::massmatrix(V, F, igl::MassMatrixType::MASSMATRIX_TYPE_BARYCENTRIC, iglmassMatrix);
        iglmassMatrix *= 10;
        std::cout << method_name << " Mass Matrix -> Difference with igl [Note: The lumped massmatrix of igl is the same as our by a factor of 10. They adopt a different convention.]: " << (Eigen::MatrixXd(iglmassMatrix - massMatrix).cwiseAbs().array() < 0.001).all() << std::endl;
    }

    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(V,F);
    viewer.data().compute_normals();

    int selectedcolumn=0;
    viewer.callback_key_down = [&](igl::opengl::glfw::Viewer & viewer,unsigned char key,int)->bool
    {
        switch(key)
        {
            default:
                return false;
            case ' ':
            {
                selectedcolumn = (eVec.cols() + selectedcolumn - 1) % eVec.cols();
                
                viewer.data().set_data(eVec.col(selectedcolumn));
                return true;
            }
        }
    };
    viewer.callback_key_down(viewer,' ',0);
    viewer.data().show_lines = false;
    viewer.launch();
}
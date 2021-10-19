#include <igl/eigs.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/read_triangle_mesh.h>
#include <Eigen/Sparse>
#include <iostream>
#include <vector>
#include <set>

#include<stdio.h>
#include<stdlib.h>

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
    double dt = argc > 2 ? strtod(argv[2],NULL) : 0.01;
    int nbIteration = argc > 3 ? strtol(argv[3], NULL,0) : 1000;


    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    Eigen::MatrixXd eigenvectors;

    igl::read_triangle_mesh(model_path, V, F);
    std::cout << V.rows() << std::endl;
    std::cout << F.rows() << std::endl;

    //1. We compute the Cotangent Laplacian and Lumped Mass Matrix
    Eigen::SparseMatrix<double> cotL, lumpedMassMatrix, lumpedMassMatrixInverted;
    computeCotangentLaplacian(V, F, cotL);
    computeLumpedMassMatrix(V, F, lumpedMassMatrix);
    lumpedMassMatrixInverted = lumpedMassMatrix.cwiseInverse();

    //2. We compute a subset of the eigenvectors.
    Eigen::MatrixXd eVec;
    Eigen::VectorXd eVal;
    igl::eigs(cotL, lumpedMassMatrix, 3, igl::EIGS_TYPE_SM, eVec, eVal);

    //Compute Heat Equation
    double lambda = 1;
    Eigen::SparseMatrix<double> I(V.rows(),V.rows());
    I.setIdentity();

    std::vector<Eigen::MatrixXd> U;
    U.resize(nbIteration + 1);

    U[0] = V;

    Eigen::SparseMatrix<double> deltaIntegrationStep;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> sLDLT;

    deltaIntegrationStep = ( (3.* lumpedMassMatrix) + (lambda * cotL * 0.5 * dt)); // Those constant has been added to have the same definition as the paper.
    deltaIntegrationStep.makeCompressed();
    sLDLT.compute(deltaIntegrationStep);

    for (int step = 0; step < U.size() - 1; step++) {
        //We solve ( 3*M + lambda * cotL * 1/2 * \Delta t) u_{n+1} = (3*lumpedMassMatrix) * u_n
        U[step + 1] = sLDLT.solve(3*lumpedMassMatrix * U[step]);
    }

    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(U[nbIteration - 1], F);
    viewer.data().set_data(eVec.col(2));

    int selectedcolumn=-1; 
    viewer.callback_key_down = [&](igl::opengl::glfw::Viewer & viewer,unsigned char key,int)->bool
    {
        switch(key)
        {
            default:
                return false;
            case ' ':
            {
                selectedcolumn++;
                if (selectedcolumn%2){
                    viewer.data().set_vertices(U[0]);
                }
                else
                {
                   viewer.data().set_vertices(U[nbIteration - 1]);
                }

                return true;
            }
        }
    };
    viewer.callback_key_down(viewer,' ',0);
    viewer.data().show_lines = false;
    viewer.launch();
}

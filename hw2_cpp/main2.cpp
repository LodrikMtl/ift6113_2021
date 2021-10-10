#include <igl/eigs.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/read_triangle_mesh.h>
#include <Eigen/Sparse>
#include <iostream>
#include <vector>
#include <set>

double cotan(Eigen::Vector3d a, Eigen::Vector3d b){
    return a.dot(b) / (a.cross(b)).norm();
}


Eigen::SparseMatrix<double> computeCotangentLaplacian(Eigen::MatrixXd V, Eigen::MatrixXi F){
    Eigen::SparseMatrix<double> cotL(V.rows(),V.rows());
    cotL.reserve(Eigen::VectorXi::Constant(V.rows(), 14));
    for (int fi = 0; fi < F.rows(); fi++) {
        Eigen::MatrixXi faces_vertices = F.row(fi);
        for (int ei = 0; ei < faces_vertices.cols(); ei++) {
            int i,j,k;
            i = faces_vertices(ei);
            j = faces_vertices((ei + 1) % 3);
            k = faces_vertices((ei + 2) % 3);
            double w = cotan(Eigen::Vector3d(V.row(k) - V.row(i)),Eigen::Vector3d(V.row(k) - V.row(j)));

            //The cotan angle will only be add two times because the edge i,j is only connected by two faces.
            cotL.coeffRef(i, j) -= w;
            cotL.coeffRef(j, i) -= w;

            //For each edge of the triangles, the vertice index will be present two times.
            //The angle for the edge i-j will be the same value as alpha for i and as beta for j.
            //It is easier to imagine with a little sketch.
            cotL.coeffRef(i, i) += w;
            cotL.coeffRef(j, j) += w;
        }
    }
    cotL.makeCompressed();
    return cotL/2.;
}

Eigen::SparseMatrix<double> computeMassMatrix(Eigen::MatrixXd V, Eigen::MatrixXi F) {
    Eigen::SparseMatrix<double> massMatrix(V.rows(), V.rows());
    massMatrix.reserve(Eigen::VectorXi::Constant(V.rows(), 14));
    for (int fi = 0; fi < F.rows(); fi++) {
        Eigen::MatrixXi faces_vertices = F.row(fi);
        for (int ei = 0; ei < faces_vertices.cols(); ei++) {
            int i, j, k;
            i = faces_vertices(ei);
            j = faces_vertices((ei + 1) % 3);
            k = faces_vertices((ei + 2) % 3);
            double area = 1/.2 * (Eigen::Vector3d(V.row(k) - V.row(i)).cross(Eigen::Vector3d(V.row(k) - V.row(j)))).norm();

            //The cotan angle will only be add two times because the edge i,j is only connected by two faces.
            massMatrix.coeffRef(i, j) += area/12.;
            massMatrix.coeffRef(j, i) += area/12.;

            //For each edge of the triangles, the vertice index will be present two times.
            //The angle for the edge i-j will be the same value as alpha for i and as beta for j.
            //It is easier to imagine with a little sketch.
            massMatrix.coeffRef(i, i) += area/6.;
        }
    }
    massMatrix.makeCompressed();
    return massMatrix / 2.;
}


int main(int argc, char * argv[])
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    Eigen::MatrixXd eigenvectors;
    if(!igl::read_triangle_mesh(
            argc>1?argv[1]: "../../input/camel-1.obj",V,F))
    {
        std::cout<<"failed to load mesh"<<std::endl;
    }
    std::cout << V.rows() << std::endl;
    std::cout << F.rows() << std::endl;

    //1.
    Eigen::SparseMatrix<double> cotL = computeCotangentLaplacian(V, F);
    Eigen::SparseMatrix<double> massMatrix = computeMassMatrix(V, F);
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> sol;
    sol.compute(massMatrix);
    Eigen::SparseMatrix<double> I(massMatrix.cols(),massMatrix.rows());
    I.setIdentity();
    Eigen::SparseMatrix<double> invMassMatrix = sol.solve(I);
    std::cout << (massMatrix * invMassMatrix).isApprox(I) << std::endl;

    //2
    Eigen::MatrixXd eVec;
    Eigen::VectorXd eVal;
    igl::eigs(cotL, massMatrix, 3, igl::EIGS_TYPE_SM, eVec, eVal);

     //1.1 Do some checkup
    std::cout << "Sparse (NonZeros/(Cols * Rows)): " << cotL.nonZeros() << "/" << cotL.cols() * cotL.rows() << std::endl;
    std::cout << "Symmetric: " << (cotL.transpose().isApprox(cotL)) << std::endl;
    std::cout << "Non-negative (Semite Definite Positive): " << (eVal.array() >= 0).all() << std::endl;

    //4.1 Verify it's property
    std::cout << "Each sum of row must be equal to one-ring area/3: " << ((massMatrix * Eigen::VectorXd::Ones(massMatrix.cols()) - massMatrix.diagonal() * 2).array() < 0.001).all() << std::endl;

    //Compute Heat Equation
    double dt = 1;
    double c = 1;
    Eigen::MatrixXd U(V.rows(), 10);
    U.col(0) = Eigen::VectorXd::Random(V.rows()) + Eigen::VectorXd::Ones(V.rows());
    for (int step = 0; step < 9; step++) {
        U.col(step + 1) = c * invMassMatrix * cotL * U.col(step) * dt + U.col(step);
        std::cout << "Step: " << step + 1 << "Min: " << U.col(step + 1).minCoeff() << "Max: " << U.col(step + 1).maxCoeff() << "Average: " << U.col(step + 1).mean() << std::endl;
    }

    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(V,F);
    viewer.data().compute_normals();

    //3.
    int selectedcolumn=-1; // it will be updated to 0 once we run viewer.callback_key_down(...)
    viewer.callback_key_down = [&](igl::opengl::glfw::Viewer & viewer,unsigned char key,int)->bool
    {
        switch(key)
        {
            default:
                return false;
            case ' ':
            {
                selectedcolumn = (selectedcolumn + 1) % U.cols();

                viewer.data().set_data(U.col(selectedcolumn));
                return true;
            }
        }
    };
    viewer.callback_key_down(viewer,' ',0);
    viewer.data().show_lines = false;
    viewer.launch();
}

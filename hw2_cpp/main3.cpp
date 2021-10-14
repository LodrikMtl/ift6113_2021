#include <igl/eigs.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/read_triangle_mesh.h>
#include <Eigen/Sparse>
#include <iostream>
#include <vector>
#include <set>

#include<stdio.h>
#include<stdlib.h>

double cotan(Eigen::Vector3d a, Eigen::Vector3d b) {
    return a.dot(b) / (a.cross(b)).norm();
}

void computeCotangentLaplacian(Eigen::MatrixXd V, Eigen::MatrixXi F, Eigen::SparseMatrix<double>& cotL) {
    cotL.resize(V.rows(), V.rows());
    cotL.data().squeeze();
    //cotL.reserve(Eigen::VectorXi::Constant(V.rows(), 14));
    for (int fi = 0; fi < F.rows(); fi++) {
        Eigen::MatrixXi faces_vertices = F.row(fi);
        for (int ei = 0; ei < faces_vertices.cols(); ei++) {
            int i, j, k;
            i = faces_vertices(ei);
            j = faces_vertices((ei + 1) % 3);
            k = faces_vertices((ei + 2) % 3);
            double w = cotan(Eigen::Vector3d(V.row(k) - V.row(i)), Eigen::Vector3d(V.row(k) - V.row(j)));

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
    cotL = cotL / 2.;
}

void computeMassMatrix(Eigen::MatrixXd V, Eigen::MatrixXi F, Eigen::SparseMatrix<double>& massMatrix) {
    massMatrix.resize(V.rows(), V.rows());
    massMatrix.reserve(Eigen::VectorXi::Constant(V.rows(), 10));
    for (int fi = 0; fi < F.rows(); fi++) {
        Eigen::MatrixXi faces_vertices = F.row(fi);
        for (int ei = 0; ei < faces_vertices.cols(); ei++) {
            int i, j, k;
            i = faces_vertices(ei);
            j = faces_vertices((ei + 1) % 3);
            k = faces_vertices((ei + 2) % 3);
            double area = 1 / .2 * (Eigen::Vector3d(V.row(k) - V.row(i)).cross(Eigen::Vector3d(V.row(k) - V.row(j)))).norm();

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
}

void computeLumpedMassMatrix(Eigen::MatrixXd V, Eigen::MatrixXi F, Eigen::SparseMatrix<double>& lumpedMassMatrix) {
    computeMassMatrix(V, F, lumpedMassMatrix);

    lumpedMassMatrix = Eigen::SparseMatrix<double>((lumpedMassMatrix.diagonal() * 2.).asDiagonal());
    lumpedMassMatrix.makeCompressed();
}

int main(int argc, char * argv[])
{
    std::string input_path = argc > 1 ? argv[1] : "../../input/camel-1.obj";
    double dt = argc > 2 ? strtod(argv[2],NULL) : 0.01;
    int iteration = argc > 3 ? strtol(argv[3], NULL,0) : 1000;


    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    Eigen::MatrixXd eigenvectors;

    igl::read_triangle_mesh(input_path, V, F);
    std::cout << V.rows() << std::endl;
    std::cout << F.rows() << std::endl;

    //1.
    Eigen::SparseMatrix<double> cotL, massMatrix, lumpedMassMatrix, lumpedMassMatrixInverted;
    computeCotangentLaplacian(V, F, cotL);
    computeLumpedMassMatrix(V, F, lumpedMassMatrix);

    //2
    Eigen::MatrixXd eVec;
    Eigen::VectorXd eVal;
    igl::eigs(cotL, lumpedMassMatrix, 3, igl::EIGS_TYPE_SM, eVec, eVal);

    //Compute Heat Equation
    double lambda = 1;
    Eigen::SparseMatrix<double> I(V.rows(),V.rows());
    I.setIdentity();

    std::vector<Eigen::MatrixXd> U;// f
    U.resize(iteration);

    U[0] = V;

    Eigen::SparseMatrix<double> deltaIntegrationStep;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> sLDLT;

    deltaIntegrationStep = ( (3.* lumpedMassMatrix) - (lambda * cotL * 0.5 * dt));
    deltaIntegrationStep.makeCompressed();
    sLDLT.compute(deltaIntegrationStep);

    for (int step = 0; step < iteration - 1; step++) {
        U[step + 1] = sLDLT.solve(3*lumpedMassMatrix * U[step]);
    }

    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(U[iteration - 1], F);
    viewer.data().set_data(eVec.col(1));

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
                selectedcolumn++;
                if (selectedcolumn%2){
                    viewer.data().set_mesh(U[0], F);
                }
                else
                {
                   //viewer.data().set_mesh(U[iteration - 1], F);
                }

                return true;
            }
        }
    };
    viewer.callback_key_down(viewer,' ',0);
    viewer.data().show_lines = false;
    viewer.launch();
}

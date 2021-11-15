#include <igl/eigs.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/parula.h>
#include <igl/read_triangle_mesh.h>
#include <igl/readDMAT.h>
#include <igl/slice_into.h>
#include <igl/harmonic.h>
#include <igl/ismember.h>
#include <igl/setxor.h>
#include <igl/colon.h>
#include <Eigen/Sparse>
#include "Eigen/Core"
#include<Eigen/SparseQR>	
#include <iostream>
#include <fstream>
#include <queue>

int main(int argc, char *argv[])
{
    {
        using namespace std;

        std::string model_path = argc > 1 ? argv[1] : "../../input/bar2.off";
        bool use_igl_solution = argc > 2 ? stoi(argv[2]) : 0;

        Eigen::MatrixXd V, rhsD, D,U;
        Eigen::MatrixXi F;

        Eigen::MatrixXd hV;
        Eigen::VectorXi hI;

        size_t botDirPos = model_path.find_last_of("/");
        size_t botExtPos = model_path.find_last_of(".");
        std::string filename = model_path.substr(botDirPos+1, model_path.size());
        igl::readDMAT("../../input/fixedindices-" + filename + ".txt", hI);
        igl::readDMAT("../../input/fixedpositions-" + filename + ".txt", hV);
        hV.transposeInPlace();

        igl::read_triangle_mesh(model_path, V, F);

        int nVertices = V.rows();
        int nHandles = hI.rows();

        
        // min J(d) = min d^T K d
        // J is the function to minimize the energy
        // (K + K^T) d = 0 
        // dK d = 0
        Eigen::SparseMatrix<double> dK{ nVertices, nVertices };
        // L is the cotangant laplacian, M is the mass matrix.
        Eigen::SparseMatrix<double> L, M;
        igl::cotmatrix(V, F, L);
        L = (-L).eval();
        igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);
        dK = (L * M.cwiseInverse() * L) + Eigen::SparseMatrix<double>{(L* M.cwiseInverse()* L).transpose()};

        //Build the corresponding rhsD from the handles contraint like discuss in the report.
        rhsD.resize(nVertices, 3);
        rhsD.setZero();
        for (int hi = 0; hi < nHandles; hi++) {
            int vi = hI.coeff(hi);
            rhsD.row(vi) = hV.row(hi) - V.row(vi);
            dK.row(vi) *= 0;
            dK.coeffRef(vi, vi) = 1;
        }
        dK.data().squeeze();

        //Find the corresponding vector that satisfy the contraint discuss above.
        Eigen::SparseLU<Eigen::SparseMatrix<double>> sLU;
        sLU.compute(dK);
        D = sLU.solve(rhsD);
        std::cout << "Did LU converge? " << ((dK * D - rhsD).array() < 0.01).all() << std::endl;

        if (use_igl_solution) {
            D.setZero();
            igl::harmonic(V, F, hI, hV - igl::slice(V,hI,1), 1, D);
        }

        U = D + V;

        int i = 0;
        igl::opengl::glfw::Viewer viewer;
        viewer.callback_key_down = [&](igl::opengl::glfw::Viewer& viewer, unsigned char key, int)->bool
        {
            switch (key)
            {
            default:
                return false;
            case ' ':
            {
                i = (i + 1) % 2;
                if (i == 0) {
                    viewer.data().set_mesh(U, F);
                }
                else if (i == 1) {
                    viewer.data().set_mesh(V, F);
                }
            }
            }
        };
        viewer.data().add_points(hV, Eigen::RowVector3d(1, 0, 0));
        viewer.data().set_mesh(V,F);
        viewer.data().show_lines = false;
        viewer.callback_key_down(viewer, ' ', 0);
        viewer.launch();
    }
}
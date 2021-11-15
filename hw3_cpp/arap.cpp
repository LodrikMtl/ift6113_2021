#include <igl/eigs.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/parula.h>
#include <igl/read_triangle_mesh.h>
#include <igl/readDMAT.h>
#include <igl/slice_into.h>
#include <igl/arap.h>
#include <igl/lbs_matrix.h>
#include <igl/deform_skeleton.h>
#include <igl/dqs.h>
#include <Eigen/Sparse>
#include "Eigen/Core"
#include <iostream>
#include <fstream>
#include <queue>


//Standard function to compute the cotan from two vectors.

int main(int argc, char *argv[])
{
    {
        using namespace std;

        std::string model_path = argc > 1 ? argv[1] : "../../input/bar2.off";
        int nb_iter = argc > 2 ? stoi(argv[2]) : 1;
        bool use_igl_solution = argc > 3 ? stoi(argv[3]): 0;

        Eigen::MatrixXd V,VPrime;
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

        //Make an initial guess from vertices layout 
        VPrime = V;
        for (int hi = 0; hi < nHandles; hi++) {
            int vi = hI.coeff(hi);
            VPrime.row(vi) = hV.row(hi);
        }

        std::vector<std::vector<int>> neighs;
        igl::adjacency_list(F, neighs);

        // W is the Weight Matrix as define as ARAP, L is the cotangant laplaciant
        Eigen::SparseMatrix<double> W,L;
        igl::cotmatrix(V, F, W);

        //Intermediate vertices through each iteration of the solver
        std::vector<Eigen::MatrixXd> Vintermediate;
        Vintermediate.push_back(VPrime);
        
        for (int i = 0; i < nb_iter; i++)
        {
            std::vector<Eigen::MatrixXd> rotations;
            
            //Solving the Local Step
            for (int vi = 0; vi < nVertices; vi++) {

                // We are building the corresponding S_vi = P_vi^T * D_vi * nP_vi as describe in the paper
                Eigen::MatrixXd P_vi{ neighs[vi].size(),3 };
                Eigen::MatrixXd nP_vi{neighs[vi].size(), 3};
                Eigen::MatrixXd D_vi{ neighs[vi].size(), neighs[vi].size() };
                D_vi.setZero();

                for (int j = 0; j < neighs[vi].size(); j++) {
                    int vj = neighs[vi][j];
                    P_vi.row(j) = V.row(vi) - V.row(vj);
                    nP_vi.row(j) = VPrime.row(vi) - VPrime.row(vj);
                    D_vi.coeffRef(j,j) = W.coeff(vi,vj);
                }

                // We use SVD decomposition to find the resulting rotation.
                Eigen::BDCSVD<Eigen::MatrixXd> svd{ P_vi.transpose() * D_vi * nP_vi, Eigen::ComputeThinU | Eigen::ComputeThinV };
                                
                //We use the same formula as describe in the paper
                Eigen::MatrixXd U_vi = svd.matrixU();
                Eigen::MatrixXd R_vi = svd.matrixV() * U_vi.transpose();
                U_vi.col(U_vi.cols()-1) *= R_vi.determinant(); // We flip the last column of U to advert flipping an Axis while Rotating
                R_vi =svd.matrixV() * U_vi.transpose();
                rotations.push_back(R_vi); //We store them for the Global Step
            }

            //Solving the Global Step

            //We compute the resulting right hand side from Global Step
            Eigen::MatrixXd rhsGlobal{ nVertices, 3 };
            rhsGlobal.setZero();
            for (int vi = 0; vi < nVertices; vi++) {
                for (int j = 0; j < neighs[vi].size(); j++) {
                    int vj = neighs[vi][j];
                    rhsGlobal.row(vi) += W.coeff(vi,vj) * (rotations[vi] + rotations[vj]) * (V.row(vi) - V.row(vj)).transpose();
                }
                rhsGlobal.row(vi) *= 0.5;
            }

            L = (-W).eval();

            //I'm using the same technique as biharmonic to impose the fixed constraint with the Cotangant Laplacian Matrix. 
            for (int hi = 0; hi < nHandles; hi++) {
                int vi = hI.coeff(hi);
                rhsGlobal.row(vi) = hV.row(hi);
                L.row(vi) *= 0;
                L.coeffRef(vi, vi) = 1;
            }

            //We are finding the new vertices.
            Eigen::SparseLU<Eigen::SparseMatrix<double>> sLU;
            sLU.compute(L);
            VPrime = sLU.solve(rhsGlobal);
            std::cout << "Did LU converge? " << ((L * VPrime - rhsGlobal).array() < 0.01).all() << std::endl;

            Vintermediate.push_back(VPrime);
        }

        if (use_igl_solution) {
            Vintermediate.clear();
            igl::ARAPData arap_data;
            arap_data.energy = igl::ARAP_ENERGY_TYPE_SPOKES;
            igl::arap_precomputation(V, F, V.cols(), hI, arap_data);
            Eigen::MatrixXd U = V;
            igl::arap_solve(hV, arap_data, U);
            Vintermediate.push_back(V);
            Vintermediate.push_back(U);
        }

        igl::opengl::glfw::Viewer viewer;
        int i = Vintermediate.size() - 1;
        viewer.callback_key_down = [&](igl::opengl::glfw::Viewer & viewer,unsigned char key,int)->bool
        {
            switch(key)
            {
                default:
                    return false;
                case ' ':
                {
                    i = (i + 1) % Vintermediate.size();
                    viewer.data().set_mesh(Vintermediate[i], F);
                }
            }
        };
        viewer.data().set_mesh(Vintermediate[i],F);
        viewer.data().add_points(hV, Eigen::RowVector3d(1, 0, 0));
        viewer.callback_key_down(viewer,' ',0);
        viewer.data().show_lines = false;
        viewer.launch();
    }
}
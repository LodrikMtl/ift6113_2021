#include <igl/eigs.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/read_triangle_mesh.h>
#include <Eigen/Sparse>
#include <iostream>


int main(int argc, char * argv[])
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    Eigen::MatrixXd eigenvectors;
    if(!igl::read_triangle_mesh(
            argc>1?argv[1]: "../../input/horse-1.obj",V,F))
    {
        std::cout<<"failed to load mesh"<<std::endl;
    }

    eigenvectors = V.block(0,0,V.rows(), 2);

    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(V,F);
//    viewer.data().compute_normals();

    int selectedcolumn=-1; // it will be updated to 0 once we run viewer.callback_key_down(...)
    viewer.callback_key_down = [&](igl::opengl::glfw::Viewer & viewer,unsigned char key,int)->bool
    {
        switch(key)
        {
            default:
                return false;
            case ' ':
            {
                selectedcolumn = (selectedcolumn + 1) % eigenvectors.cols();
                viewer.data().set_data(eigenvectors.col(selectedcolumn));
                return true;
            }
        }
    };
    viewer.callback_key_down(viewer,' ',0);
    viewer.data().show_lines = false;
    viewer.launch();
}
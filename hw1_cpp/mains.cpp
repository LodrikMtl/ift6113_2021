#include <igl/opengl/glfw/Viewer.h>
#include <igl/read_triangle_mesh.h>
#include "trimesh.h"
#include <algorithm>
#include <filesystem>
#include <iostream>
#include <filesystem>
#include <chrono>
using namespace std;

#ifndef  PI
#define PI 3.14159265359
#endif // ! PI

double sqrt3alpha(double n) {
    return (4. - 2.*cos(2.*PI/n))/9.;
}

int main(int argc, char* argv[])
{
    //Convert argument in more useful types
    int n = stoi(string(argv[1]));
    string input_path = argv[2];

    cout << " n: " << n << " Input Path: " << input_path << "\n";

    // Modified from Stack Overflow https://stackoverflow.com/questions/35530092/c-splitting-an-absolute-file-path
    size_t botDirPos = input_path.find_last_of("/");
    size_t botExtPos = input_path.find_last_of(".");
    // get filename
    std::string filename = input_path.substr(botDirPos, botExtPos - botDirPos);

    //Read Input and Output the number of verties and faces.
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    igl::read_triangle_mesh(input_path, V, F);

    std::vector<trimesh::triangle_t> prevTriangles;

    int prevNumVertices = V.rows();
    int prevNumFaces = F.rows();

    std::cout << "Number of Vertices: " << prevNumVertices << "\n";
    std::cout << "Number of Faces: " << prevNumFaces << "\n";

    //Convert Eigen Faces to Trimesh Triangles
    prevTriangles.resize(prevNumFaces);
    for (int i = 0; i < prevNumFaces; ++i) {
        prevTriangles[i].v[0] = F(i, 0);
        prevTriangles[i].v[1] = F(i, 1);
        prevTriangles[i].v[2] = F(i, 2);
    }
    Eigen::MatrixXd prevV;

    //Build Mesh
    std::vector< trimesh::edge_t > prevEdges;
    trimesh::unordered_edges_from_triangles(prevTriangles.size(), &prevTriangles[0], prevEdges);

    trimesh::trimesh_t prevMesh;
    prevMesh.build(prevNumVertices, prevTriangles.size(), &prevTriangles[0], prevEdges.size(), &prevEdges[0]);

    trimesh::trimesh_t currMesh;
    std::vector< trimesh::edge_t > currEdges;
    std::vector< trimesh::triangle_t > currTriangles;

    // Change the configuration between butterfly and loop
    std::vector<double> stencil_new = { 1./3.,1./3.,1./3. };
    double (*stencilOld0Ring)(int n) = [](int n) {return 1 - sqrt3alpha(n); };
    double (*stencilOld1Ring)(int n) = [](int n) {return sqrt3alpha(n) / ((double)n); };
  
    //Number of times that you apply the subdivision scheme
    for (int iter = 0; iter < n; iter++) {
        cout << "Iteration: " << iter << "\n";

        // 1. Insert new Vertices in V ~ We don't need to actually insert the vertices. We just need to allocate the corresponding memory.
        Eigen::MatrixXd prevV = V; //deep copy
        int currNumVertices = prevNumVertices + prevTriangles.size();
        V.conservativeResize(currNumVertices, Eigen::NoChange);

        // 2. Split face and redo triangulation with new vertices inserted
        trimesh::triangle_t abc,abi,bci,cai;
        currTriangles.resize(prevTriangles.size() * 3);

        //For each triangles, we split triangles. In this cases, we split in 3.
        for (int fi = 0; fi < prevTriangles.size(); fi++) {
            // Start with one triangle
            abc = prevTriangles[fi];

            //We retrieve the half-edge index.
            int i = fi + prevNumVertices;

            //We split the face with the following index.
            //      c
            //     /|\
            //    / | \  
            //   /  i  \
            //  / /   \ \
            // a ------- b
            abi.v[0] = abc.v[0];
            abi.v[1] = abc.v[1];
            abi.v[2] = i;

            bci.v[0] = abc.v[1];
            bci.v[1] = abc.v[2];
            bci.v[2] = i;

            cai.v[0] = abc.v[2];
            cai.v[1] = abc.v[0];
            cai.v[2] = i;

            //We add the new triangles to the list.
            currTriangles[0 + fi * 3] = abi;
            currTriangles[1 + fi * 3] = bci;
            currTriangles[2 + fi * 3] = cai;
        }
        //We rebuild our current mesh.
        trimesh::unordered_edges_from_triangles(currTriangles.size(), &currTriangles[0], currEdges);
        currMesh.build(currNumVertices, currTriangles.size(), &currTriangles[0], currEdges.size(), &currEdges[0]);
        
        trimesh::triangle_t aji, bij;

        //2.1 We redo the triangulation to balance the valence between each vertex.
        for (int ei = 0; ei < currEdges.size(); ei++) {
            // Start with one edge
            trimesh::edge_t e = currEdges[ei];

            if (e.v[0] >= prevNumVertices || e.v[1] >= prevNumVertices) {
                continue;
            }

            trimesh::trimesh_t::halfedge_t he_ab = currMesh.halfedge(currMesh.directed_edge2he_index(e.v[0],e.v[1]));
            trimesh::trimesh_t::halfedge_t he_ba = currMesh.halfedge(currMesh.directed_edge2he_index(e.v[1], e.v[0]));

            trimesh::index_t i = currMesh.halfedge(he_ab.next_he).to_vertex;
            trimesh::index_t j = currMesh.halfedge(he_ba.next_he).to_vertex;

            //We split the face with the following index. See figure 2 into sqrt(3)-subdivision
            //      
            //    i 
            //  /   \ 
            // a --- b 
            //  \   /
            //    j 
            //   ->
            //    i 
            //  / | \ 
            // a  |  b
            //  \ | /
            //    j

            aji.v[0] = e.v[0];
            aji.v[1] = j;
            aji.v[2] = i;

            bij.v[0] = e.v[1];
            bij.v[1] = i;
            bij.v[2] = j;

            //We add the new form from the existing vertex to the list of triangles.
            currTriangles[he_ab.face] = aji;
            currTriangles[he_ba.face] = bij;
        }

        trimesh::unordered_edges_from_triangles(currTriangles.size(), &currTriangles[0], currEdges);
        currMesh.build(currNumVertices, currTriangles.size(), &currTriangles[0], currEdges.size(), &currEdges[0]);

        std::vector<trimesh::index_t> neighs;
        //3. Update Old Vertices
        for (int i = 0; i < prevNumVertices; i++) {
            prevMesh.vertex_vertex_neighbors(i, neighs);
            int n = neighs.size();
            double alpha = stencilOld0Ring(n); // Alpha represent the value affecting the 0-ring neighborhood ( current vertex)
            double beta = stencilOld1Ring(n); //  Beta represent the value affecting the 1-ring neighborhood

            //Do the actual computation to update "Old" Vertex
            V.row(i) = alpha * prevV.row(i);
            if (beta != 0)
            {
                for (auto neigh : neighs) {
                    V.row(i) += beta * prevV.row(neigh);
                }
            }
        }
        neighs.clear();


        //4. Update New Vertices

        //Stencil Orders
        //               2         
        //              / \       
        //             /   \
        //            /     \   
        //           /       \ 
        //          0 ------- 1    
        trimesh::index_t vertices_of_stencil[3];

        for (int i = 0; i < currNumVertices - prevNumVertices; i++) {
            int nvi = i + prevNumVertices;

            vertices_of_stencil[0] = prevTriangles[i].v[0];
            vertices_of_stencil[1] = prevTriangles[i].v[1];
            vertices_of_stencil[2] = prevTriangles[i].v[2];

            //We do the actual computation to update the new vertices.
            V.row(nvi) = Eigen::MatrixXd::Zero(1, 3);
            for (int j = 0; j < sizeof(vertices_of_stencil) / sizeof(vertices_of_stencil[0]); j++) {
                V.row(nvi) += stencil_new[j] * prevV.row(vertices_of_stencil[j]);
            }
        }

        prevMesh = std::move(currMesh); //Move reference without recreating the object
        prevEdges = std::move(currEdges); //Move reference without recreating the object
        prevTriangles = std::move(currTriangles); //Move reference without recreating the object
        prevNumVertices = currNumVertices;
    }
    std::cout << "\n Done \n";

    F = prevMesh.get_faces();

    std::cout << "Number of Vertices: " << V.rows() << "\n";
    std::cout << "Number of Faces: " << F.rows() << "\n";

    // Plot the mesh
    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(V, F);

    // Compute per-face normals
    Eigen::MatrixXd N_faces;
    igl::per_face_normals(V, F, N_faces);
    viewer.data().set_normals(N_faces);

    //Write mesh computed into output
    std::string output = string("../../output") + filename + "_sqrt3_" + to_string(n) + string(".obj");
std:cout << "Output to " << output;
    igl::writeOBJ(output, V, F);

    // launch viewer
    viewer.launch();

    //Clear memory
    currTriangles.clear();
    currMesh.clear();
    currEdges.clear();
    prevTriangles.clear();
    prevMesh.clear();
    prevEdges.clear();
    V.resize(0, 0);
    F.resize(0, 0);
}

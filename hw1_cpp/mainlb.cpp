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

double loopWeightedValence(double n) {
    return (64. * n) / (40. - pow(3. + 2. * std::cos((2. * PI) / n), 2)) - n;
}

int main(int argc, char* argv[])
{
    //Convert argument in more useful types
    string method_name = argv[1];
    int n = stoi(string(argv[2]));
    string input_path = argv[3];

    cout << "Methode Name: " << method_name << " n: " << n << " Input Path: " << input_path << "\n";

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
    std::vector<double> stencil_new;
    double (*stencilOld0Ring)(int n);
    double (*stencilOld1Ring)(int n);

    if (!method_name.compare("butterfly")) {
        stencil_new = { 8. / 16., -1. / 16., 2. / 16., -1. / 16., 8. / 16., -1. / 16., 2. / 16. , -1. / 16. };
        stencilOld0Ring = [](int n) {return 1.0;};
        stencilOld1Ring = [](int n) {return 0.;};
    }
    else if (!method_name.compare("loop")) {
        stencil_new = { 3./8, 0, 1./8, 0, 3./8, 0, 1./8 , 0 };
        stencilOld0Ring = [](int n) {double w_n = loopWeightedValence(n);  return w_n / (n + w_n); };
        stencilOld1Ring = [](int n) {return 1. / (n + loopWeightedValence(n)); };
    }
    else {
        std::cout << "Wrong Method Name (butterfly or loop): " << method_name;
        return 0;
    }

    //Number of times that you apply the subdivision scheme
    for (int iter = 0; iter < n; iter++) {
        cout << "Iteration: " << iter << "\n";

        // 1. Insert new Vertices in V ~ We don't need to actually insert the vertices. We just need to allocate the corresponding memory.
        Eigen::MatrixXd prevV = V; //deep copy
        int currNumVertices = prevNumVertices + prevEdges.size();
        V.conservativeResize(currNumVertices, Eigen::NoChange);
        
        // 2. Split face and redo triangulation with new vertices inserted
        trimesh::triangle_t abc, aik, bji, ckj, ijk;
        trimesh::index_t hei, hej, hek;

        currTriangles.resize(prevTriangles.size() * 4);
       
        //For each triangles, we split triangles. In this cases, we split in four.
        for (int fi = 0; fi < prevTriangles.size(); fi++) {
            // Start with one triangle
            abc = prevTriangles[fi];

            //We retrieve the half-edge index.
            hei = prevMesh.directed_edge2he_index(abc.v[0], abc.v[1]);
            hej = prevMesh.halfedge(hei).next_he;
            hek = prevMesh.halfedge(prevMesh.halfedge(hei).next_he).next_he;
            //We retrieve the index of each edge with he.../2 and we shift by the number of vertices in the previous mesh because they are as many new vertices than edges in the previous mesh
            int i = hei/2 + prevNumVertices;
            int j = hej/2 + prevNumVertices;
            int k = hek/2 + prevNumVertices;

            //We split the face with the following index.
            //     c
            //    / \
            //   k - j
            //  / \ / \
            // a --i-- b
            aik.v[0] = abc.v[0];
            aik.v[1] = i;
            aik.v[2] = k;

            bji.v[0] = abc.v[1];
            bji.v[1] = j;
            bji.v[2] = i;

            ckj.v[0] = abc.v[2];
            ckj.v[1] = k;
            ckj.v[2] = j;

            ijk.v[0] = i;
            ijk.v[1] = j;
            ijk.v[2] = k;

            //We add the new triangles to the list.
            currTriangles[0 + fi * 4] = aik;
            currTriangles[1 + fi * 4] = bji;
            currTriangles[2 + fi * 4] = ckj;
            currTriangles[3 + fi * 4] = ijk;
        }
        prevTriangles.clear();

        //We rebuild our current mesh.
        trimesh::unordered_edges_from_triangles(currTriangles.size(), &currTriangles[0], currEdges);
        currMesh.build(currNumVertices, currTriangles.size(), &currTriangles[0], currEdges.size(), &currEdges[0]);

        std::vector<trimesh::index_t> neighs;
        //3. Update Old Vertices
        for (int i = 0; i < prevNumVertices; i++) {
            prevMesh.vertex_vertex_neighbors(i,neighs);
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
        //     7        6         5
        //      \       / \       /
        //       \14   /  9\10 11/
        //        \ 13/12   \   /
        //         \ /   8   \ /
        //          0 ------- 4
        //         / \   1   / \
        //        /4 3\2   5/6 7\
        //       /     \   /     \
        //      /       \ /       \
        //     1         2          3        
        trimesh::index_t vertices_of_stencil[8];
        trimesh::trimesh_t::halfedge_t he_1, he_2, he_5, he_8, he_9, he_12;
        
        for (int i = 0; i < prevEdges.size(); i++) {
            int nvi = i + prevNumVertices;

            trimesh::edge_t edge = prevEdges[i];
            int hei = prevMesh.directed_edge2he_index(edge.v[0], edge.v[1]); //Half-Edge Index
            // Retrieve the following half-edge he_x where x represent their "ID" in the stencil orders schema
            he_1 = prevMesh.halfedge(hei); 
            he_2 = prevMesh.halfedge(he_1.next_he);
            he_5 = prevMesh.halfedge(he_2.next_he);
            he_8 = prevMesh.halfedge(he_1.opposite_he);
            he_9 = prevMesh.halfedge(he_8.next_he);
            he_12 = prevMesh.halfedge(he_9.next_he);

            //We find the vertices for the "bottom" part of the stencil 
            vertices_of_stencil[0] = he_1.to_vertex;
            vertices_of_stencil[1] = prevMesh.halfedge(prevMesh.halfedge(he_2.opposite_he).next_he).to_vertex;
            vertices_of_stencil[2] = he_2.to_vertex;
            vertices_of_stencil[3] = prevMesh.halfedge(prevMesh.halfedge(he_5.opposite_he).next_he).to_vertex;

            //We find the vertices for the "upper" part of the stencil. We repeat the same operation,but we start at he_8.
            vertices_of_stencil[4] = he_8.to_vertex;
            vertices_of_stencil[5] = prevMesh.halfedge(prevMesh.halfedge(he_9.opposite_he).next_he).to_vertex;
            vertices_of_stencil[6] = he_9.to_vertex;
            vertices_of_stencil[7] = prevMesh.halfedge(prevMesh.halfedge(he_12.opposite_he).next_he).to_vertex;
            
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
    std::string output = string("../../output") + filename + "_" + method_name + "_" + to_string(n) + string(".obj");
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

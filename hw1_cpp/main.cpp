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

int retrieveIndexForEdge(trimesh::edge_t e, std::vector<trimesh::edge_t> edges) {
    return std::distance(edges.begin(), std::find_if(edges.begin(), edges.end(), [&](const auto& val) {return (val.v[0] == e.v[0] && val.v[1] == e.v[1]) || (val.v[1] == e.v[0] && val.v[0] == e.v[1]); }));
}

double loopWeightedValence(double n) {
    return (64. * n) / (40. - pow(3. + 2. * std::cos((2. * PI) / n), 2)) - n;
}

int main(int argc, char* argv[])
{
    string method_name = argv[1];
    int n = stoi(string(argv[2]));
    string input_path = argv[3];

    cout << "Methode Name: " << method_name << " n: " << n << " Input Path: " << input_path << "\n";

    // Modified from Stack Overflow https://stackoverflow.com/questions/35530092/c-splitting-an-absolute-file-path
    size_t botDirPos = input_path.find_last_of("/");
    size_t botExtPos = input_path.find_last_of(".");
    // get file
    std::string filename = input_path.substr(botDirPos, botExtPos - botDirPos);

    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    igl::read_triangle_mesh(input_path, V, F);

    std::vector<trimesh::triangle_t> prevTriangles;

    int prevNumVertices = V.rows();
    int prevNumFaces = F.rows();

    std::cout << "Number of Vertices: " << prevNumVertices << "\n";
    std::cout << "Number of Faces: " << prevNumFaces << "\n";

    prevTriangles.resize(prevNumFaces);
    for (int i = 0; i < prevNumFaces; ++i) {
        prevTriangles[i].v[0] = F(i, 0);
        prevTriangles[i].v[1] = F(i, 1);
        prevTriangles[i].v[2] = F(i, 2);
    }
    Eigen::MatrixXd prevV;

    std::vector< trimesh::edge_t > prevEdges;
    trimesh::unordered_edges_from_triangles(prevTriangles.size(), &prevTriangles[0], prevEdges);

    trimesh::trimesh_t prevMesh;
    prevMesh.build(prevNumVertices, prevTriangles.size(), &prevTriangles[0], prevEdges.size(), &prevEdges[0]);

    trimesh::trimesh_t currMesh;
    std::vector< trimesh::edge_t > currEdges;
    std::vector< trimesh::triangle_t > currTriangles;

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

    for (int iter = 0; iter < n; iter++) {
        cout << "Iteration: " << iter << "\n";

        currMesh.clear();
        currEdges.clear();
        currTriangles.clear();

        // 1. Insert new Vertices at midpoint of each edge.
        Eigen::MatrixXd prevV = V;
        int currNumVertices = prevNumVertices + prevEdges.size();
        V.conservativeResize(currNumVertices, Eigen::NoChange);

        clock_t time = 0;

        // 2. Remove Old Connections and Remake new
        for (int fi = prevTriangles.size() - 1; fi >= 0; fi--) {
            // Start with one triangle
            trimesh::triangle_t abc = prevTriangles[fi];

            trimesh::edge_t ei, ej, ek;
            ei.v[0] = abc.v[0];
            ei.v[1] = abc.v[1];
            ej.v[0] = abc.v[1];
            ej.v[1] = abc.v[2];
            ek.v[0] = abc.v[2];
            ek.v[1] = abc.v[0];

            //TODO: Improve this part -> This is very expensive. Must do a better search!
            //prevMesh.directed_edge2he_index(ei.v[0], ei.v[1]);
            clock_t start = clock();
            int i = retrieveIndexForEdge(ei, prevEdges) + prevNumVertices;
            int j = retrieveIndexForEdge(ej, prevEdges) + prevNumVertices;
            int k = retrieveIndexForEdge(ek, prevEdges) + prevNumVertices;
            clock_t end = clock();
            time += end - start;

            trimesh::triangle_t aik, bji, ckj, ijk;
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

            currTriangles.push_back(aik);
            currTriangles.push_back(bji);
            currTriangles.push_back(ckj);
            currTriangles.push_back(ijk);
        }
        trimesh::unordered_edges_from_triangles(currTriangles.size(), &currTriangles[0], currEdges);
        currMesh.build(currNumVertices, currTriangles.size(), &currTriangles[0], currEdges.size(), &currEdges[0]);


        //3. Update Old Vertices
        for (int i = 0; i < prevNumVertices; i++) {
            std::vector<trimesh::index_t> neighs = prevMesh.vertex_vertex_neighbors(i);
            int n = neighs.size();
            double alpha = stencilOld0Ring(n);
            double beta = stencilOld1Ring(n);
            V.row(i) = alpha * prevV.row(i);
            if (beta != 0)
            {
                for (auto neigh : neighs) {
                    V.row(i) += beta * prevV.row(neigh);
                }
            }
        }

        //4. Update New Vertices

        //Stencil order values
        //     8         7         6
        //      \       / \       /
        //       \14   /  9\10 11/
        //        \ 13/12   \   /
        //         \ /   8   \ /
        //          1 ------- 5
        //         / \   1   / \
        //        /4 3\2   5/6 7\
        //       /     \   /     \
        //      /       \ /       \
        //     2         3          4        
        trimesh::index_t vertices_of_stencil[8];

        for (int i = 0; i < prevEdges.size(); i++) {
            int nvi = i + prevNumVertices;

            trimesh::edge_t edge = prevEdges[i];
            int hei = prevMesh.directed_edge2he_index(edge.v[0], edge.v[1]);
            trimesh::trimesh_t::halfedge_t he_1 = prevMesh.halfedge(hei);
            trimesh::trimesh_t::halfedge_t he_2 = prevMesh.halfedge(he_1.next_he);
            trimesh::trimesh_t::halfedge_t he_5 = prevMesh.halfedge(he_2.next_he);
            trimesh::trimesh_t::halfedge_t he_8 = prevMesh.halfedge(he_1.opposite_he);
            trimesh::trimesh_t::halfedge_t he_9 = prevMesh.halfedge(he_8.next_he);
            trimesh::trimesh_t::halfedge_t he_12 = prevMesh.halfedge(he_9.next_he);

            vertices_of_stencil[0] = he_1.to_vertex;
            vertices_of_stencil[1] = prevMesh.halfedge(prevMesh.halfedge(he_2.opposite_he).next_he).to_vertex;
            vertices_of_stencil[2] = he_2.to_vertex;
            vertices_of_stencil[3] = prevMesh.halfedge(prevMesh.halfedge(he_5.opposite_he).next_he).to_vertex;

            vertices_of_stencil[4] = he_8.to_vertex;
            vertices_of_stencil[5] = prevMesh.halfedge(prevMesh.halfedge(he_9.opposite_he).next_he).to_vertex;
            vertices_of_stencil[6] = he_9.to_vertex;
            vertices_of_stencil[7] = prevMesh.halfedge(prevMesh.halfedge(he_12.opposite_he).next_he).to_vertex;
            
            V.row(nvi) = Eigen::MatrixXd::Zero(1, 3);
            for (int j = 0; j < sizeof(vertices_of_stencil) / sizeof(vertices_of_stencil[0]); j++) {
                V.row(nvi) += stencil_new[j] * prevV.row(vertices_of_stencil[j]);
            }
        }

        prevMesh = currMesh;
        prevEdges = currEdges;
        prevTriangles = currTriangles;
        prevNumVertices = currNumVertices;

        std::cout << ((float)time) / CLOCKS_PER_SEC << "\n";
    }
    std::cout << "\n Done \n";

    F = currMesh.get_faces();

    std::cout << "Number of Vertices: " << V.rows() << "\n";
    std::cout << "Number of Faces: " << F.rows() << "\n";

    // Plot the mesh
    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(V, F);

    // Compute per-face normals
    Eigen::MatrixXd N_faces;
    igl::per_face_normals(V, F, N_faces);
    viewer.data().set_normals(N_faces);
    viewer.data().point_size = 10;
    viewer.data().set_points(prevV, Eigen::RowVector3d(1, 0, 0));
    std::string output = string("../../output/") + filename + "_" + method_name + "_" + to_string(n) + string(".obj");
    std:cout << "Output to " << output;
    igl::writeOBJ(output, V, F);

    // launch viewer
    viewer.launch();

    currTriangles.clear();
    currMesh.clear();
    currEdges.clear();
    prevTriangles.clear();
    prevMesh.clear();
    prevEdges.clear();
    V.resize(0, 0);
    F.resize(0, 0);
}

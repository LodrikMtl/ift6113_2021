/*
#include <igl/opengl/glfw/Viewer.h>
#include <igl/read_triangle_mesh.h>
#include "trimesh.h"
#include <algorithm>

int retrieveIndexForEdge(trimesh::edge_t e, std::vector<trimesh::edge_t> edges) {
    return std::distance(edges.begin(), std::find_if(edges.begin(), edges.end(), [&](const auto& val) {return (val.v[0] == e.v[0] && val.v[1] == e.v[1]) || (val.v[1] == e.v[0] && val.v[0] == e.v[1]); }));
}

int main(int argc, char *argv[])
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    igl::read_triangle_mesh("../../output/cube_butterfly_3.obj", V,F);

    std::vector< trimesh::triangle_t > prevTriangles;

    int prevNumVertices = V.rows();
    int prevNumFaces = F.rows();

    std::cout << "Number of Vertices: " << prevNumVertices << "\n";
    std::cout << "Number of Faces: " << prevNumFaces << "\n";

    prevTriangles.resize( prevNumFaces );
    for (int i=0; i<prevNumFaces; ++i){
        prevTriangles[i].v[0] = F(i,0);
        prevTriangles[i].v[1] = F(i,1);
        prevTriangles[i].v[2] = F(i,2);
    }
    std::vector< trimesh::edge_t > prevEdges;
    trimesh::unordered_edges_from_triangles(prevTriangles.size(), &prevTriangles[0], prevEdges );

    trimesh::trimesh_t prevMesh;
    prevMesh.build( prevNumVertices, prevTriangles.size(), &prevTriangles[0], prevEdges.size(), &prevEdges[0] );

    // 1. Insert new Vertices at midpoint of each edge.
    int numCurrVertices = prevNumVertices + prevEdges.size();
    V.conservativeResize(numCurrVertices,Eigen::NoChange);
    for (int nvi = prevNumVertices; nvi < numCurrVertices; nvi++)
    {
        V.row(nvi) = 0.5 * V.row(prevEdges[nvi - prevNumVertices].start()) + 0.5 * V.row(prevEdges[nvi - prevNumVertices].end());
    }

    // 2. Remove Old Connections and Remake new
    std::vector< trimesh::triangle_t > currTriangles;
    for (int fi = prevNumFaces - 1; fi >= 0; fi--) {
        trimesh::triangle_t abc = prevTriangles[fi];

        trimesh::edge_t ei, ej, ek;
        ei.v[0] = abc.v[0];
        ei.v[1] = abc.v[1];
        ej.v[0] = abc.v[1];
        ej.v[1] = abc.v[2];
        ek.v[0] = abc.v[2];
        ek.v[1] = abc.v[0];

        int i = retrieveIndexForEdge(ei, prevEdges) + prevNumVertices;
        int j = retrieveIndexForEdge(ej, prevEdges) + prevNumVertices;
        int k = retrieveIndexForEdge(ek, prevEdges) + prevNumVertices;

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
    std::vector< trimesh::edge_t > currEdges;
    trimesh::unordered_edges_from_triangles(currTriangles.size(), &currTriangles[0], currEdges);

    trimesh::trimesh_t currMesh;
    currMesh.build(numCurrVertices, currTriangles.size(), &currTriangles[0], currEdges.size(), &currEdges[0]);

    //3. Update Old Vertices
    // Nothing for butterfly
   
    //4. Update New Vertices
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

        double stencil[] = { 8, -1, 2, -1, 8, -1, 2 , -1 };
        trimesh::index_t vertices_of_stencil[8];

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
            V.row(nvi) += stencil[j] / 16 * V.row(vertices_of_stencil[j]);
        }
    }
    F = currMesh.get_faces();

    std::cout << "Number of Vertices: " << V.rows() << "\n";
    std::cout << "Number of Faces: " << F.rows() << "\n";


    // Plot the mesh
    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(V, F);
    //viewer.data().add_points(V, Eigen::RowVector3d(1, 0, 0));

    // Compute per-face normals
    Eigen::MatrixXd N_faces;
    igl::per_face_normals(V, F, N_faces);
    viewer.data().set_normals(N_faces);

    igl::writeOBJ("../../output/cube_butterfly_4.obj", V, F);

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

*/
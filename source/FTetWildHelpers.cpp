#include "comp_mesh_core/FTetWildHelpers.h"

#include "floattetwild/FloatTetwild.h"
#include "mesh/algorithms/IO.h"

#include <bitset>

namespace composite_mesher
{
mesh::Mesh composite_mesher::tetrahedralize_ftetwild(const mesh::Mesh& object_surface,
                                                     const double ideal_relative_edge_length)
{
    return tetrahedralize_ftetwild(std::vector<mesh::Mesh>({ object_surface }),
                                   ideal_relative_edge_length);
}

mesh::Mesh tetrahedralize_ftetwild(const std::vector<mesh::Mesh>& objects_surfaces,
                                   const double ideal_relative_edge_length)
{
    std::cout << "Starting fTetWild tetrahedralization." << std::endl << std::endl;

    std::vector<Eigen::MatrixXd> vertices_by_object(objects_surfaces.size());
    for (auto o : utility::range(objects_surfaces.size()))
    {
        auto& vertices = vertices_by_object[o];
        vertices.resize(objects_surfaces[o].nodes().size(), 3);
        for (int i = 0; i < vertices.rows(); ++i)
        {
            vertices(i, 0) = objects_surfaces[o].nodes()[i].x();
            vertices(i, 1) = objects_surfaces[o].nodes()[i].y();
            vertices(i, 2) = objects_surfaces[o].nodes()[i].z();
        }
    }

    std::vector<Eigen::MatrixXi> triangles_by_object(objects_surfaces.size());
    for (auto i : utility::range(objects_surfaces.size()))
    {
        auto& triangles = triangles_by_object[i];
        triangles.resize(objects_surfaces[i].elements().size(), 3);
        for (auto& e : objects_surfaces[i].elements())
        {
            for (auto j : utility::range(3))
            {
                triangles(e.local_index(), j) = e.connectivity()[j];
            }
        }
    }

    floatTetWild::Parameters params;
    params.ideal_edge_length_rel = ideal_relative_edge_length;
    params.eps_rel               = 2e-2 * params.ideal_edge_length_rel;
    params.max_its               = 4;
    params.use_input_for_wn      = true;
    params.is_quiet              = true;
    params.num_threads           = std::thread::hardware_concurrency();

    Eigen::MatrixXd boundary_vertices;
    Eigen::MatrixXi boundary_triangles;

    Eigen::MatrixXd verts_out;
    Eigen::MatrixXi tets_out;
    std::vector<int> attributes;
    floatTetWild::tetrahedralization(vertices_by_object,
                                     triangles_by_object,
                                     boundary_vertices,
                                     boundary_triangles,
                                     params,
                                     verts_out,
                                     tets_out,
                                     attributes);

    //////////////////////////////////////////////////////////////////////////////////////
    // Extract volume information into our mesh data structure

    mesh::Mesh tetrahedral_mesh;
    for (auto row : utility::range(verts_out.rows()))
    {
        tetrahedral_mesh.nodes().add({ verts_out(row, 0), verts_out(row, 1), verts_out(row, 2) });
    }
    for (auto row : utility::range(tets_out.rows()))
    {
        std::vector<mesh::LocalIndex> conn = { static_cast<mesh::LocalIndex>(tets_out(row, 0)),
                                               static_cast<mesh::LocalIndex>(tets_out(row, 1)),
                                               static_cast<mesh::LocalIndex>(tets_out(row, 2)),
                                               static_cast<mesh::LocalIndex>(tets_out(row, 3)) };
        tetrahedral_mesh.elements().add(mesh::ElementType::T3D4N, conn);
    }

    //////////////////////////////////////////////////////////////////////////////////////
    // Assign material IDs by the lowest input ID a tetrahedra lies within, which works out
    // since the matrix is always the highest
    std::vector<unsigned int> material_ids(attributes.size(),
                                           std::numeric_limits<unsigned int>::max());
    for (auto i : utility::range(attributes.size()))
    {
        using RegionBitset = std::bitset<std::numeric_limits<unsigned long long>::digits>;
        RegionBitset bitset(attributes[i]);
        for (auto surface_index : utility::range(bitset.size()))
        {
            if (bitset.test(surface_index))
            {
                material_ids[i] = surface_index;
                break;
            }
        }
        tetrahedral_mesh.elements()[i].component_id() = material_ids[i];
    }

    return tetrahedral_mesh;
}
}
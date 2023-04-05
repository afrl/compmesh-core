#include "composite_mesher/MarchingTetrahedra.h"

#include "simple_math/Utility.h"

namespace composite_mesher
{

#pragma region helpers

namespace detail
{
enum class VertexType : int8_t
{
    EDGE_INTERSECTION,
    EXISTING_VERTEX,
    NA
};

struct LookupEntry
{
    // Edge index of the tet that is intersected (-1 indicates null value)
    std::array<int, 4> intersected_edges;
    // Indicates if the vertex for each intersection is inside the edge or a vertex is the mesh
    std::array<VertexType, 4> vertex_types;
};

// Note: 0b11110001 = (int) 241, but not all cases are possible. 
using LookupTable = std::array<LookupEntry, 0b11110001>;

inline auto get_cell_lookup_table_kernel()
{
    // Table format: { {edge index of the tet that is intersected}, {corresponding vertex type} }
    LookupTable table;
    for (auto& t : table)
    {
        t = { { -1, -1, -1, -1 },
              { VertexType::NA, VertexType::NA, VertexType::NA, VertexType::NA } };
    }
    // Standard cases
    table[0b0001] = { { 0, 2, 1, -1 },
                      { VertexType::EDGE_INTERSECTION,
                        VertexType::EDGE_INTERSECTION,
                        VertexType::EDGE_INTERSECTION,
                        VertexType::NA } };
    table[0b0010] = { { 0, 3, 4, -1 },
                      { VertexType::EDGE_INTERSECTION,
                        VertexType::EDGE_INTERSECTION,
                        VertexType::EDGE_INTERSECTION,
                        VertexType::NA } };
    table[0b0011] = { { 2, 1, 3, 4 },
                      { VertexType::EDGE_INTERSECTION,
                        VertexType::EDGE_INTERSECTION,
                        VertexType::EDGE_INTERSECTION,
                        VertexType::EDGE_INTERSECTION } };
    table[0b0100] = { { 5, 3, 1, -1 },
                      { VertexType::EDGE_INTERSECTION,
                        VertexType::EDGE_INTERSECTION,
                        VertexType::EDGE_INTERSECTION,
                        VertexType::NA } };
    table[0b0101] = { { 0, 2, 5, 3 },
                      { VertexType::EDGE_INTERSECTION,
                        VertexType::EDGE_INTERSECTION,
                        VertexType::EDGE_INTERSECTION,
                        VertexType::EDGE_INTERSECTION } };
    table[0b0110] = { { 0, 1, 5, 4 },
                      { VertexType::EDGE_INTERSECTION,
                        VertexType::EDGE_INTERSECTION,
                        VertexType::EDGE_INTERSECTION,
                        VertexType::EDGE_INTERSECTION } };
    table[0b0111] = { { 2, 5, 4, -1 },
                      { VertexType::EDGE_INTERSECTION,
                        VertexType::EDGE_INTERSECTION,
                        VertexType::EDGE_INTERSECTION,
                        VertexType::NA } };
    table[0b1000] = { { 4, 5, 2, -1 },
                      { VertexType::EDGE_INTERSECTION,
                        VertexType::EDGE_INTERSECTION,
                        VertexType::EDGE_INTERSECTION,
                        VertexType::NA } };
    table[0b1001] = { { 0, 4, 5, 1 },
                      { VertexType::EDGE_INTERSECTION,
                        VertexType::EDGE_INTERSECTION,
                        VertexType::EDGE_INTERSECTION,
                        VertexType::EDGE_INTERSECTION } };
    table[0b1010] = { { 0, 3, 5, 2 },
                      { VertexType::EDGE_INTERSECTION,
                        VertexType::EDGE_INTERSECTION,
                        VertexType::EDGE_INTERSECTION,
                        VertexType::EDGE_INTERSECTION } };
    table[0b1011] = { { 1, 3, 5, -1 },
                      { VertexType::EDGE_INTERSECTION,
                        VertexType::EDGE_INTERSECTION,
                        VertexType::EDGE_INTERSECTION,
                        VertexType::NA } };
    table[0b1100] = { { 4, 3, 1, 2 },
                      { VertexType::EDGE_INTERSECTION,
                        VertexType::EDGE_INTERSECTION,
                        VertexType::EDGE_INTERSECTION,
                        VertexType::EDGE_INTERSECTION } };
    table[0b1101] = { { 0, 4, 3, -1 },
                      { VertexType::EDGE_INTERSECTION,
                        VertexType::EDGE_INTERSECTION,
                        VertexType::EDGE_INTERSECTION,
                        VertexType::NA } };
    table[0b1110] = { { 0, 1, 2, -1 },
                      { VertexType::EDGE_INTERSECTION,
                        VertexType::EDGE_INTERSECTION,
                        VertexType::EDGE_INTERSECTION,
                        VertexType::NA } };

    // Variations of case 1
    table[0b00100001] = { { 1, 2, 1, -1 },
                          { VertexType::EXISTING_VERTEX,
                            VertexType::EDGE_INTERSECTION,
                            VertexType::EDGE_INTERSECTION,
                            VertexType::NA } };
    table[0b01000001] = { { 0, 2, 2, -1 },
                          { VertexType::EDGE_INTERSECTION,
                            VertexType::EDGE_INTERSECTION,
                            VertexType::EXISTING_VERTEX,
                            VertexType::NA } };
    table[0b10000001] = { { 0, 3, 1, -1 },
                          { VertexType::EDGE_INTERSECTION,
                            VertexType::EXISTING_VERTEX,
                            VertexType::EDGE_INTERSECTION,
                            VertexType::NA } };
    table[0b01100001] = { { 1, 2, 2, -1 },
                          { VertexType::EXISTING_VERTEX,
                            VertexType::EDGE_INTERSECTION,
                            VertexType::EXISTING_VERTEX,
                            VertexType::NA } };
    table[0b10100001] = { { 1, 3, 1, -1 },
                          { VertexType::EXISTING_VERTEX,
                            VertexType::EXISTING_VERTEX,
                            VertexType::EDGE_INTERSECTION,
                            VertexType::NA } };
    table[0b11000001] = { { 0, 3, 2, -1 },
                          { VertexType::EDGE_INTERSECTION,
                            VertexType::EXISTING_VERTEX,
                            VertexType::EXISTING_VERTEX,
                            VertexType::NA } };
    table[0b11100001] = { { 1, 3, 2, -1 },
                          { VertexType::EXISTING_VERTEX,
                            VertexType::EXISTING_VERTEX,
                            VertexType::EXISTING_VERTEX,
                            VertexType::NA } };

    // Variations of case 2
    table[0b00010010] = { { 0, 3, 4, -1 },
                          { VertexType::EXISTING_VERTEX,
                            VertexType::EDGE_INTERSECTION,
                            VertexType::EDGE_INTERSECTION,
                            VertexType::NA } };
    table[0b01000010] = { { 0, 2, 4, -1 },
                          { VertexType::EDGE_INTERSECTION,
                            VertexType::EXISTING_VERTEX,
                            VertexType::EDGE_INTERSECTION,
                            VertexType::NA } };
    table[0b10000010] = { { 0, 3, 3, -1 },
                          { VertexType::EDGE_INTERSECTION,
                            VertexType::EDGE_INTERSECTION,
                            VertexType::EXISTING_VERTEX,
                            VertexType::NA } };
    table[0b01010010] = { { 0, 2, 4, -1 },
                          { VertexType::EXISTING_VERTEX,
                            VertexType::EXISTING_VERTEX,
                            VertexType::EDGE_INTERSECTION,
                            VertexType::NA } };
    table[0b11000010] = { { 0, 2, 3, -1 },
                          { VertexType::EDGE_INTERSECTION,
                            VertexType::EXISTING_VERTEX,
                            VertexType::EXISTING_VERTEX,
                            VertexType::NA } };
    table[0b10010010] = { { 0, 3, 3, -1 },
                          { VertexType::EXISTING_VERTEX,
                            VertexType::EDGE_INTERSECTION,
                            VertexType::EXISTING_VERTEX,
                            VertexType::NA } };
    table[0b11010010] = { { 0, 2, 3, -1 },
                          { VertexType::EXISTING_VERTEX,
                            VertexType::EXISTING_VERTEX,
                            VertexType::EXISTING_VERTEX,
                            VertexType::NA } };

    // Variations of case 3
    table[0b01000011] = { { 2, 2, 4, -1 },
                          { VertexType::EDGE_INTERSECTION,
                            VertexType::EXISTING_VERTEX,
                            VertexType::EDGE_INTERSECTION,
                            VertexType::NA } };
    table[0b10000011] = { { 3, 1, 3, -1 },
                          { VertexType::EXISTING_VERTEX,
                            VertexType::EDGE_INTERSECTION,
                            VertexType::EDGE_INTERSECTION,
                            VertexType::NA } };

    // Variations of case 4
    table[0b00010100] = { { 5, 3, 0, -1 },
                          { VertexType::EDGE_INTERSECTION,
                            VertexType::EDGE_INTERSECTION,
                            VertexType::EXISTING_VERTEX,
                            VertexType::NA } };
    table[0b00100100] = { { 5, 1, 1, -1 },
                          { VertexType::EDGE_INTERSECTION,
                            VertexType::EXISTING_VERTEX,
                            VertexType::EDGE_INTERSECTION,
                            VertexType::NA } };
    table[0b10000100] = { { 3, 3, 1, -1 },
                          { VertexType::EXISTING_VERTEX,
                            VertexType::EDGE_INTERSECTION,
                            VertexType::EDGE_INTERSECTION,
                            VertexType::NA } };
    table[0b00110100] = { { 5, 1, 0, -1 },
                          { VertexType::EDGE_INTERSECTION,
                            VertexType::EXISTING_VERTEX,
                            VertexType::EXISTING_VERTEX,
                            VertexType::NA } };
    table[0b10100100] = { { 3, 1, 1, -1 },
                          { VertexType::EXISTING_VERTEX,
                            VertexType::EXISTING_VERTEX,
                            VertexType::EDGE_INTERSECTION,
                            VertexType::NA } };
    table[0b10010100] = { { 3, 3, 0, -1 },
                          { VertexType::EXISTING_VERTEX,
                            VertexType::EDGE_INTERSECTION,
                            VertexType::EXISTING_VERTEX,
                            VertexType::NA } };
    table[0b10110100] = { { 3, 1, 0, -1 },
                          { VertexType::EXISTING_VERTEX,
                            VertexType::EXISTING_VERTEX,
                            VertexType::EXISTING_VERTEX,
                            VertexType::NA } };

    // Variations of case 5
    table[0b00100101] = { { 1, 2, 5, -1 },
                          { VertexType::EXISTING_VERTEX,
                            VertexType::EDGE_INTERSECTION,
                            VertexType::EDGE_INTERSECTION,
                            VertexType::NA } };
    table[0b10000101] = { { 0, 3, 3, -1 },
                          { VertexType::EDGE_INTERSECTION,
                            VertexType::EXISTING_VERTEX,
                            VertexType::EDGE_INTERSECTION,
                            VertexType::NA } };

    // Variations of case 6
    table[0b00010110] = { { 0, 5, 4, -1 },
                          { VertexType::EXISTING_VERTEX,
                            VertexType::EDGE_INTERSECTION,
                            VertexType::EDGE_INTERSECTION,
                            VertexType::NA } };
    table[0b10000110] = { { 0, 1, 3, -1 },
                          { VertexType::EDGE_INTERSECTION,
                            VertexType::EDGE_INTERSECTION,
                            VertexType::EXISTING_VERTEX,
                            VertexType::NA } };

    // No variations for case 7

    // Variations of case 8
    table[0b00011000] = { { 4, 5, 0, -1 },
                          { VertexType::EDGE_INTERSECTION,
                            VertexType::EDGE_INTERSECTION,
                            VertexType::EXISTING_VERTEX,
                            VertexType::NA } };
    table[0b00101000] = { { 1, 5, 2, -1 },
                          { VertexType::EXISTING_VERTEX,
                            VertexType::EDGE_INTERSECTION,
                            VertexType::EDGE_INTERSECTION,
                            VertexType::NA } };
    table[0b01001000] = { { 4, 2, 2, -1 },
                          { VertexType::EDGE_INTERSECTION,
                            VertexType::EXISTING_VERTEX,
                            VertexType::EDGE_INTERSECTION,
                            VertexType::NA } };
    table[0b00111000] = { { 1, 5, 0, -1 },
                          { VertexType::EXISTING_VERTEX,
                            VertexType::EDGE_INTERSECTION,
                            VertexType::EXISTING_VERTEX,
                            VertexType::NA } };
    table[0b01011000] = { { 4, 2, 0, -1 },
                          { VertexType::EDGE_INTERSECTION,
                            VertexType::EXISTING_VERTEX,
                            VertexType::EXISTING_VERTEX,
                            VertexType::NA } };
    table[0b01101000] = { { 1, 2, 2, -1 },
                          { VertexType::EXISTING_VERTEX,
                            VertexType::EXISTING_VERTEX,
                            VertexType::EDGE_INTERSECTION,
                            VertexType::NA } };
    table[0b01111000] = { { 1, 2, 0, -1 },
                          { VertexType::EXISTING_VERTEX,
                            VertexType::EXISTING_VERTEX,
                            VertexType::EXISTING_VERTEX,
                            VertexType::NA } };

    // Variations of case 9
    table[0b00101001] = { { 1, 5, 1, -1 },
                          { VertexType::EXISTING_VERTEX,
                            VertexType::EDGE_INTERSECTION,
                            VertexType::EDGE_INTERSECTION,
                            VertexType::NA } };
    table[0b01001001] = { { 0, 4, 2, -1 },
                          { VertexType::EDGE_INTERSECTION,
                            VertexType::EDGE_INTERSECTION,
                            VertexType::EXISTING_VERTEX,
                            VertexType::NA } };

    // Variations of case 10
    table[0b00011010] = { { 0, 3, 5, -1 },
                          { VertexType::EXISTING_VERTEX,
                            VertexType::EDGE_INTERSECTION,
                            VertexType::EDGE_INTERSECTION,
                            VertexType::NA } };
    table[0b01001010] = { { 0, 2, 2, -1 },
                          { VertexType::EDGE_INTERSECTION,
                            VertexType::EXISTING_VERTEX,
                            VertexType::EDGE_INTERSECTION,
                            VertexType::NA } };

    // No variations for case 11

    // Variations of case 12
    table[0b00011100] = { { 4, 3, 0, -1 },
                          { VertexType::EDGE_INTERSECTION,
                            VertexType::EDGE_INTERSECTION,
                            VertexType::EXISTING_VERTEX,
                            VertexType::NA } };
    table[0b00101100] = { { 1, 1, 2, -1 },
                          { VertexType::EXISTING_VERTEX,
                            VertexType::EDGE_INTERSECTION,
                            VertexType::EDGE_INTERSECTION,
                            VertexType::NA } };

    // No variations for case 13

    // No variations for case 14

    // New cases that are not variations
    table[0b01110000] = { { 2, 1, 0, -1 },
                          { VertexType::EXISTING_VERTEX,
                            VertexType::EXISTING_VERTEX,
                            VertexType::EXISTING_VERTEX,
                            VertexType::NA } };
    table[0b11010000] = { { 2, 0, 3, -1 },
                          { VertexType::EXISTING_VERTEX,
                            VertexType::EXISTING_VERTEX,
                            VertexType::EXISTING_VERTEX,
                            VertexType::NA } };
    table[0b10110000] = { { 1, 3, 0, -1 },
                          { VertexType::EXISTING_VERTEX,
                            VertexType::EXISTING_VERTEX,
                            VertexType::EXISTING_VERTEX,
                            VertexType::NA } };
    table[0b11100000] = { { 3, 1, 2, -1 },
                          { VertexType::EXISTING_VERTEX,
                            VertexType::EXISTING_VERTEX,
                            VertexType::EXISTING_VERTEX,
                            VertexType::NA } };

    return table;
}

const auto& get_cell_lookup_table()
{
    static auto table = get_cell_lookup_table_kernel();
    return table;
}

const auto& get_edge_lookup_table()
{
    static std::array<std::array<int, 2>, 6> table = {
        { { 0, 1 }, { 0, 2 }, { 0, 3 }, { 1, 2 }, { 1, 3 }, { 2, 3 } }
    };
    return table;
}

/// We're hashing edges to deduplicate using 64 bit ints. The upper and lower
/// 32 bits of a key are the indices of vertices in the mesh. The implication is
/// that you can only have 2^32 vertices which I have deemed sufficient for
/// anything reasonable.
int64_t make_edge_key(const std::pair<int32_t, int32_t>& p)
{
    std::int64_t ret = 0;
    ret |= p.first;
    ret |= static_cast<std::int64_t>(p.second) << 32;
    return ret;
};
}

#pragma endregion helpers

void extended_marching_tets(const Eigen::MatrixXd& tet_vertices,
                            const Eigen::MatrixXi& tet_connectivity,
                            const Eigen::VectorXd& isovalues,
                            double surface_isovalue,
                            Eigen::MatrixXd& tri_vertices,
                            Eigen::MatrixXi& tri_connectivity,
                            Eigen::VectorXi& tri_vertex_to_tet_vertex_map,
                            Eigen::MatrixXi& tri_vertex_to_tet_split_edge_map)
{
    MarchingTetState state;
    extended_marching_tets(tet_vertices, tet_connectivity, isovalues, surface_isovalue, state);
    recover_isosurface(tet_vertices,
                       state,
                       tri_vertices,
                       tri_connectivity,
                       tri_vertex_to_tet_vertex_map,
                       tri_vertex_to_tet_split_edge_map);
}

void extended_marching_tets(const Eigen::MatrixXd& tet_vertices,
                            const Eigen::MatrixXi& tet_connectivity,
                            const Eigen::VectorXd& isovalues,
                            const double& surface_isovalue,
                            MarchingTetState& state)
{
    assert(tet_vertices.cols() == 3 && tet_vertices.rows() >= 4);
    assert(tet_connectivity.cols() == 4 && tet_connectivity.rows() >= 1);
    assert(isovalues.cols() == 1);


    const auto& cell_lookup_table = detail::get_cell_lookup_table();
    const auto& edge_lookup_table = detail::get_edge_lookup_table();

    if (state.is_existing_vertex_included.empty())
    {
        state.is_existing_vertex_included.resize(tet_vertices.rows(), false);
    }

    // For each tet
    for (int tet_index = 0; tet_index < tet_connectivity.rows(); tet_index++)
    {
        // The table key is effectively made of two bitsets:
        // |__|__|__|__|__|__|__|__|
        //   node i is | node i is
        //    == iso   |   > iso
        uint16_t key = 0;
        // > iso section
        for (int v = 0; v < 4; v++)
        {
            const auto& vid     = tet_connectivity(tet_index, v);
            const uint16_t flag = isovalues(vid, 0) > surface_isovalue;
            key |= flag << v;
        }
        // == iso section
        for (int v = 0; v < 4; v++)
        {
            const auto& vid    = tet_connectivity(tet_index, v);
            const uint8_t flag = isovalues(vid, 0) == surface_isovalue;
            key |= flag << (v + 4);
        }
        const auto& cell_entry = cell_lookup_table[key];

        // If: tetrahedron is intersected in some fashion
        if (cell_entry.vertex_types[0] != detail::VertexType::NA)
        {
            // This indices the vertices that form the triangle/quad along the isosurface
            // if the tetrahedron is intersected.  Negative indices will indicate new
            // vertices from a split edge, while positive (or zero) indices refer to
            // existing vertices in tet_vertices
            std::array<int, 4> vertices_on_surface;

            // Determine any vertices if the tet intersects the level surface
            for (auto i = 0; i < 4; i++)
            {
                if (cell_entry.vertex_types[i] == detail::VertexType::EDGE_INTERSECTION)
                {
                    // Note: determining the vertex position and index of the split
                    // edge is delayed until after the loop over the tetrahedra.
                    const auto& edge_i = cell_entry.intersected_edges[i];
                    auto v_1           = tet_connectivity(tet_index, edge_lookup_table[edge_i][0]);
                    auto v_2           = tet_connectivity(tet_index, edge_lookup_table[edge_i][1]);
                    if (v_2 < v_1)
                    {
                        std::swap(v_1, v_2);
                    }

                    auto edge_hash = detail::make_edge_key({ v_1, v_2 });
                    auto find_iter = state.split_edge_to_split_edge_vertex_map.find(edge_hash);
                    if (find_iter == state.split_edge_to_split_edge_vertex_map.end())
                    { // The edge has not be split yet, so create an additional vertex for it
                        const double a = fabs(isovalues(v_1, 0) - surface_isovalue);
                        const double b = fabs(isovalues(v_2, 0) - surface_isovalue);
                        const double w = a / (a + b);

                        state.split_edge_vertices.push_back((1 - w) * tet_vertices.row(v_1) +
                                                            w * tet_vertices.row(v_2));
                        state.split_edges.push_back({ v_1, v_2 });

                        const int additional_vertex_id = -1 * state.split_edge_vertices.size();
                        state.split_edge_to_split_edge_vertex_map[edge_hash] = additional_vertex_id;
                        vertices_on_surface[i]                               = additional_vertex_id;
                    }
                    else
                    {
                        vertices_on_surface[i] = find_iter->second;
                    }
                }
                else if (cell_entry.vertex_types[i] == detail::VertexType::EXISTING_VERTEX)
                {
                    vertices_on_surface[i] = tet_connectivity(tet_index, cell_entry.intersected_edges[i]);
                    state.is_existing_vertex_included[vertices_on_surface[i]] = true;
                }
            }

            // Insert the corresponding faces
            bool is_quad = cell_entry.intersected_edges[3] != -1;
            if (is_quad)
            {
                std::array<int, 3> face_sorted = { vertices_on_surface[0],
                                                   vertices_on_surface[1],
                                                   vertices_on_surface[3] };
                std::sort(face_sorted.begin(), face_sorted.end());
                if (state.unique_faces.count(face_sorted) == 0)
                {
                    state.unique_faces.insert(face_sorted);
                    const Eigen::RowVector3i face(vertices_on_surface[0],
                                                  vertices_on_surface[1],
                                                  vertices_on_surface[3]);
                    state.faces.push_back(std::make_pair(face, tet_index));
                }

                face_sorted = { vertices_on_surface[1],
                                vertices_on_surface[2],
                                vertices_on_surface[3] };
                std::sort(face_sorted.begin(), face_sorted.end());
                if (state.unique_faces.count(face_sorted) == 0)
                {
                    state.unique_faces.insert(face_sorted);
                    const Eigen::RowVector3i face(vertices_on_surface[1],
                                                  vertices_on_surface[2],
                                                  vertices_on_surface[3]);
                    state.faces.push_back(std::make_pair(face, tet_index));
                }
            }
            else
            {
                std::array<int, 3> face_sorted = { vertices_on_surface[0],
                                                   vertices_on_surface[1],
                                                   vertices_on_surface[2] };
                std::sort(face_sorted.begin(), face_sorted.end());
                if (state.unique_faces.count(face_sorted) == 0)
                {
                    state.unique_faces.insert(face_sorted);
                    const Eigen::RowVector3i face(vertices_on_surface[0],
                                                  vertices_on_surface[1],
                                                  vertices_on_surface[2]);
                    state.faces.push_back(std::make_pair(face, tet_index));
                }
            }
        }
    }
}

void recover_isosurface(const Eigen::MatrixXd& tet_vertices,
                        const MarchingTetState& state,
                        Eigen::MatrixXd& tri_vertices,
                        Eigen::MatrixXi& tri_connectivity,
                        Eigen::VectorXi& tri_vertex_to_tet_vertex_map,
                        Eigen::MatrixXi& tri_vertex_to_tet_split_edge_map)
{
    auto n_included_vertices = std::count(state.is_existing_vertex_included.begin(),
                                          state.is_existing_vertex_included.end(),
                                          true);
    auto n_vertices          = n_included_vertices + state.split_edge_vertices.size();
    tri_vertices.resize(n_vertices, 3);
    tri_connectivity.resize(state.faces.size(), 3);
    tri_vertex_to_tet_vertex_map.resize(n_vertices);
    tri_vertex_to_tet_vertex_map.setConstant(-1);
    tri_vertex_to_tet_split_edge_map.resize(n_vertices, 2);
    tri_vertex_to_tet_split_edge_map.setConstant(-1);

    int next_vertex_index = 0;
    std::vector<int> included_vertices_to_new_index(state.is_existing_vertex_included.size(), -1);
    for (auto i = 0; i < included_vertices_to_new_index.size(); i++)
    {
        if (state.is_existing_vertex_included[i])
        {
            tri_vertices.row(next_vertex_index)               = tet_vertices.row(i);
            included_vertices_to_new_index[i]                 = next_vertex_index;
            tri_vertex_to_tet_vertex_map(next_vertex_index++) = i;
        }
    }
    for (auto i = 0; i < state.split_edge_vertices.size(); i++)
    {
        auto& split_edge_vertex                                = state.split_edge_vertices[i];
        tri_vertices.row(next_vertex_index)                    = split_edge_vertex;
        tri_vertex_to_tet_split_edge_map(next_vertex_index, 0) = state.split_edges[i][0];
        tri_vertex_to_tet_split_edge_map(next_vertex_index, 1) = state.split_edges[i][1];
        next_vertex_index++;
    }
    assert(next_vertex_index == n_vertices);

    for (int f = 0; f < state.faces.size(); f++)
    {
        const int origin_tet_index = state.faces[f].second;

        const auto& face_vertices = state.faces[f].first;
        for (int v = 0; v < 3; v++)
        {
            if (face_vertices[v] < 0)
            {
                // This is a new vertex from a split edge (note: negative indices were used for
                // them in extended_marching_tets)
                const int new_vertex_index = n_included_vertices + (-1 * face_vertices[v]) - 1;
                tri_connectivity(f, v)     = new_vertex_index;
            }
            else
            {
                assert(included_vertices_to_new_index[face_vertices[v]] >= 0);
                tri_connectivity(f, v) = included_vertices_to_new_index[face_vertices[v]];
            }
        }
    }
}

void extended_marching_tets(const Tetrahedralization& t,
                            const std::vector<mesh::LocalIndex>& vertices_with_isovalues,
                            const Eigen::VectorXd& isovalues,
                            const double surface_isovalue,
                            std::vector<std::array<mesh::LocalIndex, 3>>& surface_triangles,
                            std::vector<mesh::LocalIndex>& intersected_edge_ids,
                            std::vector<std::array<double, 3>>& point_on_intersected_edges)
{
    MarchingTetState2 state;
    extended_marching_tets(t, vertices_with_isovalues, isovalues, surface_isovalue, state);
    recover_isosurface(t,
                       state,
                       surface_triangles,
                       intersected_edge_ids,
                       point_on_intersected_edges);
}

void extended_marching_tets(const Tetrahedralization& t,
                            const std::vector<mesh::LocalIndex>& vertices_with_isovalues,
                            const Eigen::VectorXd& isovalues,
                            const double surface_isovalue,
                            MarchingTetState2& state)
{
    assert(isovalues.cols() == 1);
    assert(vertices_with_isovalues.size() == isovalues.rows());

    const auto& cell_lookup_table = detail::get_cell_lookup_table();
    const auto& edge_lookup_table = detail::get_edge_lookup_table();

    if (state.edge_to_split_edge_vertex_map.empty())
    {
        state.edge_to_split_edge_vertex_map.resize(t.edges.size(), mesh::NOT_SET_LocalIndex);
    }

    // Create a map from vertex indices in the tetrahedralization to those in
    // vertices_with_isovalues
    std::vector<mesh::LocalIndex> tet_to_iso_vertex_id_map(t.vertices.size(),
                                                           mesh::NOT_SET_LocalIndex);
    for (auto i : utility::indices(vertices_with_isovalues))
    {
        tet_to_iso_vertex_id_map[vertices_with_isovalues[i]] = i;
    }

    // For each tet
    for (int tet_index = 0; tet_index < t.tetrahedra.size(); tet_index++)
    {
        // Only look at tetrahedra for which all vertices have an isovalue
        bool all_verts_have_isovalue = true;
        for (int v = 0; v < 4; v++)
        {
            all_verts_have_isovalue &=
              tet_to_iso_vertex_id_map[t.tetrahedra[tet_index][v]] != mesh::NOT_SET_LocalIndex;
        }

        if (all_verts_have_isovalue)
        {
            // The table key is effectively made of two bitsets:
            // |__|__|__|__|__|__|__|__|
            //   node i is | node i is
            //    == iso   |   > iso
            uint16_t key = 0;
            // > iso section
            for (int v = 0; v < 4; v++)
            {
                const auto& vid     = tet_to_iso_vertex_id_map[t.tetrahedra[tet_index][v]];
                const uint16_t flag = isovalues(vid, 0) > surface_isovalue;
                key |= flag << v;
            }
            // == iso section
            for (int v = 0; v < 4; v++)
            {
                const auto& vid    = tet_to_iso_vertex_id_map[t.tetrahedra[tet_index][v]];
                const uint8_t flag = isovalues(vid, 0) == surface_isovalue;
                key |= flag << (v + 4);
            }
            const auto& cell_entry = cell_lookup_table[key];

            // If: tetrahedron is intersected in some fashion
            if (cell_entry.vertex_types[0] != detail::VertexType::NA)
            {
                // These are the indices of the vertices that form the triangle/quad along the
                // isosurface if the tetrahedron is intersected.  Negative indices will indicate new
                // vertices from a split edge, while positive (or zero) indices refer to
                // existing vertices in tet_vertices
                std::array<mesh::LocalIndex, 4> vertices_on_surface;

                // Determine any vertices if the tet intersects the level surface
                for (auto i = 0; i < 4; i++)
                {
                    if (cell_entry.vertex_types[i] == detail::VertexType::EDGE_INTERSECTION)
                    {
                        // Note: determining the vertex position and index of the split
                        // edge is delayed until after the loop over the tetrahedra.
                        const auto& edge_i = cell_entry.intersected_edges[i];
                        auto v_1           = t.tetrahedra[tet_index][edge_lookup_table[edge_i][0]];
                        auto v_2           = t.tetrahedra[tet_index][edge_lookup_table[edge_i][1]];
                        if (v_2 < v_1)
                        {
                            std::swap(v_1, v_2);
                        }
                        auto iso_v_1 = tet_to_iso_vertex_id_map[v_1];
                        auto iso_v_2 = tet_to_iso_vertex_id_map[v_2];

                        auto global_edge_id = t.tet_to_edges_map[tet_index][edge_i];

                        // A vertex has not been created yet along this edge
                        if (state.edge_to_split_edge_vertex_map[global_edge_id] ==
                            mesh::NOT_SET_LocalIndex)
                        {
                            const double a = fabs(isovalues(iso_v_1, 0) - surface_isovalue);
                            const double b = fabs(isovalues(iso_v_2, 0) - surface_isovalue);
                            const double w = a / (a + b);

                            state.edge_to_split_edge_vertex_map[global_edge_id] =
                              state.split_edge_vertices.size();

                            using namespace simple_math::point_operators;
                            state.split_edge_vertices.push_back((1 - w) * t.vertices[v_1] +
                                                                w * t.vertices[v_2]);
                        }

                        const mesh::LocalIndex additional_vertex_id =
                          t.vertices.size() + state.edge_to_split_edge_vertex_map[global_edge_id];
                        vertices_on_surface[i] = additional_vertex_id;
                    }
                    else if (cell_entry.vertex_types[i] == detail::VertexType::EXISTING_VERTEX)
                    {
                        vertices_on_surface[i] = t.tetrahedra[tet_index][cell_entry.intersected_edges[i]];
                    }
                }

                // Insert the corresponding faces
                bool is_quad = cell_entry.intersected_edges[3] != -1;
                if (is_quad)
                {
                    std::array<mesh::LocalIndex, 3> face_sorted = { vertices_on_surface[0],
                                                                    vertices_on_surface[1],
                                                                    vertices_on_surface[3] };
                    std::sort(face_sorted.begin(), face_sorted.end());
                    if (state.unique_faces.count(face_sorted) == 0)
                    {
                        state.unique_faces.insert(face_sorted);
                        const std::array<mesh::LocalIndex, 3> face = { vertices_on_surface[0],
                                                                       vertices_on_surface[1],
                                                                       vertices_on_surface[3] };
                        state.faces.push_back(std::make_pair(face, tet_index));
                    }

                    face_sorted = { vertices_on_surface[1],
                                    vertices_on_surface[2],
                                    vertices_on_surface[3] };
                    std::sort(face_sorted.begin(), face_sorted.end());
                    if (state.unique_faces.count(face_sorted) == 0)
                    {
                        state.unique_faces.insert(face_sorted);
                        const std::array<mesh::LocalIndex, 3> face = { vertices_on_surface[1],
                                                                       vertices_on_surface[2],
                                                                       vertices_on_surface[3] };
                        state.faces.push_back(std::make_pair(face, tet_index));
                    }
                }
                else
                {
                    std::array<mesh::LocalIndex, 3> face_sorted = { vertices_on_surface[0],
                                                                    vertices_on_surface[1],
                                                                    vertices_on_surface[2] };
                    std::sort(face_sorted.begin(), face_sorted.end());
                    if (state.unique_faces.count(face_sorted) == 0)
                    {
                        state.unique_faces.insert(face_sorted);
                        const std::array<mesh::LocalIndex, 3> face = { vertices_on_surface[0],
                                                                       vertices_on_surface[1],
                                                                       vertices_on_surface[2] };
                        state.faces.push_back(std::make_pair(face, tet_index));
                    }
                }
            }
        }
    }
}
void recover_isosurface(const Tetrahedralization& t,
                        const MarchingTetState2& state,
                        std::vector<std::array<mesh::LocalIndex, 3>>& surface_triangles,
                        std::vector<mesh::LocalIndex>& intersected_edge_ids,
                        std::vector<std::array<double, 3>>& point_on_intersected_edges)
{
    surface_triangles.resize(state.faces.size());
    for (auto surf_tri_id : utility::indices(state.faces))
    {
        surface_triangles[surf_tri_id] = state.faces[surf_tri_id].first;
    }

    // Note: we have to play some games to make sure the order in intersected_edge_ids matches
    // the order of intersected edges
    std::vector<mesh::LocalIndex> isosurface_intersected_edge_ids(state.split_edge_vertices.size());
    for (auto edge_id : utility::indices(state.edge_to_split_edge_vertex_map))
    {
        if (state.edge_to_split_edge_vertex_map[edge_id] != mesh::NOT_SET_LocalIndex)
        {
            isosurface_intersected_edge_ids[state.edge_to_split_edge_vertex_map[edge_id]] = edge_id;
        }
    }
    for (auto edge_id : isosurface_intersected_edge_ids)
    {
        intersected_edge_ids.push_back(edge_id);
    }

    point_on_intersected_edges = state.split_edge_vertices;
}

}
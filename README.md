# CompMesh Core (`compmesh-core`)

`compmesh-core` provides core surface reconstruction, topological remeshing, marching tetrahedra surface extraction, mesh quality optimization, and tetrahedralization algorithms for [CompMesh](https://github.com/afrl/compmesh) within the **Scalable Hierarchical Enrichment Finite-element Framework (SHEFF)**.

## Capabilities

- **Tow Surface Reconstruction & VTMS Geometry (`comp_mesh_core/FiberBundle.h`, `comp_mesh_core/TowSurface.h`, `comp_mesh_core/TowSurfaceGeneration.h`)**:
  - Reads digital fiber bundles (`.dfbr`) and VTMS tow cross-sectional/surface formats (`.ctw`, `.stw`), adaptively identifies cross-sectional planes along complex 3D fiber trajectories, fits cross-sectional boundary splines, and generates closed triangular tow surfaces (`VTMSTowSurface`).
- **Topological Data Structures & Local Remeshing (`comp_mesh_core/Triangulation.h`, `comp_mesh_core/Tetrahedralization.h`, `comp_mesh_core/SplitEdges.h`, `comp_mesh_core/CollapseEdges.h`, `comp_mesh_core/Seams.h`)**:
  - Maintains bidirectional edge-triangle-tetrahedron incidence maps (`Triangulation`, `Tetrahedralization`), feature seam tracking, axis-aligned bounding-box surface clipping (`TriangulationClip.h`), edge splitting, edge collapsing, and surface quality optimization (`Optimize.h`).
- **Extended Marching Tetrahedra (`comp_mesh_core/MarchingTetrahedra.h`)**:
  - Extracts medial and level-set zero-contour surfaces across segmented tetrahedral meshes while respecting constrained seam topologies during interpenetration removal.
- **External Geometry & Tetrahedralization Wrappers (`comp_mesh_core/IGLHelpers.h`, `comp_mesh_core/FTetWildHelpers.h`)**:
  - Integrates `libigl` (exact predicates, CGAL self-intersection remeshing, and winding numbers) and `fTetWild` (`FloatTetwild`) robust multi-region tetrahedralization with `sheff-mesh` data structures.

## Project Structure

- `include/comp_mesh_core/`: Public headers for fiber bundles, tow surface generation, triangulation/tetrahedralization topology, marching tetrahedra, optimization, and `libigl`/`fTetWild` interfaces.
- `source/`: Compiled C++ implementations of surface generation, topological operations, and meshing wrappers.
- `external/`: Git submodules for [`libigl`](https://github.com/KeithBallard/libigl) and [`fTetWild`](https://github.com/KeithBallard/fTetWild).
- `unit_tests/`: Unit tests for tow surface generation, marching tetrahedra, edge splitting/collapsing, surface clipping, and mesh optimization.
- **CMake Package / Target**: `CompMeshCore` / `CompMeshCore::comp_mesh_core` (static library).

## Dependencies

- **Build System**: [`sheff-build-tools`](https://github.com/afrl/sheff-build-tools) (`SHEFF_TOOLS`), CMake 3.17+, C++17/C++20 compiler.
- **First-Party SHEFF Modules**: [`sheff-mesh`](https://github.com/afrl/sheff-mesh) (`Mesh::mesh`), [`sheff-fea-io`](https://github.com/afrl/sheff-fea-io) (`FeaIo::fea_io`), [`sheff-math`](https://github.com/afrl/sheff-math) (`SimpleMath::simple_math`), [`sheff-utility`](https://github.com/afrl/sheff-utility) (`Utility::utility`).
- **Second-Party SHEFF Modules**: [`sheff-second-party-mesh`](https://github.com/afrl/sheff-second-party-mesh) (`SecondPartyMesh::second_party_mesh`), [`sheff-second-party-metis`](https://github.com/afrl/sheff-second-party-metis) (`SecondPartyMetis::second_party_metis`).
- **Third-Party & System Libraries**: [`Eigen3`](https://github.com/eigenteam/eigen-git-mirror), [`Boost`](https://github.com/boostorg/boost) (`headers`, `thread`, `system`), [`libigl`](https://github.com/KeithBallard/libigl), [`fTetWild`](https://github.com/KeithBallard/fTetWild), CGAL, GMP (`libgmp-dev`), and MPFR (`libmpfr-dev`).

## Building

When building as part of [CompMesh](https://github.com/afrl/compmesh), `compmesh-core` is configured and built automatically via `make.py`.

To build and test standalone (with `SHEFF_TOOLS` configured and third-party modules built):

```bash
python3 make.py --config --build --third-party --mode=RelWithDebInfo
python3 make.py --config --build --make-tests --run-tests --mode=RelWithDebInfo
```

## Public Release Information

Distribution Statement A. Approved for public release: distribution is unlimited. Case #: AFRL-2026-4513

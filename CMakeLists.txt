cmake_minimum_required(VERSION 3.8)

cmake_policy(SET CMP0048 NEW)
project(SIMPLEMATH VERSION 1.0.0 LANGUAGES CXX)

#########################################
# Variables that configure this project #
#########################################

set(PUBLIC_MODULE_NAME "CompMeshCore")
set(MAIN_COMPONENT "comp_mesh_core")


set(SHEFF_FIRST_DEPENDS "Mesh::mesh")
set(SHEFF_SECOND_DEPENDS "")

# Third-party dependencies (length of all three lists must match).
# Since we cannot assume certain pattern for names, using multiple lists gives
# the flexibility needed for various packages.
set(SHEFF_THIRD_DEPENDS_name "eigen")
set(SHEFF_THIRD_DEPENDS_package "Eigen3")
set(SHEFF_THIRD_DEPENDS_target "Eigen3::Eigen")

# External dependencies need to be found from the system (such as MPI)
# toggle_var indicates the variable name that controls whether the package is included.
# The entry can be "1" for packages that are always included.
set(SHEFF_EXTERNAL_DEPENDS_toggle_var "")
# The name of the package in the call find_package(<package name>)
set(SHEFF_EXTERNAL_DEPENDS_package "")
# Any components to get from the package (use comma delimited for multiple components in the same package)
set(SHEFF_EXTERNAL_DEPENDS_components "")
# The targets to link from the package (use comma delimited for multiple targets in the same package)
set(SHEFF_EXTERNAL_DEPENDS_targets "")
# Extra compile definitions if the package is included
set(SHEFF_EXTERNAL_DEPENDS_compile_defs "")

macro(SPECIAL_TARGET_COMMANDS)
    # LIBIGL Options
    set(LIBIGL_BUILD_TESTS           OFF CACHE BOOL "" FORCE)
    set(LIBIGL_BUILD_TUTORIALS       OFF CACHE BOOL "" FORCE)
    set(LIBIGL_INSTALL               ON  CACHE BOOL "" FORCE)
    # Components to include
    set(LIBIGL_PREDICATES            ON  CACHE BOOL "" FORCE)
    set(LIBIGL_COPYLEFT_CGAL         ON  CACHE BOOL "" FORCE)
    # Components to exclude
    set(LIBIGL_EMBREE                OFF CACHE BOOL "" FORCE)
    set(LIBIGL_GLFW                  OFF CACHE BOOL "" FORCE)
    set(LIBIGL_IMGUI                 OFF CACHE BOOL "" FORCE)
    set(LIBIGL_OPENGL                OFF CACHE BOOL "" FORCE)
    set(LIBIGL_PNG                   OFF CACHE BOOL "" FORCE)
    set(LIBIGL_XML                   OFF CACHE BOOL "" FORCE)
    set(LIBIGL_COPYLEFT_COMISO       OFF CACHE BOOL "" FORCE)
    set(LIBIGL_COPYLEFT_TETGEN       OFF CACHE BOOL "" FORCE)
    set(LIBIGL_RESTRICTED_MATLAB     OFF CACHE BOOL "" FORCE)
    set(LIBIGL_RESTRICTED_MOSEK      OFF CACHE BOOL "" FORCE)
    set(LIBIGL_RESTRICTED_TRIANGLE   OFF CACHE BOOL "" FORCE)
    add_subdirectory(external/libigl)
    target_link_libraries(comp_mesh_core PRIVATE igl_core)
    target_link_libraries(comp_mesh_core PRIVATE igl_copyleft_cgal)

    #add_subdirectory(external/ftetwild)
    #target_link_libraries(comp_mesh_core PRIVATE FloatTetwild)
endmacro()

# Options
option( SHEFF_BuildTests        "Build the unit tests." OFF )
option( SHEFF_BuildPython       "Build Python SWIG module" OFF )
option( SHEFF_BuildFortran      "Build Fortran SWIG module" OFF )
option( SHEFF_USE_MPI           "Indicates if the MPI routines should be included in the build." OFF )

# For debugging, print generator, compiler, and options
message( STATUS "Gen: ${CMAKE_GENERATOR}, CXX: ${CMAKE_CXX_COMPILER_ID}, Mode: ${CMAKE_BUILD_TYPE}" )
message( STATUS "Python Bindings: ${SHEFF_BuildPython}, Fortran Bindings: ${SHEFF_BuildPython}" )
message( STATUS "Tests: ${SHEFF_BuildTests}" )

# Add cmake dir to module path to override FindEigen3
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

include("$ENV{SHEFF_TOOLS}/cmake/SHEFFModuleMain.cmake")
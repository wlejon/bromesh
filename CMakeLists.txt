cmake_minimum_required(VERSION 3.24)
project(bromesh LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(MSVC)
    add_compile_options(/MP /permissive- /utf-8)
    # Static CRT only when building standalone (not consumed as a subdirectory)
    if(PROJECT_IS_TOP_LEVEL)
        set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
    endif()
endif()

# --- bromath (shared math library, header-only) ---
# Resolution order: an already-loaded target (a parent project got there first),
# a standalone checkout at ../bromath (the multi-repo dev layout), then the
# third_party/bromath submodule (what a fresh `git clone --recursive` gets).
set(BROMATH_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../bromath" CACHE PATH
    "Standalone bromath checkout; falls back to the third_party/bromath submodule if absent")
if(NOT TARGET bromath)
    if(EXISTS "${BROMATH_DIR}/CMakeLists.txt")
        add_subdirectory("${BROMATH_DIR}" "${CMAKE_CURRENT_BINARY_DIR}/bromath" EXCLUDE_FROM_ALL)
    elseif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/third_party/bromath/CMakeLists.txt")
        add_subdirectory(third_party/bromath EXCLUDE_FROM_ALL)
    else()
        message(FATAL_ERROR
            "bromesh: bromath not found at ${BROMATH_DIR} or third_party/bromath. "
            "Run `git submodule update --init --recursive` or set -DBROMATH_DIR=<path>.")
    endif()
endif()

# --- Third-party dependencies ---

# meshoptimizer — simplification, vertex cache optimization, quantization
if(NOT TARGET meshoptimizer)
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/third_party/meshoptimizer/CMakeLists.txt")
        set(MESHOPT_BUILD_DEMO OFF CACHE BOOL "" FORCE)
        set(MESHOPT_BUILD_TOOLS OFF CACHE BOOL "" FORCE)
        set(MESHOPT_BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
        add_subdirectory(third_party/meshoptimizer EXCLUDE_FROM_ALL)
    else()
        message(STATUS "bromesh: meshoptimizer not found, simplification/optimization disabled")
    endif()
endif()

# par_shapes — parametric primitive generation (header-only, we compile it)
if(NOT TARGET par_shapes AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/third_party/par/par_shapes.h")
    add_library(par_shapes STATIC third_party/par/par_shapes.cpp)
    target_include_directories(par_shapes PUBLIC third_party/par)
    # par_shapes' index type defaults to uint16_t. Subdivision builds an
    # UNWELDED mesh of 3*ntriangles points before welding, so a geodesicSphere
    # past ~5 subdivisions (n=6 ≈ 246k intermediate points) overflows 16 bits
    # and returns corrupt topology. Widen to 32 bits. PUBLIC so par_primitives
    # (which reads par_shapes_mesh::triangles) agrees on the width; the struct
    # is all pointers/ints, so layout is unchanged either way.
    target_compile_definitions(par_shapes PUBLIC PAR_SHAPES_T=uint32_t)
    if(MSVC)
        target_compile_options(par_shapes PRIVATE /W0)
    else()
        target_compile_options(par_shapes PRIVATE -w)
    endif()
endif()
if(NOT TARGET par_shapes)
    message(STATUS "bromesh: par_shapes not found, parametric primitives disabled")
endif()

# tinygltf — glTF/GLB loading (header-only, we compile it)
if(NOT TARGET tinygltf AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/third_party/tinygltf/tiny_gltf.h")
    add_library(tinygltf STATIC third_party/impl/tinygltf_impl.cpp)
    target_include_directories(tinygltf PUBLIC third_party/tinygltf)
    # stb_image is bundled in third_party/tinygltf and is compiled into
    # tinygltf_impl.cpp so tinygltf can decode embedded textures; loadGLTF
    # then surfaces them via GltfScene::images. External images (.glb with
    # separate PNGs alongside) stay disabled for sandboxing.
    target_compile_definitions(tinygltf PUBLIC
        TINYGLTF_NO_INCLUDE_STB_IMAGE
        TINYGLTF_NO_INCLUDE_STB_IMAGE_WRITE
        TINYGLTF_NO_EXTERNAL_IMAGE
        TINYGLTF_NO_STB_IMAGE_WRITE
    )
    if(MSVC)
        target_compile_options(tinygltf PRIVATE /W0)
    else()
        target_compile_options(tinygltf PRIVATE -w)
    endif()
endif()
if(NOT TARGET tinygltf)
    message(STATUS "bromesh: tinygltf not found, glTF I/O disabled")
endif()

# V-HACD — convex decomposition (header-only, we compile it)
if(NOT TARGET vhacd AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/third_party/v-hacd/include/VHACD.h")
    add_library(vhacd STATIC third_party/impl/vhacd_impl.cpp)
    target_include_directories(vhacd PUBLIC third_party/v-hacd/include)
    if(MSVC)
        target_compile_options(vhacd PRIVATE /W0)
    else()
        target_compile_options(vhacd PRIVATE -w)
    endif()
endif()
if(NOT TARGET vhacd)
    message(STATUS "bromesh: V-HACD not found, convex decomposition disabled")
endif()

# xatlas — automatic UV unwrapping and atlas generation
if(NOT TARGET xatlas AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/third_party/xatlas/source/xatlas/xatlas.h")
    add_library(xatlas STATIC third_party/xatlas/source/xatlas/xatlas.cpp)
    target_include_directories(xatlas PUBLIC third_party/xatlas/source/xatlas)
    if(MSVC)
        target_compile_options(xatlas PRIVATE /W0)
    else()
        target_compile_options(xatlas PRIVATE -w)
    endif()
endif()
if(NOT TARGET xatlas)
    message(STATUS "bromesh: xatlas not found, UV unwrapping disabled")
endif()

# manifold — Boolean/CSG operations (built from source, no optional deps)
if(NOT TARGET manifold_lib AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/third_party/manifold/src/manifold.cpp")
    add_library(manifold_lib STATIC
        third_party/manifold/src/boolean3.cpp
        third_party/manifold/src/boolean_result.cpp
        third_party/manifold/src/constructors.cpp
        third_party/manifold/src/csg_tree.cpp
        third_party/manifold/src/edge_op.cpp
        third_party/manifold/src/face_op.cpp
        third_party/manifold/src/impl.cpp
        third_party/manifold/src/manifold.cpp
        third_party/manifold/src/minkowski.cpp
        third_party/manifold/src/polygon.cpp
        third_party/manifold/src/properties.cpp
        third_party/manifold/src/quickhull.cpp
        third_party/manifold/src/sdf.cpp
        third_party/manifold/src/smoothing.cpp
        third_party/manifold/src/sort.cpp
        third_party/manifold/src/subdivision.cpp
        third_party/manifold/src/tree2d.cpp
    )
    target_include_directories(manifold_lib PUBLIC third_party/manifold/include)
    target_include_directories(manifold_lib PRIVATE third_party/manifold/src)
    target_compile_features(manifold_lib PUBLIC cxx_std_17)
    target_compile_definitions(manifold_lib PUBLIC MANIFOLD_PAR=-1)
    if(MSVC)
        target_compile_options(manifold_lib PRIVATE /W0 /DNOMINMAX /bigobj /fp:precise)
    else()
        target_compile_options(manifold_lib PRIVATE -w)
    endif()
endif()
if(NOT TARGET manifold_lib)
    message(STATUS "bromesh: manifold not found, Boolean/CSG disabled")
endif()

# OSQP — quadratic programming solver (bounded biharmonic weights).
# Guard and link on osqpstatic — the target bromesh actually consumes — so a
# parent project that already loaded OSQP still gets BBW enabled here.
if(NOT TARGET osqpstatic AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/third_party/osqp/CMakeLists.txt")
    set(OSQP_BUILD_SHARED_LIB OFF CACHE BOOL "" FORCE)
    set(OSQP_BUILD_STATIC_LIB ON CACHE BOOL "" FORCE)
    set(OSQP_ENABLE_PRINTING OFF CACHE BOOL "" FORCE)
    set(OSQP_ENABLE_PROFILING OFF CACHE BOOL "" FORCE)
    set(OSQP_ENABLE_INTERRUPT OFF CACHE BOOL "" FORCE)
    set(OSQP_ALGEBRA_BACKEND "builtin" CACHE STRING "" FORCE)
    add_subdirectory(third_party/osqp EXCLUDE_FROM_ALL)
endif()
if(NOT TARGET osqpstatic)
    message(STATUS "bromesh: OSQP not found, bounded biharmonic weights disabled")
endif()

# OpenFBX — FBX file loading (we compile it from source)
if(NOT TARGET openfbx AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/third_party/OpenFBX/src/ofbx.h")
    add_library(openfbx STATIC
        third_party/OpenFBX/src/ofbx.cpp
        third_party/OpenFBX/src/libdeflate.c
    )
    target_include_directories(openfbx PUBLIC third_party/OpenFBX/src)
    if(MSVC)
        target_compile_options(openfbx PRIVATE /W0)
    else()
        target_compile_options(openfbx PRIVATE -w)
    endif()
endif()
if(NOT TARGET openfbx)
    message(STATUS "bromesh: OpenFBX not found, FBX loading disabled")
endif()

# --- bromesh library ---

add_library(bromesh STATIC
    # Isosurface
    src/isosurface/marching_cubes.cpp
    src/isosurface/surface_nets.cpp
    src/isosurface/dual_contouring.cpp
    src/isosurface/transvoxel.cpp

    # Discrete voxel
    src/voxel/greedy_mesh.cpp
    src/voxel/voxel_chunk.cpp

    # Primitives
    src/primitives/primitives.cpp

    # Manipulation
    src/manipulation/normals.cpp
    src/manipulation/simplify.cpp
    src/manipulation/simplify_attributes.cpp
    src/manipulation/split_components.cpp
    src/manipulation/weld.cpp
    src/manipulation/merge.cpp
    src/manipulation/repair.cpp
    src/manipulation/subdivide.cpp
    src/manipulation/smooth.cpp
    src/manipulation/remesh.cpp
    src/manipulation/skin.cpp
    src/manipulation/skin_transfer.cpp
    src/manipulation/shrinkwrap.cpp
    src/manipulation/transform.cpp
    src/manipulation/polygon.cpp
    src/manipulation/poly_mesh.cpp

    # Analysis
    src/analysis/bbox.cpp
    src/analysis/manifold_check.cpp
    src/analysis/convex_decomposition.cpp
    src/analysis/bake.cpp
    src/analysis/bake_texture.cpp
    src/analysis/bake_transfer.cpp
    src/analysis/sample.cpp
    src/analysis/raycast.cpp
    src/analysis/bvh.cpp
    src/analysis/intersect.cpp

    # Animation
    src/animation/pose.cpp
    src/animation/ik.cpp
    src/animation/retarget.cpp
    src/animation/locomotion.cpp

    # UV
    src/uv/projection.cpp

    # Primitives (par_shapes)
    src/primitives/par_primitives.cpp

    # Optimization
    src/optimization/optimize.cpp
    src/optimization/meshlets.cpp
    src/optimization/analyze.cpp
    src/optimization/spatial.cpp
    src/optimization/encode.cpp
    src/optimization/strips.cpp
    src/optimization/progressive.cpp

    # Manipulation (sweep)
    src/manipulation/sweep.cpp
    src/manipulation/bezier_sweep.cpp

    # Procedural
    src/procedural/lsystem.cpp
    src/procedural/lsystem_turtle.cpp
    src/procedural/space_colonization.cpp
    src/procedural/branches.cpp
    src/procedural/plants.cpp
    src/procedural/leaf_scatter.cpp
    src/procedural/obstacle_field.cpp

    # I/O
    src/io/obj.cpp
    src/io/vox.cpp
    src/io/stl.cpp
    src/io/ply.cpp
    src/io/splat_ply.cpp
    src/io/fbx.cpp

    # Reconstruction
    src/reconstruction/reconstruct.cpp

    # UV unwrapping
    src/uv/unwrap.cpp

    # UV metrics
    src/uv/uv_metrics.cpp

    # Boolean/CSG
    src/csg/boolean.cpp

    # Rigging
    src/rigging/rig_spec.cpp
    src/rigging/skeleton_fit.cpp
    src/rigging/voxel_bind.cpp
    src/rigging/auto_rig.cpp
    src/rigging/landmark_detect.cpp
    src/rigging/weight_smooth.cpp
    src/rigging/mesh_laplacian.cpp
    src/rigging/bone_heat.cpp
    src/rigging/bbw.cpp
    src/rigging/weighting.cpp
    src/rigging/skin_validate.cpp
)

target_include_directories(bromesh PUBLIC include)
target_include_directories(bromesh PRIVATE src)
target_compile_features(bromesh PUBLIC cxx_std_20)

target_link_libraries(bromesh PUBLIC bromath::bromath)

# OpenMP — parallelizes the hot procedural loops (leaf_scatter, primitives,
# manipulation). The `#pragma omp parallel for` sites are inert without it, so
# without this the scatter runs single-threaded. Optional: absent OpenMP the
# pragmas are simply ignored and everything still builds and runs (serially).
find_package(OpenMP QUIET)
if(OpenMP_CXX_FOUND)
    target_link_libraries(bromesh PUBLIC OpenMP::OpenMP_CXX)
    message(STATUS "bromesh: OpenMP enabled (${OpenMP_CXX_FLAGS})")
else()
    message(STATUS "bromesh: OpenMP not found, procedural loops run serially")
endif()

if(MSVC)
    # /W4 but silence C4068 (unknown pragma) so an OpenMP-less build is quiet
    target_compile_options(bromesh PRIVATE /W4 $<$<NOT:$<BOOL:${OpenMP_CXX_FOUND}>>:/wd4068>)
else()
    target_compile_options(bromesh PRIVATE -Wall -Wextra)
endif()

# Link optional dependencies. Feature availability keys off the TARGET —
# whether bromesh created it above or a parent project loaded it first —
# so the BROMESH_HAS_* compile definitions can never drift from what is
# actually linked.
if(TARGET meshoptimizer)
    target_link_libraries(bromesh PRIVATE meshoptimizer)
    target_compile_definitions(bromesh PUBLIC BROMESH_HAS_MESHOPTIMIZER=1)
endif()

if(TARGET par_shapes)
    target_link_libraries(bromesh PRIVATE par_shapes)
    target_compile_definitions(bromesh PUBLIC BROMESH_HAS_PAR_SHAPES=1)
endif()

if(TARGET tinygltf)
    target_link_libraries(bromesh PRIVATE tinygltf)
    target_compile_definitions(bromesh PUBLIC BROMESH_HAS_GLTF=1)
    target_sources(bromesh PRIVATE src/io/gltf.cpp)
endif()

if(TARGET vhacd)
    target_link_libraries(bromesh PRIVATE vhacd)
    target_compile_definitions(bromesh PUBLIC BROMESH_HAS_VHACD=1)
endif()

if(TARGET xatlas)
    target_link_libraries(bromesh PRIVATE xatlas)
    target_compile_definitions(bromesh PUBLIC BROMESH_HAS_XATLAS=1)
endif()

if(TARGET manifold_lib)
    target_link_libraries(bromesh PRIVATE manifold_lib)
    target_compile_definitions(bromesh PUBLIC BROMESH_HAS_MANIFOLD=1)
endif()

if(TARGET openfbx)
    target_link_libraries(bromesh PRIVATE openfbx)
    target_compile_definitions(bromesh PUBLIC BROMESH_HAS_OPENFBX=1)
endif()

if(TARGET osqpstatic)
    target_link_libraries(bromesh PRIVATE osqpstatic)
    target_compile_definitions(bromesh PUBLIC BROMESH_HAS_OSQP=1)
endif()

# Tests
option(BROMESH_TESTS "Build bromesh tests" ON)
if(BROMESH_TESTS AND PROJECT_IS_TOP_LEVEL)
    enable_testing()
    add_subdirectory(tests)
endif()

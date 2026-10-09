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

# --- Dependencies (cmake/bro_deps.cmake) ---
# Like every repo in the bro ecosystem: an existing target (a parent project
# got there first, and its choice wins) first; otherwise the ../<name> working
# tree beside the top-level project; otherwise the head of its main, fetched
# at configure. Third-party code is pinned to exact commits instead.
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/bro_deps.cmake)

# --- bromath (shared math library, header-only) ---
bro_dependency(bromath GITHUB wlejon/bromath)

# --- bronze, the JavaScript runtime the api binds ---
# It brings brass, its code generator, which bromesh's JIT links. The api is
# bound across bronze's C++ embed boundary, so it has to build against the
# bronze the surrounding tree builds against.
bro_dependency(bronze GITHUB wlejon/bronze TARGET bronze_runtime_shared
               OPTIONS BRONZE_BUILD_SHARED_RUNTIME=ON BRONZE_BUILD_TESTS=OFF)

# --- Third-party dependencies ---
# Pinned at the commits the former third_party/ submodules recorded. The
# header-only / no-CMake ones are fetched as source (SOURCE_ONLY) and compiled
# into targets of our own below.

# meshoptimizer — simplification, vertex cache optimization, quantization
bro_dependency(meshoptimizer GITHUB zeux/meshoptimizer REF dc9d09ed83e1004aef47a1c3c597e0ec64848a37 THIRD_PARTY
               OPTIONS MESHOPT_BUILD_DEMO=OFF MESHOPT_BUILD_TOOLS=OFF MESHOPT_BUILD_SHARED_LIBS=OFF)
bro_dependency(tinygltf GITHUB syoyo/tinygltf REF d8f3bd93f7bdcb7304017f4354d58f1ab85a5578 THIRD_PARTY SOURCE_ONLY)
bro_dependency(vhacd GITHUB kmammou/v-hacd REF f900e42361491f525262d4825e758845e4969897 THIRD_PARTY SOURCE_ONLY)
bro_dependency(xatlas GITHUB jpcy/xatlas REF f700c7790aaa030e794b52ba7791a05c085faf0c THIRD_PARTY SOURCE_ONLY)
bro_dependency(manifold GITHUB elalish/manifold REF d93943fc7ef4418f736330c00fe32411b136b1d1 THIRD_PARTY SOURCE_ONLY)
bro_dependency(OpenFBX GITHUB nem0/OpenFBX REF 9974840ed4f56cfa3b26f4a53be5e019254462de THIRD_PARTY SOURCE_ONLY)

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
if(NOT TARGET tinygltf AND EXISTS "${tinygltf_SOURCE_DIR}/tiny_gltf.h")
    add_library(tinygltf STATIC third_party/impl/tinygltf_impl.cpp)
    target_include_directories(tinygltf PUBLIC ${tinygltf_SOURCE_DIR})
    # stb_image is bundled in the tinygltf tree and is compiled into
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
if(NOT TARGET vhacd AND EXISTS "${vhacd_SOURCE_DIR}/include/VHACD.h")
    add_library(vhacd STATIC third_party/impl/vhacd_impl.cpp)
    target_include_directories(vhacd PUBLIC ${vhacd_SOURCE_DIR}/include)
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
if(NOT TARGET xatlas AND EXISTS "${xatlas_SOURCE_DIR}/source/xatlas/xatlas.h")
    add_library(xatlas STATIC ${xatlas_SOURCE_DIR}/source/xatlas/xatlas.cpp)
    target_include_directories(xatlas PUBLIC ${xatlas_SOURCE_DIR}/source/xatlas)
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
if(NOT TARGET manifold_lib AND EXISTS "${manifold_SOURCE_DIR}/src/manifold.cpp")
    add_library(manifold_lib STATIC
        ${manifold_SOURCE_DIR}/src/boolean3.cpp
        ${manifold_SOURCE_DIR}/src/boolean_result.cpp
        ${manifold_SOURCE_DIR}/src/constructors.cpp
        ${manifold_SOURCE_DIR}/src/csg_tree.cpp
        ${manifold_SOURCE_DIR}/src/edge_op.cpp
        ${manifold_SOURCE_DIR}/src/face_op.cpp
        ${manifold_SOURCE_DIR}/src/impl.cpp
        ${manifold_SOURCE_DIR}/src/manifold.cpp
        ${manifold_SOURCE_DIR}/src/minkowski.cpp
        ${manifold_SOURCE_DIR}/src/polygon.cpp
        ${manifold_SOURCE_DIR}/src/properties.cpp
        ${manifold_SOURCE_DIR}/src/quickhull.cpp
        ${manifold_SOURCE_DIR}/src/sdf.cpp
        ${manifold_SOURCE_DIR}/src/smoothing.cpp
        ${manifold_SOURCE_DIR}/src/sort.cpp
        ${manifold_SOURCE_DIR}/src/subdivision.cpp
        ${manifold_SOURCE_DIR}/src/tree2d.cpp
    )
    target_include_directories(manifold_lib PUBLIC ${manifold_SOURCE_DIR}/include)
    target_include_directories(manifold_lib PRIVATE ${manifold_SOURCE_DIR}/src)
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
set(CMAKE_WARN_DEPRECATED OFF)
bro_dependency(osqp GITHUB osqp/osqp REF 236713ce9a56c182ac3230d52108f952afce1523 THIRD_PARTY TARGET osqpstatic
               OPTIONS OSQP_BUILD_SHARED_LIB=OFF OSQP_BUILD_STATIC_LIB=ON OSQP_ENABLE_PRINTING=OFF
                       OSQP_ENABLE_PROFILING=OFF OSQP_ENABLE_INTERRUPT=OFF OSQP_ALGEBRA_BACKEND=builtin)
if(NOT TARGET osqpstatic)
    message(STATUS "bromesh: OSQP not found, bounded biharmonic weights disabled")
endif()

# OpenFBX — FBX file loading (we compile it from source)
if(NOT TARGET openfbx AND EXISTS "${OpenFBX_SOURCE_DIR}/src/ofbx.h")
    add_library(openfbx STATIC
        ${OpenFBX_SOURCE_DIR}/src/ofbx.cpp
        ${OpenFBX_SOURCE_DIR}/src/libdeflate.c
    )
    target_include_directories(openfbx PUBLIC ${OpenFBX_SOURCE_DIR}/src)
    if(MSVC)
        target_compile_options(openfbx PRIVATE /W0)
    else()
        target_compile_options(openfbx PRIVATE -w)
    endif()
endif()
if(NOT TARGET openfbx)
    message(STATUS "bromesh: OpenFBX not found, FBX loading disabled")
endif()

# draco — Google's mesh compression. The DECODER is what .drc files and
# glTF's KHR_draco_mesh_compression need; the encoder rides along in the same
# static lib. Draco's own CMake builds it (the vendored tree is upstream
# 1.5.7 minus testdata/plugins); EXCLUDE_FROM_ALL keeps its tools out.
if(NOT TARGET draco::draco AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/third_party/draco/CMakeLists.txt")
    set(DRACO_TESTS OFF CACHE BOOL "" FORCE)
    set(DRACO_INSTALL OFF CACHE BOOL "" FORCE)
    set(DRACO_TRANSCODER_SUPPORTED OFF CACHE BOOL "" FORCE)
    set(CMAKE_POLICY_DEFAULT_CMP0148 OLD)
    # Draco names its library `draco` on MSVC and `draco_static` elsewhere;
    # the draco::draco alias is the one spelling both branches provide.
    add_subdirectory(third_party/draco EXCLUDE_FROM_ALL)
endif()
if(NOT TARGET draco::draco)
    message(STATUS "bromesh: draco not found, Draco decode disabled")
endif()

# --- bromesh library ---

add_library(bromesh STATIC
    # Isosurface
    src/isosurface/marching_cubes.cpp
    src/isosurface/surface_nets.cpp
    src/isosurface/dual_contouring.cpp
    src/isosurface/transvoxel.cpp
    src/isosurface/jit/sdf_node.cpp
    src/isosurface/jit/sdf_compiler.cpp
    src/isosurface/jit/jit_mesher.cpp

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
    src/manipulation/poly_mesh_ops.cpp
    src/manipulation/poly_mesh_tess.cpp
    src/manipulation/lathe.cpp
    src/manipulation/splat_ops.cpp

    # Analysis
    src/analysis/bbox.cpp
    src/analysis/manifold_check.cpp
    src/analysis/convex_decomposition.cpp
    src/analysis/convex_hull.cpp
    src/analysis/bake.cpp
    src/analysis/bake_texture.cpp
    src/analysis/bake_transfer.cpp
    src/analysis/sample.cpp
    src/analysis/raycast.cpp
    src/analysis/bvh.cpp
    src/analysis/intersect.cpp
    src/analysis/sdf.cpp

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
    target_sources(bromesh PRIVATE
        src/io/gltf_load.cpp
        src/io/gltf_save.cpp
    )
endif()

if(TARGET draco::draco)
    target_link_libraries(bromesh PRIVATE draco::draco)
    target_compile_definitions(bromesh PUBLIC BROMESH_HAS_DRACO=1)
    target_sources(bromesh PRIVATE src/io/draco.cpp)
    # draco declares its include dirs PRIVATE, and draco_features.h is
    # generated as <top-level build>/draco/draco_features.h (draco's CMake
    # uses CMAKE_BINARY_DIR, deliberately, so a superbuild parent shares
    # one) — name both roots here.
    target_include_directories(bromesh PRIVATE
        third_party/draco/src
        "${CMAKE_BINARY_DIR}")
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

if(TARGET brass)
    target_link_libraries(bromesh PUBLIC brass)
    target_compile_definitions(bromesh PUBLIC BROMESH_HAS_BRASS_JIT=1)
endif()

# Standalone Bronze JavaScript API
option(BROMESH_ENABLE_API "Build Bronze JavaScript API" ON)
if(BROMESH_ENABLE_API)
    add_subdirectory(src/api)
endif()

# Tests
option(BROMESH_TESTS "Build bromesh tests" ON)
if(BROMESH_TESTS AND PROJECT_IS_TOP_LEVEL)
    enable_testing()
    add_subdirectory(tests)
endif()

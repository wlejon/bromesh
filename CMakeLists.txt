cmake_minimum_required(VERSION 3.24)
project(bromesh LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(MSVC)
    add_compile_options(/MP /permissive- /utf-8)
    # Static CRT only when building standalone (not consumed as a subdirectory)
    if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
        set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
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
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/third_party/par/par_shapes.h")
    add_library(par_shapes STATIC third_party/par/par_shapes.cpp)
    target_include_directories(par_shapes PUBLIC third_party/par)
    if(MSVC)
        target_compile_options(par_shapes PRIVATE /W0)
    else()
        target_compile_options(par_shapes PRIVATE -w)
    endif()
    set(BROMESH_HAS_PAR_SHAPES ON)
else()
    message(STATUS "bromesh: par_shapes not found, parametric primitives disabled")
    set(BROMESH_HAS_PAR_SHAPES OFF)
endif()

# tinygltf — glTF/GLB loading (header-only, we compile it)
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/third_party/tinygltf/tiny_gltf.h")
    add_library(tinygltf STATIC third_party/impl/tinygltf_impl.cpp)
    target_include_directories(tinygltf PUBLIC third_party/tinygltf)
    target_compile_definitions(tinygltf PUBLIC
        TINYGLTF_NO_INCLUDE_STB_IMAGE
        TINYGLTF_NO_INCLUDE_STB_IMAGE_WRITE
        TINYGLTF_NO_EXTERNAL_IMAGE
        TINYGLTF_NO_STB_IMAGE
        TINYGLTF_NO_STB_IMAGE_WRITE
    )
    if(MSVC)
        target_compile_options(tinygltf PRIVATE /W0)
    else()
        target_compile_options(tinygltf PRIVATE -w)
    endif()
    set(BROMESH_HAS_GLTF ON)
else()
    message(STATUS "bromesh: tinygltf not found, glTF I/O disabled")
    set(BROMESH_HAS_GLTF OFF)
endif()

# V-HACD — convex decomposition (header-only, we compile it)
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/third_party/v-hacd/include/VHACD.h")
    add_library(vhacd STATIC third_party/impl/vhacd_impl.cpp)
    target_include_directories(vhacd PUBLIC third_party/v-hacd/include)
    if(MSVC)
        target_compile_options(vhacd PRIVATE /W0)
    else()
        target_compile_options(vhacd PRIVATE -w)
    endif()
    set(BROMESH_HAS_VHACD ON)
else()
    message(STATUS "bromesh: V-HACD not found, convex decomposition disabled")
    set(BROMESH_HAS_VHACD OFF)
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

    # Primitives
    src/primitives/primitives.cpp

    # Manipulation
    src/manipulation/normals.cpp
    src/manipulation/simplify.cpp
    src/manipulation/split_components.cpp
    src/manipulation/weld.cpp

    # Analysis
    src/analysis/bbox.cpp
    src/analysis/manifold_check.cpp
    src/analysis/convex_decomposition.cpp

    # UV
    src/uv/projection.cpp

    # Optimization
    src/optimization/optimize.cpp

    # I/O
    src/io/obj.cpp
    src/io/vox.cpp
    src/io/stl.cpp
)

target_include_directories(bromesh PUBLIC include)
target_include_directories(bromesh PRIVATE src)
target_compile_features(bromesh PUBLIC cxx_std_20)

if(MSVC)
    target_compile_options(bromesh PRIVATE /W4)
else()
    target_compile_options(bromesh PRIVATE -Wall -Wextra)
endif()

# Link optional dependencies
if(TARGET meshoptimizer)
    target_link_libraries(bromesh PRIVATE meshoptimizer)
    target_compile_definitions(bromesh PUBLIC BROMESH_HAS_MESHOPTIMIZER=1)
endif()

if(BROMESH_HAS_PAR_SHAPES)
    target_link_libraries(bromesh PRIVATE par_shapes)
    target_compile_definitions(bromesh PUBLIC BROMESH_HAS_PAR_SHAPES=1)
endif()

if(BROMESH_HAS_GLTF)
    target_link_libraries(bromesh PRIVATE tinygltf)
    target_compile_definitions(bromesh PUBLIC BROMESH_HAS_GLTF=1)
    target_sources(bromesh PRIVATE src/io/gltf.cpp)
endif()

if(BROMESH_HAS_VHACD)
    target_link_libraries(bromesh PRIVATE vhacd)
    target_compile_definitions(bromesh PUBLIC BROMESH_HAS_VHACD=1)
endif()

# Tests
option(BROMESH_TESTS "Build bromesh tests" ON)
if(BROMESH_TESTS AND CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    enable_testing()
    add_subdirectory(tests)
endif()

add_library(bromesh_api STATIC
    api.cpp
    host_class.cpp
    native_mesh_core.cpp
    native_mesh_ops.cpp
    native_mesh_analysis.cpp
    native_mesh_plants.cpp
    native_mesh_io.cpp
    native_mesh_sdf.cpp
    native_rigging_core.cpp
    native_rigging_anim.cpp
)

target_include_directories(bromesh_api PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../../include>
)

target_link_libraries(bromesh_api PUBLIC
    bromesh
    bromath::bromath
    bronze_runtime_shared
)

target_compile_features(bromesh_api PUBLIC cxx_std_20)

if(MSVC)
    target_compile_options(bromesh_api PRIVATE /W4)
else()
    target_compile_options(bromesh_api PRIVATE -Wall -Wextra)
endif()

add_executable(bromesh_test
    test_main.cpp
    test_primitives.cpp
    test_isosurface.cpp
    test_voxel.cpp
    test_manipulation.cpp
    test_polymesh.cpp
    test_analysis.cpp
    test_bake.cpp
    test_uv.cpp
    test_optimization.cpp
    test_csg.cpp
    test_reconstruction.cpp
    test_io.cpp
    test_io_gltf.cpp
    test_io_draco.cpp
    test_skinning.cpp
    test_animation.cpp
    test_rigging.cpp
    test_rigging_multiped.cpp
    test_procedural.cpp
    test_coverage.cpp
    test_coverage2.cpp
    test_leaf_scatter.cpp
    test_splat.cpp
    test_sdf.cpp
)
target_link_libraries(bromesh_test PRIVATE bromesh)
# Tests poke at edge cases by hand-crafting glTF models via tinygltf directly.
if(TARGET tinygltf)
    target_link_libraries(bromesh_test PRIVATE tinygltf)
endif()

add_test(NAME bromesh_test COMMAND bromesh_test)

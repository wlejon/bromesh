add_executable(bromesh_test
    test_main.cpp
    test_primitives.cpp
    test_isosurface.cpp
    test_voxel.cpp
    test_manipulation.cpp
    test_polymesh.cpp
    test_lathe.cpp
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
    test_sdf_jit.cpp
)
target_link_libraries(bromesh_test PRIVATE bromesh)
# Tests poke at edge cases by hand-crafting glTF models via tinygltf directly.
if(TARGET tinygltf)
    target_link_libraries(bromesh_test PRIVATE tinygltf)
endif()
# Where the checked-in fixtures are (OpenFBX's sample b.fbx). The build tree
# is not necessarily inside the source tree — CI builds beside the checkout —
# so the tests must not guess at them from the working directory.
target_compile_definitions(bromesh_test PRIVATE
    BROMESH_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/..")

add_test(NAME bromesh_test COMMAND bromesh_test)

if(BROMESH_ENABLE_API AND TARGET bromesh_api)
    add_executable(test_mesh_api test_mesh_api.cpp test_mesh_api_restored.cpp test_mesh_api_validation.cpp)
    target_link_libraries(test_mesh_api PRIVATE bromesh_api bromesh bromath::bromath bronze_runtime_shared)
    # bronze's shared runtime is built into BRONZE_SHARED_RUNTIME_DIR, not
    # beside this executable, and the PE loader only looks beside the .exe
    # (0xc0000135 + a modal "dll was not found" otherwise). Stage it next to
    # the test as a post-build step and put both dirs on ctest's PATH.
    if(WIN32)
        add_custom_command(TARGET test_mesh_api POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
                $<TARGET_FILE:bronze_runtime_shared>
                $<TARGET_FILE_DIR:test_mesh_api>
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
                $<TARGET_LINKER_FILE:bronze_runtime_shared>
                $<TARGET_FILE_DIR:test_mesh_api>
        )
    endif()
    # It runs JavaScript, which brass compiles to machine code for the target;
    # on one brass has no code generator for (BRASS_HOST_BACKEND, decided in
    # brass's CMakeLists; x86_64 and AArch64 have one) that code cannot run. The
    # executable still builds there, which is the check that the embed links;
    # the test is not registered rather than registered to bus-error.
    if(NOT BRASS_HOST_BACKEND)
        message(STATUS "bromesh: test_mesh_api not registered (brass has no code generator for this target)")
    else()
        add_test(NAME test_mesh_api COMMAND test_mesh_api)
        # The same suite with bronze's collector flipping at every allocation:
        # a raw Value held across an allocating embed call in a binding
        # crashes or misreads here where the normal run happens to pass.
        add_test(NAME test_mesh_api_gcstress COMMAND test_mesh_api)
        if(WIN32)
            set_tests_properties(test_mesh_api PROPERTIES
                ENVIRONMENT "PATH=$<TARGET_FILE_DIR:bronze_runtime_shared>;$<TARGET_FILE_DIR:test_mesh_api>;$ENV{PATH}")
            set_tests_properties(test_mesh_api_gcstress PROPERTIES
                ENVIRONMENT "BRONZE_GC_STRESS=1;PATH=$<TARGET_FILE_DIR:bronze_runtime_shared>;$<TARGET_FILE_DIR:test_mesh_api>;$ENV{PATH}")
        else()
            set_tests_properties(test_mesh_api_gcstress PROPERTIES ENVIRONMENT "BRONZE_GC_STRESS=1")
        endif()
    endif()
endif()

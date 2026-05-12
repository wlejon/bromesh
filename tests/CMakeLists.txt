add_executable(bromesh_test
    test_main.cpp
    test_mesh.cpp
    test_animation.cpp
    test_rigging.cpp
    test_procedural.cpp
    test_coverage.cpp
    test_coverage2.cpp
)
target_link_libraries(bromesh_test PRIVATE bromesh)
# Tests poke at edge cases by hand-crafting glTF models via tinygltf directly.
if(TARGET tinygltf)
    target_link_libraries(bromesh_test PRIVATE tinygltf)
endif()

add_test(NAME bromesh_test COMMAND bromesh_test)

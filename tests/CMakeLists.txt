add_executable(bromesh_test
    test_main.cpp
    test_mesh.cpp
    test_animation.cpp
    test_rigging.cpp
    test_procedural.cpp
)
target_link_libraries(bromesh_test PRIVATE bromesh)

add_test(NAME bromesh_test COMMAND bromesh_test)

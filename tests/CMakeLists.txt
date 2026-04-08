add_executable(bromesh_test
    test_main.cpp
)
target_link_libraries(bromesh_test PRIVATE bromesh)

add_test(NAME bromesh_test COMMAND bromesh_test)

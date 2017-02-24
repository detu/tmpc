add_executable(test_matrix
    TestMatrix.cpp
    TestDynamicSubmatrix.cpp
)

target_include_directories(test_matrix
    PRIVATE "${CMAKE_SOURCE_DIR}/include"
)

add_dependencies(test_matrix googletest)

target_link_libraries(test_matrix
    gtest
    gtest_main
    gmock
    ${CMAKE_THREAD_LIBS_INIT}    # Needed for gtest
)
add_executable(test_matrix
    MatrixTest.cpp
    DynamicSubmatrixTest.cpp
    DynamicSubvectorTest.cpp
)

target_include_directories(test_matrix
    PRIVATE "${CMAKE_SOURCE_DIR}/include"
)

target_link_libraries(test_matrix
    GTest::GTest
    GTest::Main
    gmock
    ${CMAKE_THREAD_LIBS_INIT}    # Needed for gtest
)
find_package(OpenMP REQUIRED)

# Generate C code from CasADi functions
tmpc_generate_casadi_function(SRCS HDRS "sample_functions.py")

file(GLOB files "*.cpp")

foreach (file ${files})
    get_filename_component(file_we ${file} NAME_WE)

    add_executable(example-${file_we}
        ${file}
        ${SRCS} ${HDRS}
    )

    target_include_directories(example-${file_we}
        PRIVATE ${CMAKE_CURRENT_BINARY_DIR}
    )

    target_link_libraries(example-${file_we}
        tmpc-casadi
    )
endforeach ()

target_link_libraries(example-generated_function_parallel
    OpenMP::OpenMP_CXX
)
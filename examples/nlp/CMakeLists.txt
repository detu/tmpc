add_executable(sqp_rosenbrock
    sqp_rosenbrock.cpp
)

target_link_libraries(sqp_rosenbrock
    PRIVATE tmpc
)


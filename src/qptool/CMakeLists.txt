# Find Nlohmann json library
find_package(nlohmann_json REQUIRED)

add_executable(qptool qptool.cpp)
target_link_libraries(qptool 
    PRIVATE tmpc
    PRIVATE nlohmann_json
)

# Find Nlohmann json library
find_package(nlohmann_json REQUIRED)

add_library(tmpc-json INTERFACE)

target_link_libraries(tmpc-json
    INTERFACE tmpc
    INTERFACE nlohmann_json
)
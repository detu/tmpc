include_directories("json/src")

file(GLOB examples "qp/*.cpp" "ad/*.cpp")

foreach(example ${examples})
    get_filename_component(project_name ${example} NAME_WE)

    add_executable(${project_name} ${example})
    target_link_libraries(${project_name} PRIVATE tmpc)
endforeach()

#add_subdirectory(qp)
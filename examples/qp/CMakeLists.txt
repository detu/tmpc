set(files
    "ocp_tree.cpp"
    "ocp_size.cpp"
    "ocp_qp.cpp"
)

foreach (file ${files})
    get_filename_component(file_we ${file} NAME_WE)

    add_executable(example-${file_we}
        ${file}
        ${SRCS} ${HDRS}
    )

target_link_libraries(example-${file_we}
        tmpc
    )
endforeach ()


if (TMPC_WITH_HPIPM)
    set(files
        "solve_qp.cpp"
        # "soft_constraints.cpp"
    )

    foreach (file ${files})
        get_filename_component(file_we ${file} NAME_WE)

        add_executable(example-${file_we}
            ${file}
            ${SRCS} ${HDRS}
        )

        target_link_libraries(example-${file_we}
            tmpc-hpipm
        )
    endforeach ()
endif ()
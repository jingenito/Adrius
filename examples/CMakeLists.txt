foreach(example cf_expansion bosma_smeets)
    add_executable(${example} ${example}/main.cpp)
    target_link_libraries(${example} PRIVATE adrius::adrius)
    target_compile_options(${example} PRIVATE
        $<$<CXX_COMPILER_ID:MSVC>:/W4>
        $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wall -Wextra -Wpedantic>
    )
endforeach()

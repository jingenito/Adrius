include(GoogleTest)

set(ADRIUS_TEST_TARGETS
    test_concepts
    test_gram_schmidt
    test_lll
    test_continued_fraction
    test_rational
    test_illl
)

foreach(target IN LISTS ADRIUS_TEST_TARGETS)
    add_executable(${target} ${target}.cpp)
    target_link_libraries(${target} PRIVATE adrius::adrius GTest::gtest_main)
    target_compile_options(${target} PRIVATE
        $<$<CXX_COMPILER_ID:MSVC>:/W4 /WX>
        $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wall -Wextra -Wpedantic -Werror>
    )
    gtest_add_tests(TARGET ${target} AUTO)
endforeach()

if(ADRIUS_ENABLE_BOOST_MULTIPRECISION)
    add_executable(test_boost_backend test_boost_backend.cpp)
    target_link_libraries(test_boost_backend PRIVATE adrius::adrius GTest::gtest_main)
    target_compile_options(test_boost_backend PRIVATE
        $<$<CXX_COMPILER_ID:MSVC>:/W4 /WX>
        $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wall -Wextra -Wpedantic -Werror>
    )
    gtest_add_tests(TARGET test_boost_backend AUTO)
endif()

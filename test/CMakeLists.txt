# Fetch googletest
find_package(GTest REQUIRED)
include(GoogleTest)

#add_executable(twofloat_test pair-arithmetic.test.cpp double-word-arithmetic.test.cpp)
add_executable(twofloat_test double-word-arithmetic.test.cpp)
target_link_libraries(twofloat_test twofloat GTest::gtest_main)


gtest_discover_tests(twofloat_test)

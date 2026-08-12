# Fetch googletest
find_package(GTest REQUIRED)
include(GoogleTest)

add_executable(twofloat_test algorithms.test.cpp pair-arithmetic.test.cpp
                             double-word-arithmetic.test.cpp)
target_link_libraries(twofloat_test twofloat GTest::gtest_main)


gtest_discover_tests(twofloat_test)

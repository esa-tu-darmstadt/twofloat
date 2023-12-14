# Fetch googletest
INCLUDE(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG "v1.14.0"
)
FetchContent_MakeAvailable(googletest)

add_executable(twofloat_test pair-arithmetic.test.cpp double-word-arithmetic.test.cpp)
target_link_libraries(twofloat_test twofloat GTest::gtest_main)

include(GoogleTest)
gtest_discover_tests(twofloat_test)

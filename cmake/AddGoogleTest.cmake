# Adapted from the google test documentation
# Fetchs a specific version and make it locally available
set(GTEST_ARCHIVE_URL https://githUb.com/google/googletest/archive/)
# version - `release-1.12.1`
set(VERSION_COMMIT_HASH 609281088cfefc76f9d0ce82e1ff6c30cc3591e5)
set(EXTENSION .zip)

include(FetchContent)
FetchContent_Declare(
    googletest
    URL ${GTEST_ARCHIVE_URL}${VERSION_COMMIT_HASH}${EXTENSION} 
)

FetchContent_MakeAvailable(googletest)
enable_testing()

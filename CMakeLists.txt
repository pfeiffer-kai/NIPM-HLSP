cmake_minimum_required(VERSION 3.4)

project(NIPM_HLSP)

# set(CMAKE_BUILD_TYPE Release CACHE STRING "Release RelWithDebInfo Debug DEFAULT: Release" FORCE)
set(CMAKE_CXX_STANDARD 14 CACHE STRING "C++ version selection")

option(COMPILE_TESTS "compile the tests" ON)

add_compile_options(
  -O3 -g -Wall -Wpedantic -fPIC -std=c++14 -Wno-terminate
  )

################################################################################

file(GLOB SRC_Nipmhlsp nipm_hlsp/*.cpp)
file(GLOB INC_Nipmhlsp nipm_hlsp/*.h)

add_library(nipm_hlsp SHARED
    ${SRC_Nipmhlsp}
  )

################################################################################

find_package (Eigen3 3.2.10 REQUIRED NO_MODULE)
target_link_libraries(nipm_hlsp Eigen3::Eigen)

################################################################################

set_target_properties(nipm_hlsp PROPERTIES SOVERSION 0)

install(
	TARGETS nipm_hlsp
    EXPORT "${TARGETS_EXPORT_NAME}"
    LIBRARY DESTINATION "lib"
    ARCHIVE DESTINATION "lib"
    RUNTIME DESTINATION "bin"
    INCLUDES DESTINATION "${INCLUDE_INSTALL_DIR}"
)
install(FILES ${INC_Nipmhlsp} DESTINATION "include/nipm_hlsp")

if(${COMPILE_TESTS})
         add_subdirectory(tests)
endif()

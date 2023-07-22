# CMakeList.txt : CMake project for astrocytes_in_Cpp, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

project ("main")


file(GLOB source_files CONFIGURE_DEPENDS
            "*.h"
            "*.cpp"
            "src/*.cpp"
            "src/*.h"
            "src/algorithm/*.cpp"
            "src/algorithm/*.h"
            "src/cfgsys/*.cpp"
            "src/cfgsys/*.h"
            "util/*.cpp"
            "util/*.h"
)

# Add source to this project's executable.

add_executable (main ${source_files})

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET main PROPERTY CXX_STANDARD 20)
endif()

# TODO: Add tests and install targets if needed.
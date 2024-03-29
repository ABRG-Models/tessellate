#
# This is an example CMakeLists.txt file to compile a program against
# a morphologica that is 'in-tree'
#

# cmake version 3.1 provides the set(CMAKE_CXX_STANDARD 17) feature
cmake_minimum_required(VERSION 3.1)

# Give your project a name
project(Schnakenberg)

# From CMAKE_SYSTEM work out which of __OSX__, __GLN__, __NIX__ are required
if(CMAKE_SYSTEM MATCHES Linux.*)
  set(OS_FLAG "-D__GLN__")
elseif(CMAKE_SYSTEM MATCHES BSD.*)
  set(OS_FLAG "-D__NIX__")
elseif(APPLE)
  set(OS_FLAG "-D__OSX__")
else()
  message(ERROR "Operating system not supported: " ${CMAKE_SYSTEM})
endif()

# morphologica uses c++-17 language features
set(CMAKE_CXX_STANDARD 17)

# Add the host definition to CXXFLAGS along with other switches, depending on OS/Compiler
if (APPLE)
  set(CMAKE_CXX_FLAGS "${OS_FLAG} -Wall -Wfatal-errors -g -O3")
else()
  # This assumes a gcc compiler (or a gcc mimic like Clang)
  set(CMAKE_CXX_FLAGS "${OS_FLAG} -Wall -Wfatal-errors -g -O3 -Wno-unused-result -Wno-unknown-pragmas -march=native")
endif()

# Tell clang to be quiet about brace initialisers:
if(CMAKE_CXX_COMPILER_ID MATCHES Clang)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-missing-braces")
endif()

# Add OpenMP flags here, if necessary
find_package(OpenMP)
if(OpenMP_FOUND)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# Additional GL compiler flags.
#
set(OpenGL_GL_PREFERENCE "GLVND") # Following `cmake --help-policy CMP0072`
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DGL3_PROTOTYPES -DGL_GLEXT_PROTOTYPES")
if(APPLE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DGL_SILENCE_DEPRECATION")
endif()

# Tell the program where the morph fonts are, to compile them into the binary
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMORPH_FONTS_DIR=\"\\\"${PROJECT_SOURCE_DIR}/morphologica/fonts\\\"\"")

# Find the libraries which will be needed
find_package(jsoncpp REQUIRED)
find_package(HDF5 REQUIRED)
find_package(Armadillo REQUIRED)
find_package(OpenCV REQUIRED)
find_package(LAPACK REQUIRED)
find_package(OpenGL REQUIRED)
find_package(glfw3 3.3 REQUIRED)
find_package(Freetype REQUIRED)

# Define collections of includes that have to be made for morphologica
get_target_property(JSON_INC_PATH jsoncpp_lib INTERFACE_INCLUDE_DIRECTORIES)
set(MORPH_INC_CORE ${JSON_INC_PATH} ${ARMADILLO_INCLUDE_DIR} ${ARMADILLO_INCLUDE_DIRS} ${HDF5_INCLUDE_DIR})
set(MORPH_INC_GL ${OpenCV_INCLUDE_DIRS} ${OPENGL_INCLUDE_DIR} ${GLFW3_INCLUDE_DIR} ${FREETYPE_INCLUDE_DIRS})
include_directories(${MORPH_INC_CORE} ${MORPH_INC_GL})
# Assuming that you installed morphologica in-tree (i.e. 'next to' schnakenberg.cpp).
#set(MORPH_INCLUDE_PATH "$HOME/Neuroscience/DirichletRD/Threadbeast/sim/morphologica" CACHE PATH "The path to morphologica")
set(MORPH_INCLUDE_PATH "${PROJECT_SOURCE_DIR}/morphologica" CACHE PATH "The path to morphologica")
include_directories(BEFORE ${MORPH_INCLUDE_PATH}/include) # Allows GL3/gl3.h to be found
include_directories(BEFORE ${MORPH_INCLUDE_PATH})         # Allows morph/Header.h to be found

#
# All the library finding is now done. It's time to compile the program.
#
# Morphologica code requires a number of libraries, collected into 'CORE' and 'GL'.
set(MORPH_LIBS_GL ${OpenCV_LIBS} OpenGL::GL Freetype::Freetype glfw)
set(MORPH_LIBS_CORE ${ARMADILLO_LIBRARY} ${ARMADILLO_LIBRARIES} ${HDF5_C_LIBRARIES} ${LAPACK_LIBRARIES} jsoncpp_lib)


# The singleCircle system is compiled from a single c++ file
add_executable(schNakenberg schNakenberg.cpp)
target_compile_definitions(schNakenberg PUBLIC FLT=float SINGLE COMPILE_PLOTTING)
target_link_libraries(schNakenberg ${MORPH_LIBS_CORE} ${MORPH_LIBS_GL} )


# The displayCircles system is compiled from a single c++ file
add_executable(dCircles dCircles.cpp)
target_compile_definitions(dCircles PUBLIC FLT=float SINGLE COMPILE_PLOTTING)
target_link_libraries(dCircles ${MORPH_LIBS_CORE} ${MORPH_LIBS_GL} /usr/local/lib/libmorphdisplay0.so )


# The singleCircle system is compiled from a single c++ file
add_executable(sch sch.cpp)
target_compile_definitions(sch PUBLIC FLT=float SINGLE COMPILE_PLOTTING)
target_link_libraries(sch ${MORPH_LIBS_CORE} ${MORPH_LIBS_GL})


# The singleCircle system is compiled from a single c++ file
add_executable(setCentres setCentres.cpp)
target_compile_definitions(setCentres PUBLIC FLT=float SINGLE COMPILE_PLOTTING)
target_link_libraries(setCentres ${MORPH_LIBS_CORE} ${MORPH_LIBS_GL})

if(OpenMP_FOUND)
    target_link_libraries(schNakenberg OpenMP::OpenMP_CXX)
endif()
# For debugging of variables:
option(DEBUG_VARIABLES OFF)
if(DEBUG_VARIABLES)
  get_cmake_property(_variableNames VARIABLES)
  foreach (_variableName ${_variableNames})
    message(STATUS "${_variableName}=${${_variableName}}")
  endforeach()
endif(DEBUG_VARIABLES)

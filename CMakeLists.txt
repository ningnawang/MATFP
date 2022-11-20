cmake_minimum_required(VERSION 3.0.0)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

set(CGAL_DO_NOT_WARN_ABOUT_CMAKE_BUILD_TYPE TRUE)
add_definitions(-w -Wno-dev)

# macos
if(${CMAKE_HOST_SYSTEM_NAME} MATCHES "Darwin")
	set(MACOSX TRUE)
endif()

# used for MATFPDownloadExternal.cmake
option(AUTO_DOWNLOAD "Auto download dependencies" ON)

# ###############################################################################
# Options
# ###############################################################################
# matfp
option(MATFP_WITH_HUNTER "Use Hunter to download and configure Boost" OFF)

# libigl library
option(LIBIGL_USE_STATIC_LIBRARY "Use libigl as static library" OFF)
option(LIBIGL_WITH_ANTTWEAKBAR "Use AntTweakBar" OFF)
option(LIBIGL_WITH_CGAL "Use CGAL" ON)
option(LIBIGL_WITH_COMISO "Use CoMiso" OFF)
option(LIBIGL_WITH_CORK "Use Cork" OFF)
option(LIBIGL_WITH_LIM "Use LIM" OFF)
option(LIBIGL_WITH_MATLAB "Use Matlab" OFF)
option(LIBIGL_WITH_MOSEK "Use MOSEK" OFF)
option(LIBIGL_WITH_PNG "Use PNG" OFF)
option(LIBIGL_WITH_PYTHON "Use Python" OFF)
option(LIBIGL_WITH_TETGEN "Use Tetgen" OFF)
option(LIBIGL_WITH_TRIANGLE "Use Triangle" OFF)
option(LIBIGL_WITH_XML "Use XML" OFF)

# gui
option(LIBIGL_WITH_EMBREE "Use Embree" OFF)
option(LIBIGL_WITH_NANOGUI "Use Nanogui menu" OFF)
option(LIBIGL_WITH_OPENGL "Use OpenGL" OFF)
option(LIBIGL_WITH_OPENGL_GLFW "Use GLFW" OFF)
option(LIBIGL_WITH_VIEWER "Use OpenGL viewer" OFF)

# geogram
option(GEOGRAM_WITH_TRIANGLE "Use Triangle" OFF)

# ###############################################################################
# Project name
# ###############################################################################
project(MATFP)

# ###############################################################################
# Configure the compiler (after project)
# ###############################################################################
set(CMAKE_CXX_STANDARD 14)

# expected behaviour
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /std:c++lastest")

# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic-errors -Wno-missing-field-initializers -Wno-missing-braces -Wall -pthread -fsanitize=address -g -fno-omit-frame-pointer")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -Wall -pthread -D_DEBUG")

# ## Configure the compiler
# This is a basic, decent setup that should do something sane on most compilers
if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	# using Clang (linux or apple) or GCC
	message("Using clang/gcc compiler flags")
	SET(BASE_CXX_FLAGS "-std=c++14 -Wall -Wextra -g3 -pthread -fsanitize=address") # no use of -Werror here, too many errors
	SET(CMAKE_CXX_FLAGS_DEBUG "")
	SET(DISABLED_WARNINGS "-Wno-unused-parameter -Wno-unused-variable -Wno-unused-function -Wno-deprecated-declarations -Wno-missing-braces -Wno-non-virtual-dtor")
	SET(TRACE_INCLUDES " -H -Wno-error=unused-command-line-argument")

	if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
		message("Setting clang-specific options")
		SET(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} -ferror-limit=3 -fcolor-diagnostics")
		SET(CMAKE_CXX_FLAGS_DEBUG "-fsanitize=address -fno-limit-debug-info")
	elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
		message("Setting gcc-specific options")
		SET(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} -fmax-errors=5")
		SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} -Wno-maybe-uninitialized -Wno-format-zero-length -Wno-unused-but-set-parameter -Wno-unused-but-set-variable -Wno-sign-compare -Wno-reorder")
	endif()

	SET(CMAKE_CXX_FLAGS "${BASE_CXX_FLAGS} ${DISABLED_WARNINGS}")
	SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g3")

	# SET(CMAKE_CXX_FLAGS_RELEASE        "-O3 -march=native -DNDEBUG")
	# '-march=native' would cause M1 printing 'error: unknown target CPU 'vortex''
	SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
	# using Visual Studio C++
	message("Using Visual Studio compiler flags")
	set(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} /W4")
	set(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} /MP") # parallel build
	SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} /wd\"4267\"") # ignore conversion to smaller type (fires more aggressively than the gcc version, which is annoying)
	SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} /wd\"4244\"") # ignore conversion to smaller type (fires more aggressively than the gcc version, which is annoying)
	SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} /wd\"4305\"") # ignore truncation on initialization
	SET(CMAKE_CXX_FLAGS "${BASE_CXX_FLAGS} ${DISABLED_WARNINGS}")
	set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MD")
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MDd")

	add_definitions(/D "_CRT_SECURE_NO_WARNINGS")
	add_definitions(-DNOMINMAX)
	add_definitions(-D_USE_MATH_DEFINES)
else()
	# unrecognized
	message(FATAL_ERROR "Unrecognized compiler [${CMAKE_CXX_COMPILER_ID}]")
endif()

# ###############################################################################
# Settings
# ###############################################################################
if(NOT CMAKE_BUILD_TYPE)
	message(STATUS "No build type selected, default to Release")
	set(CMAKE_BUILD_TYPE "Release")
endif()

set(MATFP_EXTERNAL "${CMAKE_CURRENT_SOURCE_DIR}/extern")

# Color output
include(UseColors)

# Use folder in Visual Studio
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Extra warnings
include(Warnings)

# Export compile flags (used for autocompletion of the C++ code)
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

# Generate position independent code
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# ###############################################################################
# 3rd party libraries
# ###############################################################################
include(MATFPDependencies)
find_package(OpenGL REQUIRED COMPONENTS OpenGL)
find_package(GLUT REQUIRED)

# ###############################################################################
# MATFP
# ###############################################################################
# Build static library for executable
add_library(libMATFP
	include/matfp/Args.h
	include/matfp/Logger.h
	src/matfp/Common.cpp
	src/matfp/Common.h
	src/matfp/Logger.cpp

	src/matfp/geogram/predicates.cpp
	src/matfp/geogram/predicates.h
	src/matfp/geogram/mesh/mesh_sampling.h

	src/matfp/geogram/RPD.cpp
	src/matfp/geogram/RPD.h
	src/matfp/geogram/generic_RPD.h
	src/matfp/geogram/generic_RPD_vertex.h
	src/matfp/geogram/generic_RPD_polygon.cpp
	src/matfp/geogram/generic_RPD_polygon.h
	src/matfp/geogram/generic_RPD_cell.cpp
	src/matfp/geogram/generic_RPD_cell.h
	src/matfp/geogram/generic_RPD_utils.h
	src/matfp/geogram/RPD_mesh_builder.cpp
	src/matfp/geogram/RPD_mesh_builder.h
	src/matfp/geogram/RPD_callback.cpp
	src/matfp/geogram/RPD_callback.h

	src/matfp/NonManifoldMesh/Nonmanifoldmesh.cpp
	src/matfp/NonManifoldMesh/Nonmanifoldmesh.h
	src/matfp/NonManifoldMesh/Primitives.cpp
	src/matfp/NonManifoldMesh/Primitives.h

	src/matfp/MeshProcessor.cpp
	src/matfp/MeshProcessor.h
	src/matfp/LinearAlgebra/Wm4Math.cpp
	src/matfp/LinearAlgebra/Wm4Math.h
	src/matfp/LinearAlgebra/Wm4Matrix.cpp
	src/matfp/LinearAlgebra/Wm4Matrix.h
	src/matfp/LinearAlgebra/Wm4Vector.cpp
	src/matfp/LinearAlgebra/Wm4Vector.h
	src/matfp/GeometryObjects/GeometryObjects.cpp
	src/matfp/GeometryObjects/GeometryObjects.h
	src/matfp/GuiWindow.cpp
	src/matfp/GuiWindow.h
	src/matfp/ThreeDimensionalShape.cpp
	src/matfp/ThreeDimensionalShape.h
	src/matfp/WindingFilter.cpp
	src/matfp/WindingFilter.h
	src/matfp/AABBWrapper.cpp
	src/matfp/AABBWrapper.h
	src/matfp/MeshIO.cpp
	src/matfp/MeshIO.h
	src/matfp/external/predicates.c
	src/matfp/external/Predicates.hpp
	src/matfp/external/Predicates.cpp
	src/matfp/external/FastWindingNumber.hpp
	src/matfp/external/FastWindingNumber.cpp
	src/matfp/external/WindingNumber.h
	src/matfp/Parameters.h
	src/matfp/RPDGenerator.cpp
	src/matfp/RPDGenerator.h
	src/matfp/Triangulation.h
	src/matfp/Triangulation.cpp
	src/matfp/InscribedSpheres.h
	src/matfp/InscribedSpheres.cpp
	src/matfp/Types/DelaunayTriangulation.h
	src/matfp/Types/RegularTriangulation.h
	src/matfp/Types/CommonTypes.h
	src/matfp/Types/OtherTypes.h
	src/matfp/Types/OtherTypes.cpp
	src/matfp/MedialMeshGenerator.cpp
	src/matfp/MedialMeshGenerator.h
	src/matfp/Thinning.cpp
	src/matfp/Thinning.h
	src/matfp/UpdateSpheres.cpp
	src/matfp/UpdateSpheres.h
	src/matfp/ShrinkSpheres.cpp
	src/matfp/ShrinkSpheres.h
	src/matfp/IterateSpheres.cpp
	src/matfp/IterateSpheres.h
	src/matfp/FeaturePreservation.h
	src/matfp/FeaturePreservation.cpp
	src/matfp/InternalFeatureAddition.h
	src/matfp/InternalFeatureAddition.cpp
	src/matfp/MedialSpheresProcessor.h
	src/matfp/MedialSpheresProcessor.cpp
	src/matfp/LFS.h
	src/matfp/LFS.cpp
	src/matfp/MeshFeatureProcessor.h
	src/matfp/MeshFeatureProcessor.cpp
)

# ###############################################################################
# Include headers
# ###############################################################################
target_include_directories(libMATFP
	PRIVATE
	src
	PUBLIC
	include
)

# ###############################################################################
# Link libraries
# ###############################################################################
target_link_libraries(libMATFP
	PUBLIC
	geogram
	igl::core
	spdlog::spdlog
	polyscope
	PRIVATE
	igl::cgal
	warnings::all
)

set_target_properties(libMATFP PROPERTIES OUTPUT_NAME "matfp")
target_link_libraries(libMATFP PRIVATE OpenGL::GLU)

# ###############################################################################
# MATFP binary
# ###############################################################################

# Building executable
add_executable(MATFP src/main.cpp)
target_link_libraries(MATFP
	libMATFP
	CLI11::CLI11
	igl::cgal
	warnings::all
)
target_include_directories(MATFP PRIVATE src)
igl_copy_cgal_dll(MATFP)

# Install
install(TARGETS MATFP RUNTIME DESTINATION bin)

# ###############################################################################
# MATFP preprocess binary
# ###############################################################################
add_library(libMATFP_PRE
	include/matfp/Logger.h
	src/matfp/Logger.cpp
	src/pre/pre_types.h
	src/pre/pre_meshIO.h
	src/pre/pre_meshIO.cpp
	src/pre/sharp_feature_detection.h
	src/pre/sharp_feature_detection.cpp
	src/pre/pre_guiwindow.h
	src/pre/pre_guiwindow.cpp
)
target_include_directories(libMATFP_PRE PUBLIC include)
target_link_libraries(libMATFP_PRE
	PUBLIC
	geogram
	igl::core
	spdlog::spdlog
	polyscope
	PRIVATE
	warnings::all
)
add_executable(MATFP_PRE src/pre/pre_main.cpp)
target_link_libraries(MATFP_PRE
	libMATFP_PRE
	CLI11::CLI11
	warnings::all
)

# target_include_directories(MATFP_PRE PRIVATE src)

# ###############################################################################
# Folders for Visual Studio/XCode IDEs
# ###############################################################################

# geogram
set_target_properties(geogram PROPERTIES FOLDER extern/geogram)
set_target_properties(geogram_third_party PROPERTIES FOLDER extern/geogram)
set_target_properties(uninstall PROPERTIES FOLDER extern/geogram)
# ###############################################################################
# CMake download helpers
# ###############################################################################

# download external dependencies
include(MATFPDownloadExternal)

# ###############################################################################
# Required dependencies
# ###############################################################################

# openvdb
# if (NOT TARGET openvdb)
# matfp_download_openvdb()
# add_subdirectory(${MATFP_EXTERNAL}/openvdb)
# endif()

# geogram
if(NOT TARGET geogram)
	matfp_download_geogram()
	include(geogram)
endif()

# Boost
if(MATFP_WITH_HUNTER)
	hunter_add_package(Boost COMPONENTS thread system)
endif()

# fmt
if(NOT TARGET fmt::fmt)
	matfp_download_fmt()
	add_subdirectory(${MATFP_EXTERNAL}/fmt)
endif()

# spdlog
if(NOT TARGET spdlog::spdlog)
	matfp_download_spdlog()
	add_library(spdlog INTERFACE)
	add_library(spdlog::spdlog ALIAS spdlog)
	target_include_directories(spdlog INTERFACE ${MATFP_EXTERNAL}/spdlog/include)
	target_compile_definitions(spdlog INTERFACE -DSPDLOG_FMT_EXTERNAL)
	target_link_libraries(spdlog INTERFACE fmt::fmt)
endif()

# libigl
if(NOT TARGET igl::core)
	matfp_download_libigl()
	find_package(LIBIGL REQUIRED)
endif()

# CL11
if(NOT TARGET CLI11::CLI11)
	matfp_download_cli11()
	add_subdirectory(${MATFP_EXTERNAL}/cli11)
endif()

# # winding number
# matfp_download_windingnumber()
# set(windingnumber_SOURCES
# ${MATFP_EXTERNAL}/windingnumber/SYS_Math.h
# ${MATFP_EXTERNAL}/windingnumber/SYS_Types.h
# ${MATFP_EXTERNAL}/windingnumber/UT_Array.cpp
# ${MATFP_EXTERNAL}/windingnumber/UT_Array.h
# ${MATFP_EXTERNAL}/windingnumber/UT_ArrayImpl.h
# ${MATFP_EXTERNAL}/windingnumber/UT_BVH.h
# ${MATFP_EXTERNAL}/windingnumber/UT_BVHImpl.h
# ${MATFP_EXTERNAL}/windingnumber/UT_FixedVector.h
# ${MATFP_EXTERNAL}/windingnumber/UT_ParallelUtil.h
# ${MATFP_EXTERNAL}/windingnumber/UT_SmallArray.h
# ${MATFP_EXTERNAL}/windingnumber/UT_SolidAngle.cpp
# ${MATFP_EXTERNAL}/windingnumber/UT_SolidAngle.h
# ${MATFP_EXTERNAL}/windingnumber/VM_SIMD.h
# ${MATFP_EXTERNAL}/windingnumber/VM_SSEFunc.h
# )

# add_library(fast_winding_number ${windingnumber_SOURCES})
# if(MACOSX)
# target_link_libraries(fast_winding_number PRIVATE TBB::tbb)
# else() # for Linux
# target_link_libraries(fast_winding_number PRIVATE tbb)
# endif()
# target_compile_features(fast_winding_number PRIVATE ${CXX17_FEATURES})
# target_include_directories(fast_winding_number PUBLIC "${MATFP_EXTERNAL}/")

# polyscope
matfp_download_polyscope()
add_subdirectory(${MATFP_EXTERNAL}/polyscope)
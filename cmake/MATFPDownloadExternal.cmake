# ###############################################################################
include(DownloadProject)

# Shortcut function
function(matfp_download_project name)
    if(AUTO_DOWNLOAD)
        download_project(
            PROJ ${name}
            SOURCE_DIR ${MATFP_EXTERNAL}/${name}
            DOWNLOAD_DIR ${MATFP_EXTERNAL}/.cache/${name}
            ${ARGN}
        )
    endif()
endfunction()

# ###############################################################################

# # openvdb
function(matfp_download_openvdb)
    matfp_download_project(openvdb
        GIT_REPOSITORY https://github.com/AcademySoftwareFoundation/openvdb.git
        GIT_TAG v7.1.0
    )
endfunction()

# # libigl
function(matfp_download_libigl)
    matfp_download_project(libigl
        GIT_REPOSITORY https://github.com/libigl/libigl.git
        GIT_TAG v2.3.0
    )
endfunction()

# # geogram
function(matfp_download_geogram)
    matfp_download_project(geogram
        GIT_REPOSITORY https://github.com/alicevision/geogram.git

        # GIT_TAG        v1.6.7
        GIT_TAG v1.7.5
    )
endfunction()

# # fmt
function(matfp_download_fmt)
    matfp_download_project(fmt
        GIT_REPOSITORY https://github.com/fmtlib/fmt.git
        GIT_TAG 5.2.0
    )
endfunction()

# # spdlog
function(matfp_download_spdlog)
    matfp_download_project(spdlog
        GIT_REPOSITORY https://github.com/gabime/spdlog.git
        GIT_TAG v1.1.0
    )
endfunction()

# # CLI11
function(matfp_download_cli11)
    matfp_download_project(cli11
        GIT_REPOSITORY https://github.com/CLIUtils/CLI11
        GIT_TAG v1.6.1
    )
endfunction()

# ## winding number
# function(matfp_download_windingnumber)
# matfp_download_project(windingnumber
# GIT_REPOSITORY https://github.com/alecjacobson/WindingNumber.git
# GIT_TAG        1e6081e52905575d8e98fb8b7c0921274a18752f
# )
# endfunction()

# polyscope
function(matfp_download_polyscope)
    matfp_download_project(polyscope
        GIT_REPOSITORY https://github.com/nmwsharp/polyscope.git
        GIT_TAG v1.3.0
    )
endfunction()

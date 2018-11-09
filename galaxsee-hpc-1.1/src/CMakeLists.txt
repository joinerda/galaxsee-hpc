include_directories(${galaxseehpc_SOURCE_DIR}/src)
include_directories(${galaxseehpc_BINARY_DIR}/src)
link_directories(${galaxseehpc_BINARY_DIR}/src)

set (galaxseehpc_VERSION_MAJOR 1)
set (galaxseehpc_VERSION_MINOR 1)

set(EXECNAME galaxsee)
add_executable(${EXECNAME} cubeinterp.c fcr.c gal2gad2.c galaxsee.c
      mem.c nbody.c octtree.c pppm.c quaternion.c
      rand_tools.c readline.c sdlwindow.c text11.c timer.c)
option(USE_MPI "Compile mpi version of library" OFF)
if(USE_MPI)
    find_package(MPI REQUIRED)
    mark_as_advanced(CLEAR MPI_INCLUDE_PATH
        MPI_COMPILE_FLAGS MPI_LINK_FLAGS
        MPI_LIBRARY MPI_EXTRA_LIBRARY)
    SET(CMAKE_C_COMPILER mpicc)
    set(CMAKE_C_FLAGS ${CMAKE_C_FLAGS}
        ${MPI_COMPILE_FLAGS})
    add_definitions(-DHAS_MPI)
    include_directories(MPI_INCLUDE_PATH)
    target_link_libraries(${EXECNAME} ${MPI_LIBRARIES})
else(USE_MPI)
endif(USE_MPI)

configure_file (
    "${galaxseehpc_SOURCE_DIR}/src/nbody_config.h.in"
    "${galaxseehpc_BINARY_DIR}/src/nbody_config.h"
)

set(CMAKE_VERBOSE_MAKEFILE OFF)

option(USE_FFTW "FFTW3 required for PPPM options" OFF)
if(USE_FFTW)
   FIND_PATH(FFTW_INCLUDE_DIR fftw3.h
      /usr/local/include
      /usr/include
      /opt/local/include
      /share/apps/fftw/3/3.2.2/include
   )
     
   SET(FFTW_NAMES ${FFTW_NAMES} fftw3 fftw3f fftw3-3)
   FIND_LIBRARY(FFTW_LIBRARY
      NAMES ${FFTW_NAMES}
      PATHS /usr/lib /usr/local/lib /opt/local/lib
      /share/apps/fftw/3/3.2.2/lib
   )
   include_directories(${FFTW_INCLUDE_DIR})
   target_link_libraries(${EXECNAME} ${FFTW_LIBRARY})
   add_definitions(-DHAS_FFTW3)
   add_definitions(-DUSE_PPPM)
endif(USE_FFTW)


option(USE_NETCDF "NETCDF Output Options" OFF)
if(USE_NETCDF)
   FIND_PATH(NETCDF_INCLUDE_DIR netcdf.h
      /usr/local/include
      /usr/include
      /opt/local/include
      /share/apps/NETCDF/4.1.2/include
   )
     
   SET(NETCDF_NAMES ${NETCDF_NAMES} netcdf)
   FIND_LIBRARY(NETCDF_LIBRARY
      NAMES ${NETCDF_NAMES}
      PATHS /usr/lib /usr/local/lib /opt/local/lib
/share/apps/NETCDF/4.1.2/lib
   )
   include_directories(${NETCDF_INCLUDE_DIR})
   target_link_libraries(${EXECNAME} ${NETCDF_LIBRARY})
   add_definitions(-DHAS_NETCDF)
endif(USE_NETCDF)

option(USE_HDF5 "HDF5 used for Gadget2 conversion Options" OFF)
if(USE_HDF5)
   FIND_PATH(HDF5_INCLUDE_DIR hdf5.h
      /usr/local/include
      /usr/include
      /opt/local/include
      /share/apps/hdf5/gnu/include
   )
     
   SET(HDF5_NAMES ${HDF5_NAMES} hdf5)
   FIND_LIBRARY(HDF5_LIBRARY
      NAMES ${HDF5_NAMES}
      PATHS /usr/lib /usr/local/lib /opt/local/lib
/share/apps/hdf5/gnu/lib
   )
   SET(HDF5_HL_NAMES ${HDF5_HL_NAMES} hdf5_hl)
   FIND_LIBRARY(HDF5_HL_LIBRARY
      NAMES ${HDF5_HL_NAMES}
      PATHS /usr/lib /usr/local/lib /opt/local/lib
/share/apps/hdf5/gnu/lib
   )
   include_directories(${HDF5_INCLUDE_DIR})
   target_link_libraries(${EXECNAME} ${HDF5_LIBRARY} ${HDF5_HL_LIBRARY})
   add_definitions(-DHAS_HDF5)
endif(USE_HDF5)

option(USE_LIBGD "LIBGD Output Options" OFF)
if(USE_LIBGD)
   FIND_PATH(LIBGD_INCLUDE_DIR gd.h
      /usr/local/include
      /usr/include
      /opt/local/include
      /share/apps/libgd/include
   )
     
   SET(LIBGD_NAMES ${LIBGD_NAMES} gd)
   FIND_LIBRARY(LIBGD_LIBRARY
      NAMES ${LIBGD_NAMES}
      PATHS /usr/lib /usr/local/lib /opt/local/lib
      /share/apps/libgd/lib
   )
   include_directories(${LIBGD_INCLUDE_DIR})
     
   SET(LIBPNG_NAMES ${LIBPNG_NAMES} png)
   FIND_LIBRARY(LIBPNG_LIBRARY
      NAMES ${LIBPNG_NAMES}
      PATHS /usr/lib /usr/local/lib
   )
   SET(LIBZ_NAMES ${LIBZ_NAMES} z)
   FIND_LIBRARY(LIBZ_LIBRARY
      NAMES ${LIBZ_NAMES}
      PATHS /usr/lib /usr/local/lib
   )
   SET(LIBJPEG_NAMES ${LIBJPEG_NAMES} jpeg)
   FIND_LIBRARY(LIBJPEG_LIBRARY
      NAMES ${LIBJPEG_NAMES}
      PATHS /usr/lib /usr/local/lib
   )
   target_link_libraries(${EXECNAME} ${LIBGD_LIBRARY} ${LIBJPEG_LIBRARY} ${LIBPNG_LIBRARY}
${LIBZ_LIBRARY})
   add_definitions(-DHAS_LIBGD)
endif(USE_LIBGD)

option(USE_X11 "X11 Output Options" OFF)
if(USE_X11)
   FIND_PATH(X11_INCLUDE_DIR X11/Xlib.h
      /usr/local/include
      /usr/include
   )
     
   SET(X11_NAMES ${X11_NAMES} X11)
   FIND_LIBRARY(X11_LIBRARY
      NAMES ${X11_NAMES}
      PATHS /usr/lib /usr/local/lib /usr/lib64 /usr/local/lib64
   )
   include_directories(${X11_INCLUDE_DIR})
   target_link_libraries(${EXECNAME} ${X11_LIBRARY})
   add_definitions(-DHAS_X11)
endif(USE_X11)

option(USE_SDL "SDL Output Options" OFF)
if(USE_SDL)
   FIND_PATH(SDL_INCLUDE_DIR SDL/SDL.h
      /usr/local/include
      /usr/include
   )
     
   SET(SDL_NAMES ${SDL_NAMES} SDL SDL-1.2)
   FIND_LIBRARY(SDL_LIBRARY
      NAMES ${SDL_NAMES}
      PATHS /usr/lib /usr/local/lib /usr/lib64 /usr/local/lib64
   )
   SET(GL_NAMES ${GL_NAMES} GL)
   FIND_LIBRARY(GL_LIBRARY
      NAMES ${GL_NAMES}
      PATHS /usr/lib /usr/local/lib /usr/lib64 /usr/local/lib64
   )
   SET(GLU_NAMES ${GLU_NAMES} GLU)
   FIND_LIBRARY(GLU_LIBRARY
      NAMES ${GLU_NAMES}
      PATHS /usr/lib /usr/local/lib /usr/lib64 /usr/local/lib64
   )
   include_directories(${SDL_INCLUDE_DIR})
   target_link_libraries(${EXECNAME} ${SDL_LIBRARY} ${GL_LIBRARY}
${GLU_LIBRARY})
   add_definitions(-DHAS_SDL)
endif(USE_SDL)



if(UNIX)
  target_link_libraries(${EXECNAME} m)
  target_link_libraries(${EXECNAME} rt)
endif(UNIX)

install(TARGETS ${EXECNAME}
    RUNTIME DESTINATION bin)

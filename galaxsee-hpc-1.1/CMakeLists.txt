cmake_minimum_required (VERSION 2.6)
project(galaxseehpc)
add_subdirectory(src)
add_subdirectory(utils)
IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  SET(CMAKE_INSTALL_PREFIX
    "/usr/local/galaxsee" CACHE PATH "FOO install prefix" FORCE
    )
ENDIF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

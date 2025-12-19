cmake_minimum_required( VERSION 3.12 )
project( find_targets
  DESCRIPTION "Find package targets"
  VERSION 1.0 )

get_property(_before GLOBAL PROPERTY TARGETS)

find_package( ${PACKAGE} REQUIRED)
# find_package( PkgConfig REQUIRED )
# pkg_check_modules( ${PACKAGE} REQUIRED ${PACKAGE} )

get_property(_after GLOBAL PROPERTY TARGETS)

foreach(t ${_after})
  list(FIND _before ${t} _idx)
  if(_idx EQUAL -1)
    message(STATUS "New target from MyPkg: ${t}")
  endif()
endforeach()


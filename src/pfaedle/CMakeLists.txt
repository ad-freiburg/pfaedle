file(GLOB_RECURSE pfaedle_SRC *.cpp)

set(pfaedle_main PfaedleMain.cpp)

list(REMOVE_ITEM pfaedle_SRC ${pfaedle_main})

include_directories(
	${PFAEDLE_INCLUDE_DIR}
)

configure_file (
  "_config.h.in"
  "_config.h"
)

add_executable(pfaedle ${pfaedle_main})
add_library(pfaedle_dep ${pfaedle_SRC})

target_include_directories(pfaedle_dep PUBLIC ${PROJECT_SOURCE_DIR}/src/xml/include/ ${PROJECT_SOURCE_DIR}/src/cppgtfs/src)

target_link_libraries(pfaedle_dep PRIVATE pb_util pb_util_geo pb_util_graph pb_util_xml pb_util_json configparser ad_cppgtfs -lpthread)

target_link_libraries(pfaedle PRIVATE pfaedle_dep pb_util pb_util_geo pb_util_graph pb_util_xml pb_util_json configparser ad_cppgtfs -lpthread)

if (LIBZIP_FOUND)
	target_link_libraries(pfaedle PRIVATE ${LIBZIP_LIBRARY})
endif( LIBZIP_FOUND )

if (BZIP2_FOUND)
	target_link_libraries(pfaedle PRIVATE ${BZIP2_LIBRARY})
endif( BZIP2_FOUND )

if(CMAKE_TESTING_ENABLED)
	add_subdirectory(tests)
endif()

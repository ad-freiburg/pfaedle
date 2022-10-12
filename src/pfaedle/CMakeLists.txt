file(GLOB_RECURSE pfaedle_SRC *.cpp)

set(pfaedle_main PfaedleMain.cpp)

list(REMOVE_ITEM pfaedle_SRC ${pfaedle_main})

include_directories(
	${PFAEDLE_INCLUDE_DIR}
	SYSTEM ${LIBZIP_INCLUDE_DIR}
	SYSTEM ${LIBZIP_CONF_INCLUDE_DIR}
)

configure_file (
  "_config.h.in"
  "_config.h"
)

add_executable(pfaedle ${pfaedle_main})
add_library(pfaedle_dep ${pfaedle_SRC})

include_directories(pfaedle_dep PUBLIC ${PROJECT_SOURCE_DIR}/src/cppgtfs/src)
target_link_libraries(pfaedle pfaedle_dep util configparser ad_cppgtfs -lpthread ${LIBZIP_LIBRARY})

add_subdirectory(tests)

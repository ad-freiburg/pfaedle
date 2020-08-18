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

include_directories(pfaedle_dep PUBLIC ${PROJECT_SOURCE_DIR}/src/cppgtfs/src)
if (OSMIUM_FOUND)
	include_directories(${OSMIUM_INCLUDE_DIRS})
	target_link_libraries(pfaedle pfaedle_dep util configparser ad_cppgtfs -lpthread
			${OSMIUM_IO_LIBRARIES})
	message("linking with osmium")
else ()
	target_link_libraries(pfaedle pfaedle_dep util configparser ad_cppgtfs -lpthread)
endif()


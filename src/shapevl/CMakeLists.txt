file(GLOB_RECURSE shapevl_SRC *.cpp)

set(shapevl_main ShapevlMain.cpp)

list(REMOVE_ITEM shapevl_SRC ${shapevl_main})

include_directories(
	${PFAEDLE_INCLUDE_DIR}
	SYSTEM ${LIBZIP_INCLUDE_DIR}
	SYSTEM ${LIBZIP_CONF_INCLUDE_DIR}
)

add_executable(shapevl ${shapevl_main})
add_library(shapevl_dep ${shapevl_SRC})
target_link_libraries(shapevl_dep PUBLIC ad_cppgtfs)

include_directories(shapevl_dep PUBLIC ${PROJECT_SOURCE_DIR}/src/cppgtfs/src)
target_link_libraries(shapevl shapevl_dep util -lpthread ${LIBZIP_LIBRARY})

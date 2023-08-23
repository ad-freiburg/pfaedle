set(PFAEDLE_INCLUDE_DIR ${CMAKE_CURRENT_LIST_DIR} ${CMAKE_CURRENT_BINARY_DIR})

include_directories(
	${PFAEDLE_INCLUDE_DIR}
)

if (ZLIB_FOUND)
	include_directories( ${ZLIB_INCLUDE_DIRS} )
endif( ZLIB_FOUND )

add_subdirectory(util)
add_subdirectory(pfaedle)
add_subdirectory(cppgtfs)
add_subdirectory(configparser)
add_subdirectory(shapevl)
add_subdirectory(xml)

set(PFAEDLE_INCLUDE_DIR ${CMAKE_CURRENT_LIST_DIR} ${CMAKE_CURRENT_BINARY_DIR})

include_directories(
	${PFAEDLE_INCLUDE_DIR}
)

if (ZLIB_FOUND)
	include_directories( ${ZLIB_INCLUDE_DIRS} )
endif( ZLIB_FOUND )

if (LIBZIP_FOUND)
	include_directories( SYSTEM ${LIBZIP_INCLUDE_DIR} )
	include_directories( SYSTEM ${LIBZIP_CONF_INCLUDE_DIR} )
endif( LIBZIP_FOUND )

if (BZIP2_FOUND)
	include_directories( SYSTEM ${BZIP2_INCLUDE_DIR} )
endif( BZIP2_FOUND )

add_subdirectory(util)
add_subdirectory(pfaedle)
add_subdirectory(cppgtfs)
add_subdirectory(configparser)
add_subdirectory(shapevl)
add_subdirectory(xml)

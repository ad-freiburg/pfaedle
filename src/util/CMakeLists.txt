file(GLOB_RECURSE util_SRC *.cpp)
list(REMOVE_ITEM util_SRC TestMain.cpp)
add_library(util ${util_SRC})

if (!ZLIB_FOUND)
    find_package( ZLIB )
    if (ZLIB_FOUND)
		add_definitions( -DZLIB_FOUND=${ZLIB_FOUND} )
	endif()
endif()

if (ZLIB_FOUND)
	include_directories( ${ZLIB_INCLUDE_DIRS} )
	target_link_libraries( util ${ZLIB_LIBRARIES} )
endif( ZLIB_FOUND )

add_subdirectory(tests)

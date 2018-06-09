set(PFAEDLE_INCLUDE_DIR ${CMAKE_CURRENT_LIST_DIR} ${CMAKE_CURRENT_BINARY_DIR})

cpplint_add_subdirectory(pfaedle)

include_directories(
	${PFAEDLE_INCLUDE_DIR}
)

add_subdirectory(util)
add_subdirectory(pfaedle)
add_subdirectory(cppgtfs)
add_subdirectory(xml)
add_subdirectory(configparser)

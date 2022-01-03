set(PFAEDLE_INCLUDE_DIR ${CMAKE_CURRENT_LIST_DIR} ${CMAKE_CURRENT_BINARY_DIR})

if (COMMAND cpplint_add_subdirectory)
	cpplint_add_subdirectory(pfaedle)
endif()

include_directories(
	${PFAEDLE_INCLUDE_DIR}
)

add_subdirectory(util)
add_subdirectory(pfaedle)
add_subdirectory(cppgtfs)
add_subdirectory(configparser)
add_subdirectory(shapevl)

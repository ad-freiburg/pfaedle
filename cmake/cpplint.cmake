#
# CMake module to C++ static analysis against 
# Google C++ Style Guide (https://google.github.io/styleguide/cppguide.html)
#
# For more detials please follow links:
#
# - https://github.com/google/styleguide
# - https://pypi.python.org/pypi/cpplint
# - https://github.com/theandrewdavis/cpplint
#
# Copyright (c) 2016 Piotr L. Figlarek
# 
# Usage
# -----
# Include this module via CMake include(...) command and then add each source directory 
# via introduced by this module cpplint_add_subdirectory(...) function. Added directory 
# will be recursivelly scanned and all available files will be checked. 
#
# Example
# -------
# # include CMake module
# include(cmake/cpplint.cmake)
#
# # add all source code directories
# cpplint_add_subdirectory(core)
# cpplint_add_subdirectory(modules/c-bind)
#
# License (MIT)
# -------------
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


# select files extensions to check
option(CPPLINT_TEST_C_FILES     "Check *.c files"     ON)
option(CPPLINT_TEST_H_FILES     "Check *.h files"     ON)
option(CPPLINT_TEST_CPP_FILES   "Check *.cpp files"   ON)
option(CPPLINT_TEST_HPP_FILES   "Check *.hpp files"   ON)
option(CPPLINT_TEST_TPP_FILES   "Check *.tpp files"   ON)

# target to run cpplint.py for all configured sources
set(CPPLINT_TARGET lint CACHE STRING "Name of C++ style checker target") 

# project root directory
set(CPPLINT_PROJECT_ROOT "${PROJECT_SOURCE_DIR}" CACHE STRING "Project ROOT directory")


# find cpplint.py script
if(CPPLINT)
    message(STATUS "cpplint parser: ${CPPLINT}")
else()
    message(FATAL_ERROR "cpplint script: NOT FOUND! "
		"Please set the CPPLINT variable.")
endif()


# common target to concatenate all cpplint.py targets
add_custom_target(${CPPLINT_TARGET} ALL)


# use cpplint.py to check source code files inside DIR directory
function(cpplint_add_subdirectory DIR)
    # create relative path to the directory
    set(ABSOLUTE_DIR ${CMAKE_CURRENT_LIST_DIR}/${DIR})

    # add *.c files
    if(CPPLINT_TEST_C_FILES)
        set(EXTENSIONS       ${EXTENSIONS}c,)
        set(FILES_TO_CHECK   ${FILES_TO_CHECK} ${ABSOLUTE_DIR}/*.c)
    endif()

    # add *.h files
    if(CPPLINT_TEST_H_FILES)
        set(EXTENSIONS       ${EXTENSIONS}h,)
        set(FILES_TO_CHECK   ${FILES_TO_CHECK} ${ABSOLUTE_DIR}/*.h)
    endif()

    # add *.cpp files
    if(CPPLINT_TEST_CPP_FILES)
        set(EXTENSIONS       ${EXTENSIONS}cpp,)
        set(FILES_TO_CHECK   ${FILES_TO_CHECK} ${ABSOLUTE_DIR}/*.cpp)
    endif()

    # add *.hpp files
    if(CPPLINT_TEST_HPP_FILES)
        set(EXTENSIONS       ${EXTENSIONS}hpp,)
        set(FILES_TO_CHECK   ${FILES_TO_CHECK} ${ABSOLUTE_DIR}/*.hpp)
    endif()

    # add *.tpp files
    if(CPPLINT_TEST_TPP_FILES)
        set(EXTENSIONS       ${EXTENSIONS}tpp,)
        set(FILES_TO_CHECK   ${FILES_TO_CHECK} ${ABSOLUTE_DIR}/*.tpp)
    endif()

    # find all source files inside project
    file(GLOB_RECURSE LIST_OF_FILES ${FILES_TO_CHECK})

    # create valid target name for this check
    string(REGEX REPLACE "/" "." TEST_NAME ${DIR})
    set(TARGET_NAME ${CPPLINT_TARGET}.${TEST_NAME})

    # perform cpplint check
    add_custom_target(${TARGET_NAME}
        COMMAND ${CPPLINT} "--extensions=${EXTENSIONS}"
                           "--root=${CPPLINT_PROJECT_ROOT}"
						   "--quiet"
                           ${LIST_OF_FILES}
        DEPENDS ${LIST_OF_FILES}
        COMMENT "cpplint: Checking source code style"
    )

    # run this target when root cpplint.py test is triggered
    add_dependencies(${CPPLINT_TARGET} ${TARGET_NAME})

    # add this test to CTest
    add_test(${TARGET_NAME} ${CMAKE_MAKE_PROGRAM} ${TARGET_NAME})
endfunction()

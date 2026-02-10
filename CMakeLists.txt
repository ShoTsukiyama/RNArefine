cmake_minimum_required(VERSION 3.15)
project(RNArefine LANGUAGES CXX)

if(APPLE)
    # Apple Silicon (arm64) 用にビルドすることを強制する
    set(CMAKE_OSX_ARCHITECTURES "arm64")

    # OpenMPの場所を直接指定
    set(OpenMP_omp_LIBRARY "/opt/homebrew/opt/libomp/lib/libomp.dylib")
    set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include")
    set(OpenMP_CXX_LIB_NAMES "omp")
endif()
find_package(OpenMP REQUIRED)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(OpenMP REQUIRED)

file(GLOB_RECURSE SRC_FILES "src/*.cpp")

add_executable(RNArefine ${SRC_FILES})

file(GLOB_RECURSE INC_DIRS LIST_DIRECTORIES true "include/*")
foreach(dir ${INC_DIRS})
    if(IS_DIRECTORY ${dir})
        target_include_directories(RNArefine PRIVATE ${dir})
    endif()
endforeach()
target_include_directories(RNArefine PRIVATE include include/MC)

target_link_libraries(RNArefine PRIVATE OpenMP::OpenMP_CXX)

set_target_properties(RNArefine PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

install(TARGETS RNArefine DESTINATION bin)
install(DIRECTORY init_param/ DESTINATION share/RNArefine/init_param)

target_compile_definitions(RNArefine PRIVATE 
    RNAREFINE_DATA_PATH="${CMAKE_INSTALL_PREFIX}/share/RNArefine"
)

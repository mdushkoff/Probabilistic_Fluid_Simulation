cmake_minimum_required(VERSION 3.17.0)

# Version Number
set(PROJECT_MAJOR_VERSION 0)
set(PROJECT_MINOR_VERSION 0)
set(PROJECT_DEBUG_VERSION 1)
set(PROJECT_VERSION ${PROJECT_MAJOR_VERSION}.${PROJECT_MINOR_VERSION}.${PROJECT_DEBUG_VERSION})

# Project Name
project(FluidSim VERSION ${PROJECT_VERSION} LANGUAGES C CXX)
string(TOLOWER ${PROJECT_NAME} EXECUTABLE_NAME)
set(EXECUTABLE_NAME_CPU "${EXECUTABLE_NAME}_cpu")

# Packages
find_package(PNG REQUIRED)
find_library(MATH_LIBRARY m)

# C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Options
option(USE_CUDA "Whether or not to use CUDA" ON)
if(USE_CUDA)
    find_package(CUDAToolkit REQUIRED)
    set(EXECUTABLE_NAME_CUDA "${EXECUTABLE_NAME}_cuda")
endif()

# Sources
set(COMMON_SOURCE_FILES
    "src/main.cpp"
)
set(CPU_SOURCE_FILES
    "src/fluid.cpp"
)
set(GPU_SOURCE_FILES
    "src/fluid.cu"
)

# CPU Targets
add_executable(${EXECUTABLE_NAME_CPU}
    "${COMMON_SOURCE_FILES}"
    "${CPU_SOURCE_FILES}"
)
target_include_directories(${EXECUTABLE_NAME_CPU} PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/includes"
)
target_link_libraries(${EXECUTABLE_NAME_CPU}
    ${PNG_LIBRARY}
    ${MATH_LIBRARY}
)

# GPU Targets
if(USE_CUDA)
    enable_language("CUDA")
    add_executable(${EXECUTABLE_NAME_CUDA}
        "${COMMON_SOURCE_FILES}"
        "${GPU_SOURCE_FILES}"
    )
    set_target_properties(${EXECUTABLE_NAME_CUDA} PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
    )
    set_target_properties(${EXECUTABLE_NAME_CUDA} PROPERTIES
        CUDA_ARCHITECTURES "75"
    )
    target_include_directories(${EXECUTABLE_NAME_CUDA} PRIVATE
        "${CUDAToolkit_INCLUDE_DIRS}"
        "${CMAKE_CURRENT_SOURCE_DIR}/includes"
    )
    target_compile_definitions(${EXECUTABLE_NAME_CUDA} PRIVATE
        -DUSE_CUDA
    )
    target_link_libraries(${EXECUTABLE_NAME_CUDA} PUBLIC
        CUDA::cudart
        ${PNG_LIBRARY}
        ${MATH_LIBRARY}
    )
endif()

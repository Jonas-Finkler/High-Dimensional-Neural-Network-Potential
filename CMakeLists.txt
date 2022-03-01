cmake_minimum_required(VERSION 3.6)
project(hdnnp Fortran)
enable_language(Fortran)


option(INTEL "compile with ifort and mkl" OFF)
option(OPENMP "compile with openmp support" OFF)
option(DEBUG "compile with debug flags" ON)

if(INTEL)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mkl -no-wrap-margin")
    if(DEBUG)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -traceback -O0 -fpe0 ") # -Ofast -ipo -check all
    else()
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -xHost -Ofast -ipo")
    endif()
    if(OPENMP)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qopenmp")
    endif()
else()
    if(DEBUG)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Og -fcheck=all -Wall -Wconversion -g -fbacktrace -ffpe-trap=zero,invalid,overflow -flto")
    else()
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Ofast -march=native -mtune=native -mfma -fomit-frame-pointer -ffast-math -flto")
    endif()
    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)
    set(libraries ${libraries} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
    if(OPENMP)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp")
    endif()
endif()


add_executable(
        hdnnp.x


        ${sources}

        src/potentials/hdnnpExample.f90 # <- main

        src/potentials/hdnnps.f90
        src/potentials/symmetryFunctions.f90
        src/potentials/neighborLists.f90

        src/neuralnetwork/neuralNetworks.f90

        src/math/linalg.f90
        src/math/quaternions.f90
        src/math/random.f90
        src/math/mathUtil.f90

        src/util/doubleLinkedList.f90
        src/util/util.f90
        src/util/dynamicArray.f90
        src/util/sorting.f90
        #src/util/argumentParser.f90

        src/constants/constants.f90

        src/configuration/parameters.f90

        src/atoms/atomicStructure.f90
        src/atoms/pointsets.f90

        src/precision.f90
)

target_link_libraries(hdnnp.x ${libraries})

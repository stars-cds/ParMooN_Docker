# ===================================================================
# This is a user configuration file for ParMooN Version 1.1
# written by Sashikumaar Ganesan, CDS, IISc Bangalore, India
# date: 05 June 2015
# ===================================================================

# controlling the output messages
set(CMAKE_VERBOSE_MAKEFILE FALSE)

# selection of dimension (2D 3D)
# set(AParMooN_GEO "3D" CACHE STRING "Change AParMooN_GEO, to select the Dimensio of the problem")
set(AParMooN_GEO "2D" CACHE STRING "Change AParMooN_GEO, to select the Dimension of the problem")

#...................................................................................................................................................
# select this line accordingly to include your main program
set(AParMooN_MODEL "${PROJECT_SOURCE_DIR}/ProblemSets/Scalar_CD2D.C" CACHE STRING "Enter to select the Main file of the model") 
# set(AParMooN_MODEL "${PROJECT_SOURCE_DIR}/2DPrograms/CD2D_ParMooN.C" CACHE STRING "Enter to select the Main file of the model") 
# set(AParMooN_MODEL "${PROJECT_SOURCE_DIR}/2DPrograms/TCD2D_ParMooN.C" CACHE STRING "Enter to select the Main file of the model") 
# set(AParMooN_MODEL "${PROJECT_SOURCE_DIR}/2DPrograms/TCD2D_ParMooN_ALE.C" CACHE STRING "Enter to select the Main file of the model") 
# set(AParMooN_MODEL "${PROJECT_SOURCE_DIR}/2DPrograms/NSE2D_ParMooN.C" CACHE STRING "Enter to select the Main file of the model")
# set(AParMooN_MODEL "${PROJECT_SOURCE_DIR}/2DPrograms/TNSE2D_ParMooN.C" CACHE STRING "Enter to select the Main file of the model")
# set(AParMooN_MODEL "${PROJECT_SOURCE_DIR}/3DPrograms/CD3D_ParMooN.C" CACHE STRING "Enter to select the Main file of the model")
# set(AParMooN_MODEL "${PROJECT_SOURCE_DIR}/3DPrograms/TCD3D_ParMooN.C" CACHE STRING "Enter to select the Main file of the model")
# set(AParMooN_MODEL "${PROJECT_SOURCE_DIR}/3DPrograms/NSE3D_ParMooN.C" CACHE STRING "Enter to select the Main file of the model")
# set(AParMooN_MODEL "${PROJECT_SOURCE_DIR}/3DPrograms/TNSE3D_ParMooN.C" CACHE STRING "Enter to select the Main file of the model")
# set(AParMooN_MODEL "${PROJECT_SOURCE_DIR}/Main_Users/Thivin/2D_Programs/deepikaji.cpp" CACHE STRING "Enter to select the Main file of the model")
# set(AParMooN_MODEL "${PROJECT_SOURCE_DIR}/Main_Users/Thivin/TNSE3D/thivin_TNSE3D.cpp" CACHE STRING "Enter to select the Main file of the model")
# set(AParMooN_MODEL "${PROJECT_SOURCE_DIR}/Main_Users/Thivin/Sample_mesh_move_2d.cpp" CACHE STRING "Enter to select the Main file of the model")

# set(AParMooN_MODEL "${PROJECT_SOURCE_DIR}/Assignment2/Prob2.C" CACHE STRING "Enter to select the Main file of the model") 

# selection of architect type (LINUX64 MAC64 INTEL64 TYRONE64 CRAY64)
set(AParMooN_ARCH "LINUX64" CACHE STRING "select the machine type")

#  selection of program type (SEQUENTIAL SMPI MPI OMPONLY HYBRID SCUDA)
set(AParMooN_PARALLEL_TYPE "SEQUENTIAL" CACHE STRING "select the parallel type")

#  selection of program type (MPICH OPENMPI INTELMPI CRAYMPI MACMPI)
set(AParMooN_MPI_IMPLEMENTATION "INTELMPI" CACHE STRING "select the MPI Implementation type")
 
# set the path to save the exe file ....................................................................................
#.......................................................................................................................
# set(AParMooN_OUTPUT_DIR_PATH "${CMAKE_SOURCE_DIR}/../ParMooN_Output/aletnse3d" CACHE STRING "select the model")
 set(AParMooN_OUTPUT_DIR_PATH "${CMAKE_SOURCE_DIR}/../ParMooN_Output/CD2D" CACHE STRING "select the model")
# set(AParMooN_OUTPUT_DIR_PATH "${CMAKE_SOURCE_DIR}/../ParMooN_Output/NSE2D" CACHE STRING "select the model")
# set(AParMooN_OUTPUT_DIR_PATH "${CMAKE_SOURCE_DIR}/../ParMooN_Output/TNSE2D" CACHE STRING "select the model")
# set(AParMooN_OUTPUT_DIR_PATH "${CMAKE_SOURCE_DIR}/../ParMooN_Output/CD3D" CACHE STRING "select the model")
# set(AParMooN_OUTPUT_DIR_PATH "${CMAKE_SOURCE_DIR}/../ParMooN_Output/TCD3D" CACHE STRING "select the model")
# set(AParMooN_OUTPUT_DIR_PATH "${CMAKE_SOURCE_DIR}/../ParMooN_Output/NSE3D" CACHE STRING "select the model")
# set(AParMooN_OUTPUT_DIR_PATH "${CMAKE_SOURCE_DIR}/../ParMooN_Output/TNSE3D" CACHE STRING "select the model")
# set(AParMooN_OUTPUT_DIR_PATH "${CMAKE_SOURCE_DIR}/../ParMooN_Output/thivin3d" CACHE STRING "select the model")
# set(AParMooN_OUTPUT_DIR_PATH "${CMAKE_SOURCE_DIR}/../ParMooN_Output/mesh2d" CACHE STRING "select the model")
# set(AParMooN_OUTPUT_DIR_PATH "${CMAKE_SOURCE_DIR}/../ParMooN_Output/deepikaji" CACHE STRING "select the model")

# set(AParMooN_OUTPUT_DIR_PATH "${CMAKE_SOURCE_DIR}/../ParMooN_Output/Assignment2" CACHE STRING "select the model")


set(USE_PARMOON_DEFINE -D__PRIVATE__)

# CMAKE_BUILD_TYPE [ DEBUG | RELEASE | RELWITHDEBINFO | MINSIZEREL ]
set(EXE_BUILD_TYPE RELEASE)


# set FALSE, if you want to use libs provided in PARMOON 
# if you set TRUE, it will search in all you lib paths and if not found, PARMOON libs will be used
set(AParMooN_USE_SYSTEM_MKLBLAS TRUE)
set(AParMooN_USE_SYSTEM_UMFPACK TRUE)
set(AParMooN_USE_SYSTEM_LAPACK TRUE)
set(AParMooN_USE_SYSTEM_MUMPS TRUE)
set(AParMooN_USE_SYSTEM_GRIDGEN TRUE)
set(AParMooN_USE_SYSTEM_TETGEN TRUE)

# set(USE_PARMOON_DEFINE -D__PRIVATE__)
# SET(DAI_WITH_BP ON CACHE BOOL "Belief Propagation" FORCE)
# ========================================================================================================================
# no need to change anyting after this line
# used only when ccmake or cmake-gui is used
# ========================================================================================================================
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${AParMooN_OUTPUT_DIR_PATH}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${AParMooN_OUTPUT_DIR_PATH}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${AParMooN_OUTPUT_DIR_PATH})

set_property(CACHE AParMooN_GEO PROPERTY STRINGS 2D 3D ) 

# selection of all main programs
if("${AParMooN_GEO}" STREQUAL "2D")
  file(GLOB_RECURSE MAIN_SOURCES "${PROJECT_SOURCE_DIR}/2DPrograms/*.C")
  set_property(CACHE AParMooN_MODEL PROPERTY STRINGS  ${MAIN_SOURCES})   
elseif("${AParMooN_GEO}" STREQUAL "3D")
  file(GLOB_RECURSE MAIN_SOURCES "${PROJECT_SOURCE_DIR}/3DPrograms/*.C")
  set_property(CACHE AParMooN_MODEL PROPERTY STRINGS  ${MAIN_SOURCES})  
endif()  

# selection of all architects
set_property(CACHE AParMooN_ARCH PROPERTY STRINGS LINUX64 MAC64 INTEL64 TYRONE64 CRAY64)

# selection of all program types
set_property(CACHE AParMooN_PARALLEL_TYPE PROPERTY STRINGS SEQUENTIAL SMPI MPI OMPONLY HYBRID SCUDA)

# selection of all program types
set_property(CACHE AParMooN_MPI_IMPLEMENTATION PROPERTY STRINGS MPICH OPENMPI INTELMPI CRAYMPI MACMPI)

# selection of all output folder
if("${AParMooN_GEO}" STREQUAL "2D")
set_property(CACHE AParMooN_OUTPUT_DIR_PATH PROPERTY STRINGS  ${CMAKE_SOURCE_DIR}/OutPut/cd2d  ${CMAKE_SOURCE_DIR}/OutPut/tcd2d  ${CMAKE_SOURCE_DIR}/OutPut/nse2d  ${CMAKE_SOURCE_DIR}/OutPut/tnse2d)
 elseif("${AParMooN_GEO}" STREQUAL "3D")
set_property(CACHE AParMooN_OUTPUT_DIR_PATH PROPERTY STRINGS  ${CMAKE_SOURCE_DIR}/OutPut/cd3d  ${CMAKE_SOURCE_DIR}/OutPut/tcd3d  ${CMAKE_SOURCE_DIR}/OutPut/nse3d  ${CMAKE_SOURCE_DIR}/OutPut/tnse3d)
endif()

# ======================================================================
# general settings
# ======================================================================
 if("${AParMooN_PARALLEL_TYPE}" STREQUAL "MPI")
   set(PARMOON_PRG_DEFINE "-D_PAR -D_MPIONLY -D_MPI")
 elseif("${AParMooN_PARALLEL_TYPE}" STREQUAL "OMPONLY")
   set(PARMOON_PRG_DEFINE "-D_PAR -D_OMPONLY -D_OMP")
 elseif("${AParMooN_PARALLEL_TYPE}" STREQUAL "HYBRID")
   set(PARMOON_PRG_DEFINE "-D_PAR -D_MPI -D_HYBRID")
 elseif("${AParMooN_PARALLEL_TYPE}" STREQUAL "SEQUENTIAL")
    set(PARMOON_PRG_DEFINE "-D_SEQ ")
 elseif("${AParMooN_PARALLEL_TYPE}" STREQUAL "SCUDA")
    set(PARMOON_PRG_DEFINE "-D_SEQ -D_CUDA")
  elseif("${AParMooN_PARALLEL_TYPE}" STREQUAL "SMPI")
    set(PARMOON_PRG_DEFINE "-D_SEQ -D_SMPI")
 endif()

 if("${AParMooN_ARCH}" STREQUAL "LINUX64")
   set(PARMOON_CXX_DEF "${PARMOON_CXX_DEF} -fopenmp -std=c++11 ")
   set(PARMOON_C_DEF "  ${PARMOONx_C_DEF}  -DREDUCED -DNO_TIMER -m64   ")
   set(PARMOON_CUDA_DEF "${PARMOON_CUDA_DEF} -arch=sm_70 -Xcompiler -fopenmp -w ")
   
 elseif("${AParMooN_ARCH}" STREQUAL "MAC64")
   set(PARMOON_CXX_DEF " -stdlib=libc++ ${PARMOON_CXX_DEF}  -fapple-pragma-pack -Wdeprecated-register  ")
   set(PARMOON_C_DEF "  ${PARMOON_C_DEF}  -DREDUCED -DNO_TIMER -DMKL_ILP64 -m64 -fapple-pragma-pack ")
 elseif("${AParMooN_ARCH}" STREQUAL "INTEL64")
   set(PARMOON_CXX_DEF "${PARMOON_CXX_DEF} -std=c++11   -qopenmp ")
   set(PARMOON_C_DEF "  ${PARMOON_C_DEF}   -DREDUCED -DNO_TIMER -DMKL_ILP64 -m64  ")
   set(PARMOON_CUDA_DEF "${PARMOON_CUDA_DEF} -arch=sm_70 -Xcompiler -fopenmp -w ")

 elseif("${AParMooN_ARCH}" STREQUAL "TYRONE64")
   set(PARMOON_CXX_DEF "${PARMOON_CXX_DEF} -DREDUCED -DNO_TIMER")
   set(PARMOON_C_DEF "  ${PARMOON_C_DEF}  -DREDUCED -DNO_TIMER  -DMPICH_IGNORE_CXX_SEEK ")  
 elseif("${AParMooN_ARCH}" STREQUAL "CRAY64")
   set(PARMOON_CXX_DEF "${PARMOON_CXX_DEF} -DMPICH_IGNORE_CXX_SEEK ")  
   set(PARMOON_C_DEF "  ${PARMOON_C_DEF}  -DREDUCED -DNO_TIMER -DMPICH_IGNORE_CXX_SEEK  ")
   set(PARMOON_CUDA_DEF "${PARMOON_CUDA_DEF} -arch=sm_70 -Xcompiler -fopenmp -w ")

 endif()

 #### SET UP MKL BLAS VARIABLES ######
if ("${MKLBLASFOUND}" STREQUAL "TRUE")
if("${AParMooN_ARCH}" STREQUAL "LINUX64")  
  set(PARMOON_CXX_DEF "${PARMOON_CXX_DEF} -lmkl_rt -lmkl_intel_thread -lpthread -liomp5 ")
  set(PARMOON_C_DEF "  ${PARMOONx_C_DEF}  -lmkl_rt -lmkl_intel_thread")
endif("${AParMooN_ARCH}" STREQUAL "LINUX64") 

if("${AParMooN_ARCH}" STREQUAL "INTEL64")  
  set(PARMOON_CXX_DEF "${PARMOON_CXX_DEF} -lmkl_rt -lmkl_intel_thread")
  set(PARMOON_C_DEF "  ${PARMOONx_C_DEF}  -lmkl_rt -lmkl_intel_thread")
endif("${AParMooN_ARCH}" STREQUAL "INTEL64") 

endif ("${MKLBLASFOUND}" STREQUAL "TRUE")
####  --END-- SET UP MKL BLAS VARIABLES ######
 
set(PARMOON_C_DEF " ${PARMOON_C_DEF} -D__${AParMooN_GEO}__ -D__${AParMooN_ARCH}__  -DTRILIBRARY -DTETLIBRARY  ${PARMOON_PRG_DEFINE} ${USE_PARMOON_DEFINE}") 
set(PARMOON_CXX_DEF " ${PARMOON_CXX_DEF} -D__${AParMooN_GEO}__ -D__${AParMooN_ARCH}__  -DTRILIBRARY -DTETLIBRARY ${PARMOON_PRG_DEFINE} ${USE_PARMOON_DEFINE}")
set(PARMOON_CUDA_DEF " ${PARMOON_CUDA_DEF} -D__${AParMooN_GEO}__ -D__${AParMooN_ARCH}__  -DTRILIBRARY -DTETLIBRARY ${PARMOON_PRG_DEFINE} ${USE_PARMOON_DEFINE}")

cmake_minimum_required( VERSION 2.8.12 )

project( MAT_IncPions )

#Add extra "modules" to CMake
include(CTest)
include(AddFileDependencies)

#Compiler flags
set( GCC_Flags_For_CXX "-std=c++1y -Wall" )
set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_Flags_For_CXX}" )
set( CMAKE_CXX_FLAGS_DEBUG "-ggdb" )
set( CMAKE_CXX_FLAGS_RELEASE "-O2 -DNDEBUG" )

#Add a PROFiling CMAKE_BUILD_TYPE
set( CMAKE_CXX_FLAGS_PROF "-ggdb -pg -DNDEBUG" )
set( CMAKE_C_FLAGS_PROF "-ggdb -pg -DNDEBUG" )
set( CMAKE_EXE_LINKER_FLAGS_PROF "-ggdb -pg -DNDEBUG" )
set( CMAKE_SHARED_LINKER_FLAGS_PROF "-ggdb -pg -DNDEBUG" )

#Tell this package where it is installed and version control status
add_definitions(-DINSTALL_DIR="${CMAKE_INSTALL_PREFIX}/")

#Let directories in this package see each other
include_directories( "${PROJECT_SOURCE_DIR}" )

#Find dependencies
list( APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS RIO Net Tree Core Geom EG GenVector Minuit2 Minuit OPTIONAL_COMPONENTS Cintex)

if(${ROOT_VERSION} VERSION_LESS 6 AND NOT ${ROOT_Cintex_FOUND})
  MESSAGE(FATAL_ERROR "Cintex is optional except when it's not.  ROOT 6 has Reflex "
                      "support built in, so it doesn't have a separate component "
                      "for Cintex.  Cintex was an experimental feature in ROOT 5, "
                      "so I have to require it as a component there.  You appear to "
                      "be using ${ROOT_VERSION}, so I can't let you get away with "
                      "skipping out on Cintex support!  I need Cintex to load metadata "
                      "about MnvH1D and MnvH2D at runtime.")
endif()

if(${ROOT_VERSION} VERSION_LESS 6)
  MESSAGE("I need Cintex on this platform, so making sure to compile against it...")
else()
  MESSAGE("I don't need Cintex on this platform, so skipping it...")
  add_definitions(-DNCINTEX)
endif()

include(${ROOT_USE_FILE})

find_package(yaml-cpp 0.6.0 REQUIRED)
include_directories(${YAML_CPP_INCLUDE_DIR})

find_package(MAT REQUIRED)
include_directories(${MAT_INCLUDE_DIR})

find_package(MAT-MINERvA REQUIRED)
include_directories(${MAT-MINERvA_INCLUDE_DIR})

find_package(UnfoldUtils REQUIRED)
include_directories(${UnfoldUtils_INCLUDE_DIR})

find_package(GENIEXSecExtract REQUIRED)
include_directories(${GENIEXSecExtract_INCLUDE_DIR})

#Tell CMake to compile subdirectories before compiling main application
add_subdirectory(playlists)
add_subdirectory(util)
add_subdirectory(event)
add_subdirectory(fits)
#TODO: Split these directories' headers into .cpp files
#add_subdirectory(cuts)
#add_subdirectory(studies)
#add_subdirectory(systematics)
#add_subdirectory(plotTools)

#Build main executables
add_executable(runEventLoop runEventLoop.cpp)
target_link_libraries(runEventLoop ${ROOT_LIBRARIES} util event fits MAT MAT-MINERvA yaml-cpp) #event cuts studies systematics)
install(TARGETS runEventLoop DESTINATION bin)

add_executable(ExtractCrossSection ExtractCrossSection.cpp)
target_link_libraries(ExtractCrossSection ${ROOT_LIBRARIES} util event fits MAT UnfoldUtils)
install(TARGETS ExtractCrossSection DESTINATION bin)

add_executable(ExtractCrossSection2D ExtractCrossSection2D.cpp)
target_link_libraries(ExtractCrossSection2D ${ROOT_LIBRARIES} util event fits MAT UnfoldUtils)
install(TARGETS ExtractCrossSection2D DESTINATION bin)

add_executable(runXSecLooper runXSecLooper.cpp)
target_link_libraries(runXSecLooper ${ROOT_LIBRARIES} MAT GENIEXSecExtract)
install(TARGETS runXSecLooper DESTINATION bin)

add_executable(signalBKGStack signalBKGStack.cxx)
target_link_libraries(signalBKGStack ${ROOT_LIBRARIES} MAT MAT-MINERvA)# plotTools)
install(TARGETS signalBKGStack DESTINATION bin)

add_executable(Test2DPlots Test2DPlots.cxx)
target_link_libraries(Test2DPlots ${ROOT_LIBRARIES} MAT MAT-MINERvA)# plotTools)
install(TARGETS Test2DPlots DESTINATION bin)

add_executable(CombineTgts CombineTgts.cxx)
target_link_libraries(CombineTgts ${ROOT_LIBRARIES} MAT MAT-MINERvA)
install(TARGETS CombineTgts DESTINATION bin)

add_executable(ScaleHistos ScaleHistos.cxx)
target_link_libraries(ScaleHistos ${ROOT_LIBRARIES} MAT MAT-MINERvA)
install(TARGETS ScaleHistos DESTINATION bin)

add_executable(ScaleErrorBand ScaleErrorBand.cxx)
target_link_libraries(ScaleErrorBand ${ROOT_LIBRARIES} MAT MAT-MINERvA)
install(TARGETS ScaleErrorBand DESTINATION bin)

add_executable(ScaleMCToData ScaleMCToData.cxx)
target_link_libraries(ScaleMCToData ${ROOT_LIBRARIES} MAT MAT-MINERvA)
install(TARGETS ScaleMCToData DESTINATION bin)

add_executable(ProjectOuterPlastic ProjectOuterPlastic.cxx)
target_link_libraries(ProjectOuterPlastic ${ROOT_LIBRARIES} MAT MAT-MINERvA)
install(TARGETS ProjectOuterPlastic DESTINATION bin)

add_executable(AddUpPOT AddUpPOT.cxx)
target_link_libraries(AddUpPOT ${ROOT_LIBRARIES} MAT MAT-MINERvA)
install(TARGETS AddUpPOT DESTINATION bin)

add_executable(CombineFits2D CombineFits2D.cxx)
target_link_libraries(CombineFits2D ${ROOT_LIBRARIES} MAT MAT-MINERvA)
install(TARGETS CombineFits2D DESTINATION bin)

add_executable(PlotNPlanesByTgt PlotNPlanesByTgt.cxx)
target_link_libraries(PlotNPlanesByTgt ${ROOT_LIBRARIES} MAT MAT-MINERvA)
install(TARGETS PlotNPlanesByTgt DESTINATION bin)

add_executable(FitNPlanes FitNPlanes.cxx)
target_link_libraries(FitNPlanes ${ROOT_LIBRARIES} fits MAT MAT-MINERvA)
install(TARGETS FitNPlanes DESTINATION bin)

add_executable(FitNPlanes2D FitNPlanes2D.cxx)
target_link_libraries(FitNPlanes2D ${ROOT_LIBRARIES} fits MAT MAT-MINERvA)
install(TARGETS FitNPlanes2D DESTINATION bin)

add_executable(twoPlotComparison twoPlotComparison.cxx)
target_link_libraries(twoPlotComparison ${ROOT_LIBRARIES} MAT MAT-MINERvA)
install(TARGETS twoPlotComparison DESTINATION bin)

add_executable(BKGFitting BKGFitting.cxx)
target_link_libraries(BKGFitting ${ROOT_LIBRARIES} fits MAT MAT-MINERvA)
install(TARGETS BKGFitting DESTINATION bin)

add_executable(BKGFittingNEWSTRUCTURE BKGFittingNEWSTRUCTURE.cxx)
target_link_libraries(BKGFittingNEWSTRUCTURE ${ROOT_LIBRARIES} fits MAT MAT-MINERvA)
install(TARGETS BKGFittingNEWSTRUCTURE DESTINATION bin)

add_executable(BKGFittingOverhaul BKGFittingOverhaul.cxx)
target_link_libraries(BKGFittingOverhaul ${ROOT_LIBRARIES} fits MAT MAT-MINERvA yaml-cpp)
install(TARGETS BKGFittingOverhaul DESTINATION bin)

add_executable(BKGFitting2D BKGFitting2D.cxx)
target_link_libraries(BKGFitting2D ${ROOT_LIBRARIES} fits MAT MAT-MINERvA)
install(TARGETS BKGFitting2D DESTINATION bin)

add_executable(BKGFittingTFrac BKGFittingTFrac.cxx)
target_link_libraries(BKGFittingTFrac ${ROOT_LIBRARIES} MAT MAT-MINERvA)
install(TARGETS BKGFittingTFrac DESTINATION bin)

add_executable(TMinFitting TMinFitting.cxx)
target_link_libraries(TMinFitting ${ROOT_LIBRARIES} fits MAT MAT-MINERvA)
install(TARGETS TMinFitting DESTINATION bin)

install(FILES backgroundStack.py PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE DESTINATION bin)

install(FILES runTransWarp.sh PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE DESTINATION bin)

configure_file(setup.sh.in setup_${PROJECT_NAME}.sh @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/setup_${PROJECT_NAME}.sh DESTINATION bin)

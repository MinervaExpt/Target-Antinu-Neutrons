#add_library(event CVUniverse.cpp)
add_library(fits ScaleFactors.cpp Fit.cpp ScaleFactor.cpp Line.cpp NonFit.cpp Piecewise.cpp FitMgr.cpp)
target_link_libraries(fits ${ROOT_LIBRARIES})
install(TARGETS fits DESTINATION lib)

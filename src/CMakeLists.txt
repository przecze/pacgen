file(GLOB packgen_srcs "*.cpp")
file(GLOB packgen_headers "*.h")
add_library(packgen ${packgen_headers} ${packgen_srcs})
target_include_directories(packgen PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                                         $<INSTALL_INTERFACE:include/packgen>)
install(TARGETS packgen EXPORT packgen DESTINATION lib)
install(FILES ${packgen_headers} DESTINATION include/packgen)
install(FILES packgen-config.cmake DESTINATION lib/cmake/packgen)
install(EXPORT packgen DESTINATION lib/cmake/packgen)

FILE(GLOB src_list *.cpp *.h)

add_executable(GreenFunction ${src_list})


target_link_libraries(GreenFunction libhald.a libgfortran.a)





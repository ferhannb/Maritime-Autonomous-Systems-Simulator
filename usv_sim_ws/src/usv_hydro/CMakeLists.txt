cmake_minimum_required(VERSION 3.8)
project(usv_hydro)

find_package(ament_cmake REQUIRED)
find_package(gz-sim8 REQUIRED)
find_package(gz-plugin2 REQUIRED COMPONENTS register)

add_library(usv_hydro_system SHARED src/HydroCellSystem.cc)
target_compile_features(usv_hydro_system PUBLIC cxx_std_17)
target_link_libraries(usv_hydro_system
  PRIVATE
    gz-sim8::gz-sim8
    gz-plugin2::register
)

install(
  TARGETS usv_hydro_system
  LIBRARY DESTINATION lib
)

ament_package()

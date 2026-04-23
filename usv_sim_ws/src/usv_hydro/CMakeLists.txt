cmake_minimum_required(VERSION 3.8)
project(usv_hydro)

find_package(ament_cmake REQUIRED)
find_package(gz-sim8 REQUIRED)
find_package(gz-plugin2 REQUIRED COMPONENTS register)

add_library(usv_hydro_system SHARED
  src/BuoyancyModel.cc
  src/DragModel.cc
  src/EnvironmentModel.cc
  src/HydroCellSystem.cc
  src/HydroIntegrator.cc
  src/HydroTypes.cc
)
target_compile_features(usv_hydro_system PUBLIC cxx_std_17)
target_include_directories(usv_hydro_system
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)
target_link_libraries(usv_hydro_system
  PUBLIC
    gz-sim8::gz-sim8
  PRIVATE
    gz-plugin2::register
)

install(
  TARGETS usv_hydro_system
  LIBRARY DESTINATION lib
)

install(
  DIRECTORY include/
  DESTINATION include
)

if(BUILD_TESTING)
  find_package(ament_cmake_gtest REQUIRED)

  ament_add_gtest(hydro_models_test
    test/hydro_models_test.cc
  )
  target_link_libraries(hydro_models_test
    usv_hydro_system
  )
endif()

ament_package()

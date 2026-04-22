cmake_minimum_required(VERSION 3.8)
project(usv_description)

find_package(ament_cmake REQUIRED)

install(
  DIRECTORY launch sdf config meshes
  DESTINATION share/${PROJECT_NAME}
)

ament_package()

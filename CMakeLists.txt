cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(AEGIS)
add_subdirectory(MICROCERN)
add_subdirectory(THijing)
add_subdirectory(TDPMjet)

#Install(DIRECTORY Geometry gconfig DESTINATION share/Detectors/)

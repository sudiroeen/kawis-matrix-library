cmake_minimum_required(VERSION 2.8.3)
project(main)

add_compile_options(-std=c++14)

set(${HEADER}
   kawis/kawis.hpp
)

set(${SOURCE}
   src/Matrix.cpp
)

add_executable(${PROJECT_NAME}
   ${HEADER}
   ${SOURCE}
   main.cpp
)

cmake_minimum_required(VERSION 3.25)
project(E5_Curvas)

set(CMAKE_CXX_STANDARD 17)

add_executable(E5_Curvas main.cpp glad.c glut_ply.cpp Objeto.cpp BoundingVolume.h DTriangulation.h)

target_link_libraries(E5_Curvas -lOpenGL32 -lGLU32 -lfreeGLUT -lglfw3)

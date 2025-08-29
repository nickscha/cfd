@echo off

set DEF_FLAGS_COMPILER=-std=c89 -march=native -mtune=native -pedantic -ffast-math -funroll-loops -finline-functions -flto -Wall -Wextra -Werror -Wvla -Wconversion -Wdouble-promotion -Wsign-conversion -Wuninitialized -Winit-self -Wunused -Wunused-macros -Wunused-local-typedefs
set DEF_FLAGS_LINKER=
set SOURCE_NAME=lbm_3d

rm *.ppm
cc -s -O3 %DEF_FLAGS_COMPILER% -o %SOURCE_NAME%.exe %SOURCE_NAME%.c %DEF_FLAGS_LINKER%
%SOURCE_NAME%.exe

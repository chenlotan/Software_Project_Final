#!/bin/bash
python3 setup.py build_ext --inplace
gcc -ansi -Wall -Wextra -Werror -pedantic-errors spkmeans.c -lm -o spkmeans



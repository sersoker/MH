#! /bin/bash
g++ principal.cpp -std=c++11 libgtest.a libarff.a -o ../BIN/programa
cd ../BIN
./programa
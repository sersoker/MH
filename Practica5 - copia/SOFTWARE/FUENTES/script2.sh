#! /bin/bash
g++ principal2.cpp -std=c++11 libgtest.a libarff.a -o ../BIN/programa
cd ../BIN
./programa
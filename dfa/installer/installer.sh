#!/bin/bash

rm -rf build dfa dfa.c  dfa.cpython-36m-x86_64-linux-gnu.so
python setup.py build_ext --inplace
cp dfa/* . 
python tester.py

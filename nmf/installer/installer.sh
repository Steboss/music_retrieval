#!/bin/bash

rm -rf build nmf nmf.c  nmf.cpython-36m-x86_64-linux-gnu.so
python setup.py build_ext --inplace
cp nmf/* . 
python tester.py

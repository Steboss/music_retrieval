#!/bin/bash

rm -rf build mfcc mfcc.c  mfcc.cpython-36m-x86_64-linux-gnu.so
python setup.py build_ext --inplace
cp mfcc/* .
python tester.py

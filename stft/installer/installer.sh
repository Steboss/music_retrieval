#!/bin/bash

rm -rf build stft stft.c  stft.cpython*
python setup.py build_ext --inplace
cp nmf/* . 
python tester.py

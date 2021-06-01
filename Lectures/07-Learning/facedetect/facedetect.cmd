#!/bin/bash

build/facedetect --cascade="/usr/local/opencv-3.4.1/share/OpenCV/haarcascades/haarcascade_frontalface_alt.xml" --nested-cascade="/usr/local/opencv-3.4.1/share/OpenCV/haarcascades/haarcascade_eye.xml" --scale=1.3 0


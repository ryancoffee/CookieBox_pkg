#!/sdf/group/lcls/ds/ana/sw/conda1/inst/envs/ana-4.0.66-py3/bin/python
import numpy as np
import cv2

def massage(im):
    im = cv2.erode(im,kernel=(3,3))
    im = cv2.blur(im,ksize=(3,3))
    return im

def getshape(raw):
    nx = len(raw)
    ny = len(raw[0])
    #print(nx,ny)
    return nx,ny

def bit_count(x:int):
    return bin(x).count("1")


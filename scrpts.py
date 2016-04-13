import cv2
import numpy as np
from matplotlib import pyplot as plt


def imgChop(img, pos, shape):
    newimage=np.zeros(shape,img.dtype)
    for i in range(shape[0]):
        newimage[i,:]=img[pos[0]+i,pos[1]:(pos[1]+shape[1])]
    return newimage

def calBlockSize(lines):
    d=dict()
    for rho,theta in lines[:,0,:]:
        if theta not in d:
            d[theta]=[rho]
        else:
            d[theta].append(rho)
    ##should only have two values of theta
    boundary=np.zeros((4,2))
    blocksize=np.zeros((2,1),dtype=int)
    i=0
    j=0
    for key in d:
        d[key].sort()
        boundary[i,0]=key
        boundary[i,1]=d[key][0]
        boundary[i+1,0]=key
        boundary[i+1,1]=d[key][-1]
        dv=np.diff(d[key])
        blocksize[j]=int(np.average(dv))
        i=i+2
        j=j+1
    return boundary, blocksize
    
folder="C:\\Users\\eziyguo\\Desktop\\opencv\\"
filename="test1.jpg"

fig1=cv2.imread(folder+filename)

g_fig1=cv2.cvtColor(fig1,cv2.COLOR_BGR2GRAY)

ret,bina=cv2.threshold(imgChop(g_fig1,(15,15),(100,100)),127,255,cv2.THRESH_BINARY)
binv=(bina+1)*255 ##unit8


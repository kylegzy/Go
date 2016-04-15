import cv2
import numpy as np
from matplotlib import pyplot as plt

class ParallelLine(Exception):
    pass

def imgChop(img, pos, shape):
    newimage=np.zeros(shape,img.dtype)
    for i in range(shape[0]):
        newimage[i,:]=img[pos[0]+i,pos[1]:(pos[1]+shape[1])]
    return newimage

def getCandidateLines(d):
    max1=0
    max2=0
    maxkey1=''
    maxkey2=''
    for key in d:
        length=len(d[key])
        if  length > max1:
            max2=max1
            maxkey2=maxkey1
            max1=length
            maxkey1=key
        elif length > max2:
            max2=length
            maxkey2=key
        else:
            pass
        
    return [maxkey1, maxkey2]
    
def getEffectiveLines(l):
    
   
def calBlockSize(lines):
    d=dict()
    for rho,theta in lines[:,0,:]:
        if theta not in d:
            d[theta]=[rho]
        else:
            d[theta].append(rho)
            
    
    ##should only have two values of theta for an ideal image
    ##take the two thetas with the longest member
    
    candidate=getCandidateLines(d)
        
        
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
    
def calIntersection(rho1, theta1, rho2, theta2):
    if theta1 == theta2:
        raise ParallelLine
    else:
        x = (rho2*np.sin(theta1)-rho1*np.sin(theta2))/np.sin(theta1-theta2)
        y = (rho1*np.cos(theta2)-rho2*np.cos(theta1))/np.sin(theta1-theta2)
    
    return x,y
    
    
    
folder="C:\\Users\\eziyguo\\Desktop\\opencv\\"
filename="test1.jpg"

fig1=cv2.imread(folder+filename)

g_fig1=cv2.cvtColor(fig1,cv2.COLOR_BGR2GRAY)

ret,bina=cv2.threshold(imgChop(g_fig1,(15,15),(100,100)),127,255,cv2.THRESH_BINARY)
binv=(bina+1)*255 ##unit8


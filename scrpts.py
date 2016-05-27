# -*- coding: utf-8 -*-

import numpy as np
import cv2
import os, struct
from matplotlib import pyplot as plt
from array import array as pyarray

def imgChop(img, pos, shape,pos_type):
    if len(img.shape)==3:
        newimage=np.zeros(np.append(shape,img.shape[2]),img.dtype)
    else:
        newimage=np.zeros(shape,img.dtype)
    if pos_type == "upright":
        start1=pos[0]
        start2=pos[1]
    if pos_type == "center":
        start1=pos[0]-int(np.floor(shape[0]*0.5))
        start2=pos[1]-int(np.floor(shape[1]*0.5))
    for i in range(shape[0]):
        newimage[i,:]=img[start1+i,start2:(start2+shape[1])]
    return newimage
    

def warp(thetas, discount=2*np.pi):
    if isinstance(thetas, list):
        length=thetas.__len__()
        r=np.zeros(length,np.float64)
    elif isinstance(thetas,np.ndarray):
        length=thetas.size
        r=np.zeros(length,np.float64)
    elif isinstance(thetas,float) or isinstance(thetas,int) or isinstance(thetas,np.float64) or isinstance(thetas,np.float32):
        n=0
        if thetas > 0:
            while (thetas-discount*n)>0:
                n=n+1
            r=thetas-discount*(n-1)
            return r
        elif thetas < 0:
            while (thetas+discount*n)<0:
                n=n+1
            r=thetas+discount*(n)
            return r
        else:
            return thetas
        
    for i in range(length):
        n=0
        if thetas[i] > 0:
            while (thetas[i]-discount*n)>0:
                n=n+1
            r[i]=thetas[i]-discount*(n-1)
        elif thetas[i] < 0:
            while (thetas[i]+discount*n)<0:
                n=n+1
            r[i]=thetas[i]+discount*n
        
        
            
    return r
    
    
def quadrantCount(thetas,start):
    t=warp(thetas)
    c=[0,0,0,0]
    for i in range(t.size):
        if t[i] > start and t[i] < (start+np.pi*0.5):
            c[0]=c[0]+1
        elif t[i] > (start+np.pi*0.5) and t[i] < (start+np.pi):
            c[1]=c[1]+1
        elif t[i] > (start+np.pi) and t[i] < (start+np.pi*1.5):
            c[2]=c[2]+1
        else:
            c[3]=c[3]+1
    
    return c
            


def halfChordLength(rho, theta):
    return rho*np.sin(theta*0.5)
    
def DominantThetaEstimation(rhos, thetas, accuracy=np.pi/180, limit=45):
    n=int(np.floor(np.pi/2/accuracy))
    #print n
    d=np.full(n,np.inf,np.float64)
    for i in range(limit)+range(n-limit,n):
        phi=thetas-i*accuracy
        t=np.array([np.sin(phi), np.cos(phi)])
        m=np.min(np.absolute(t),axis=0)
        m=np.arcsin(m)
        #penalty=np.cos(np.min((thetas,np.absolute(45-thetas))))
        #m=m/penalty
        d[i]=np.sum(m)
            
    
    #print d    
    return np.deg2rad((d.argmin(),d.argmin()+90))

def majorityVote(data,lower_bound=0):
    d={}
    m=0
    r=0
    #print lower_bound
    #print data
    for i in data:
        #print i
        i=int(i)
        if i not in d:
            d[i]=0
        d[i]=d[i]+1
        if d[i]>m and i>lower_bound:
            m=d[i]
            r=i
            
    return r
    
        

def DominantRhoEstimation(rhos, thetas, dominant_thetas, min_line_gap, tolerance=np.deg2rad(2)):
    dominant_rhos = []   
    block=[]
    for t in dominant_thetas:    
        ix=np.where(np.absolute(np.sin(thetas-t))<np.sin(tolerance))
        s=np.sort(rhos[ix[0]],axis=0)
        #print s
        d=np.diff(s,axis=0)
        #print d
        common_diff=majorityVote(d,min_line_gap)
        #print common_diff
        """
        while common_diff < min_line_gap:
            ix1=np.where(d==common_diff)
            for i in range(len(ix1[0])):
                s[ix1[0][i]]=(s[ix1[0][i]]+s[ix1[0][i+1]])*0.5
            s=np.delete(s,ix1[0]+1)
            d=np.diff(s,axis=0)
            common_diff=np.median(d)
            print s
            print common_diff
            if s.shape[0]<2:
                break
        """    
        block.append(common_diff)
        ix1=np.where(d==common_diff)
        recont=[]
        for i in range(len(ix1[0])):
            if len(recont) == 0:
                recont.append(s[ix1[0][i]])
                recont.append(s[ix1[0][i]+1])
            elif recont[-1] == s[ix1[0][i]]:
                recont.append(s[ix1[0][i]+1])
            else:
                recont.append(s[ix1[0][i]])
                recont.append(s[ix1[0][i]+1])
        recont = np.array(recont)
        #print recont
        fisrtTerm = recont[0]
        termN=np.floor((recont-fisrtTerm)/common_diff)
        #print termN
        termN_size = termN.size
        #print termN_size
        if (termN_size == np.absolute(termN[-1] - termN[0] + 1)) and (termN_size >= 19):
            dominant_rhos.append(recont)
        else:
            if termN[-1] <= 19:
                i=0
                j=0
                while (i < termN[-1]):
                    #print "i = %d" % i
                    #print "j = %d" % j
                    #print "termN[i] = %d" % termN[i]
                    while (j < termN[i]):
                        #print "j = %d" % j
                        termN=np.insert(termN,j,j)
                        #print "%d" % (j*common_diff+fisrtTerm)
                        recont=np.insert(recont,j,j*common_diff+fisrtTerm)
                        j=j+1
                        i=i+1
                        #print termN
                        #print recont
                    j=j+1
                    i=i+1
            termN_size = termN.size
            dominant_rhos.append(recont)
            
    
    return dominant_rhos,block

def sequenceSum(a1, diff, n):
    an = a1 + (n-1)*diff
    if np.mod(n,2) == 0:
        return (a1+an)*n/2
    else:
        return ((a1+an)*(n-1)/2 + a1+(n-1)/2*diff)

    
def hough(fig,rho_accruracy=1, theta_accuracy=np.pi/180, rho_lower_bound=0.3):

    lines = cv2.HoughLines(fig,rho_accruracy,theta_accuracy,int(np.round(fig.shape[0]*rho_lower_bound)))        
    lines=np.float64(lines)
    ix=np.where(lines[:,:,0]<0)
    lines[:,:,1][ix[0]]=lines[:,:,1][ix[0]] + np.pi
    lines[:,:,0][ix[0]]=lines[:,:,0][ix[0]] * (-1)
    return lines
    
def myhough(fig,rho_accruracy=1, theta_accuracy=np.pi/180, rho_num=100):
    c=int(np.floor(np.sqrt(fig.shape[0]**2+fig.shape[1]**2)))
    r=int(np.floor(np.pi/theta_accuracy))
    
    
    
    h=np.zeros((c,r))
    
    result=np.zeros(rho_num,np.dtype((float,(1,2))))
    
    for i in range(fig.shape[0]):
        #print i
        for j in range(fig.shape[1]):
            #print j
            if fig[i][j] == 0 :
                for k in range(r):
                    l=int(np.floor(i*np.cos(np.deg2rad(k))+j*np.sin(np.deg2rad(k))))
                    h[l][k]=h[l][k]+1
    flatten=np.ravel(h)
    ix=np.argsort(flatten)
    l=ix.size
    
    for i in range(rho_num):
        
        coordinate=ix[l-i-1]
        k=np.mod(coordinate,c)
        if k==0:
            p=coordinate
        else:
            p=np.floor_divide(coordinate/k)
        result[i]=(p,np.deg2rad(k))
    
    return result
    
def calIntersection(rhos1, thetas1, rhos2, thetas2, img_shape):
    s = np.zeros((rhos1.size,rhos2.size),dtype=np.dtype(object))
    for i in range(rhos1.size):
        for j in range(rhos2.size):
            rho1=rhos1[i]
            theta1=thetas1[i]
            rho2=rhos2[j]
            theta2=thetas2[j]
            x = int((rho2*np.sin(theta1)-rho1*np.sin(theta2))/np.sin(theta1-theta2))
            y = int((rho1*np.cos(theta2)-rho2*np.cos(theta1))/np.sin(theta1-theta2))
            #print x,y            
            if x > 0 and x < img_shape[0] and y > 0 and y < img_shape[1]:
                s[i,j] = (x,y)
                
    return s
    
def plotPredictPoint(img,point):
    n=img.copy()
    for i in range(point.shape[0]):
        for j in range(point.shape[1]):
            n=cv2.circle(n,point[i][j],5,(0,0,255),-1)
    return n
    
def plotbox(img,cnt):
    n=img.copy()
    x,y,w,h = cv2.boundingRect(cnt)
    n = cv2.rectangle(n,(x,y),(x+w,y+h),(255,255,255),1)
    plt.imshow(n,'gray')

def sigmoid(x,deriv=False):
    if(deriv==True):
        s = sigmoid(x,False)
        return s*(1-s)
    return 1/(1+np.exp(-x))
    
def adaptiveThresh(img):
    m=int(np.mean(img))
    if m<=127:
        #black
        #ret,fa=cv2.threshold(img,m,255,cv2.THRESH_BINARY)
        ret,fa=cv2.threshold(img,180,255,cv2.THRESH_BINARY)
    if m>127:
        #white
        #ret,fa=cv2.threshold(img,m,255,cv2.THRESH_BINARY)
        ret,fa=cv2.threshold(img,80,255,cv2.THRESH_BINARY)
        fa=(fa+1)*255
    return fa 
    
def roiNormalize(img,targetSize=(28,28),margin=0.1):
    #suppose height (row) > width (column) for a digit    
    r=np.zeros(targetSize,np.uint8)
    m=margin*targetSize[0]
    alpha=(targetSize[0]-2.0*m)/img.shape[0]
    new_img=cv2.resize(img,None,fx=alpha,fy=alpha, interpolation=cv2.INTER_LINEAR)
    x=int(np.floor((targetSize[1]-new_img.shape[1])/2))
    y=int(np.floor((targetSize[0]-new_img.shape[0])/2))
    r[y:(y+new_img.shape[0]),x:(x+new_img.shape[1])]=new_img
    return r
    
def load_mnist(dataset="training", digits=None, path=None, asbytes=False, selection=None, return_labels=True, return_indices=False):
    """
    Loads MNIST files into a 3D numpy array.

    You have to download the data separately from [MNIST]_. It is recommended
    to set the environment variable ``MNIST`` to point to the folder where you
    put the data, so that you don't have to select path. On a Linux+bash setup,
    this is done by adding the following to your ``.bashrc``::

        export MNIST=/path/to/mnist

    Parameters
    ----------
    dataset : str 
        Either "training" or "testing", depending on which dataset you want to
        load. 
    digits : list 
        Integer list of digits to load. The entire database is loaded if set to
        ``None``. Default is ``None``.
    path : str 
        Path to your MNIST datafiles. The default is ``None``, which will try
        to take the path from your environment variable ``MNIST``. The data can
        be downloaded from http://yann.lecun.com/exdb/mnist/.
    asbytes : bool
        If True, returns data as ``numpy.uint8`` in [0, 255] as opposed to
        ``numpy.float64`` in [0.0, 1.0].
    selection : slice
        Using a `slice` object, specify what subset of the dataset to load. An
        example is ``slice(0, 20, 2)``, which would load every other digit
        until--but not including--the twentieth.
    return_labels : bool
        Specify whether or not labels should be returned. This is also a speed
        performance if digits are not specified, since then the labels file
        does not need to be read at all.
    return_indicies : bool
        Specify whether or not to return the MNIST indices that were fetched.
        This is valuable only if digits is specified, because in that case it
        can be valuable to know how far
        in the database it reached.

    Returns
    -------
    images : ndarray
        Image data of shape ``(N, rows, cols)``, where ``N`` is the number of images. If neither labels nor inices are returned, then this is returned directly, and not inside a 1-sized tuple.
    labels : ndarray
        Array of size ``N`` describing the labels. Returned only if ``return_labels`` is `True`, which is default.
    indices : ndarray
        The indices in the database that were returned.

    Examples
    --------
    Assuming that you have downloaded the MNIST database and set the
    environment variable ``$MNIST`` point to the folder, this will load all
    images and labels from the training set:

    >>> images, labels = ag.io.load_mnist('training') # doctest: +SKIP

    Load 100 sevens from the testing set:    

    >>> sevens = ag.io.load_mnist('testing', digits=[7], selection=slice(0, 100), return_labels=False) # doctest: +SKIP

    """

    # The files are assumed to have these names and should be found in 'path'
    files = {
        'training': ('train-images.idx3-ubyte', 'train-labels.idx1-ubyte'),
        'testing': ('t10k-images.idx3-ubyte', 't10k-labels.idx1-ubyte'),
    }

    if path is None:
        try:
            path = os.environ['MNIST']
        except KeyError:
            raise ValueError("Unspecified path requires environment variable $MNIST to be set")

    try:
        images_fname = os.path.join(path, files[dataset][0])
        labels_fname = os.path.join(path, files[dataset][1])
    except KeyError:
        raise ValueError("Data set must be 'testing' or 'training'")

    # We can skip the labels file only if digits aren't specified and labels aren't asked for
    if return_labels or digits is not None:
        flbl = open(labels_fname, 'rb')
        magic_nr, size = struct.unpack(">II", flbl.read(8))
        labels_raw = pyarray("b", flbl.read())
        flbl.close()

    fimg = open(images_fname, 'rb')
    magic_nr, size, rows, cols = struct.unpack(">IIII", fimg.read(16))
    images_raw = pyarray("B", fimg.read())
    fimg.close()

    if digits:
        indices = [k for k in range(size) if labels_raw[k] in digits]
    else:
        indices = range(size)

    if selection:
        indices = indices[selection] 
    N = len(indices)

    images = np.zeros((N, rows, cols), dtype=np.uint8)

    if return_labels:
        labels = np.zeros((N), dtype=np.int8)
    for i, index in enumerate(indices):
        images[i] = np.array(images_raw[ indices[i]*rows*cols : (indices[i]+1)*rows*cols ]).reshape((rows, cols))
        if return_labels:
            labels[i] = labels_raw[indices[i]]

    if not asbytes:
        images = images.astype(float)/255.0

    ret = (images,)
    if return_labels:
        ret += (labels,)
    if return_indices:
        ret += (indices,)
    if len(ret) == 1:
        return ret[0] # Don't return a tuple of one
    else:
        return ret
    
    

folder="C:\\Users\\eziyguo\\Documents\\GitHub\\Go\\"
filename="test1.jpg"

fig1=cv2.imread(folder+filename)

g_fig1=cv2.cvtColor(fig1,cv2.COLOR_BGR2GRAY)

#ret,bina=cv2.threshold(imgChop(g_fig1,(15,15),(100,100)),127,255,cv2.THRESH_BINARY)
#binv=(bina+1)*255 ##unit8

ret,fa=cv2.threshold(g_fig1,127,255,cv2.THRESH_BINARY)
fv=(fa+1)*255

#edges = cv2.Canny(fv,fv.shape[0]*0.5,fv.shape[0])
edges = cv2.Laplacian(fa,cv2.CV_8U)
f=hough(edges)

dthetas=DominantThetaEstimation(f[:,:,0],f[:,:,1])

drhos,block=DominantRhoEstimation(f[:,:,0],f[:,:,1],dthetas,fa.shape[0]/38)

p=calIntersection(drhos[0],[dthetas[0]]*19,drhos[1],[dthetas[1]]*19,fa.shape)
import numpy as np
import cv2
from scipy import signal
import math
import matplotlib.pyplot as plt
from PIL import Image

def rescale(matrix):
    scaled = ((matrix - matrix.min()) * (255 / (matrix.max() - matrix.min()))).astype(np.uint8)
    return scaled

def computeC2(sigma, theta, filterRadius):
    numerater = 0
    denominator = 0
    e = [math.cos(theta), math.sin(theta)]
    num = (np.pi/(2*sigma))
    for x in range(-filterRadius, filterRadius+1):
        for y in range(-filterRadius, filterRadius+1):
            u = [x,y]
            u2 = (x**2) + (y**2)
            denominator += math.exp(-u2/(2*(sigma**2)))
            numerater += complex(math.cos(num*np.dot(u, e)),math.sin(num*np.dot(u,e)))*math.exp(-u2/(2*(sigma**2)))
    c2 = numerater / denominator
    return c2

def computeC1(sigma, theta, c2, filterRadius):
    z=0
    e = [math.cos(theta), math.sin(theta)]
    num = np.pi / (2*sigma)
    for x in range(-filterRadius, filterRadius+1):
        for y in range(-filterRadius, filterRadius+1):
            u = [x,y]
            u2 = (x**2)+(y**2)
            z += ((1 - 2 * c2 * math.cos(num*np.dot(u,e))+(c2**2))*math.exp(-u2/(sigma**2)))
    c1 = sigma/(z**0.5)
    return c1

def psiFunc(x,y,c1,c2,sigma,theta):
    num = np.pi/(2*sigma)
    u2 = (x**2)+(y**2)
    u = [x,y]
    e = [math.cos(theta), math.sin(theta)]
    psi = (c1/sigma) * (complex(math.cos(num*np.dot(u,e)), math.sin(num*np.dot(u,e))) - c2)* math.exp(-u2/(2*(sigma**2)))
    return psi

def makeWavelet(sigma, theta, filterRadius):
    c2 = computeC2(sigma, theta, filterRadius)
    c1 = computeC1(sigma, theta, c2, filterRadius)
    for x in range(-filterRadius, filterRadius+ 1):
        for y in range(-filterRadius, filterRadius+ 1):
            morlet = psiFunc(x*1., y*1., c1, c2, sigma, theta)
            morletReal[x+filterRadius][y+filterRadius] = morlet.real
            morletImaginary[x+filterRadius][y+filterRadius] = morlet.imag
    return

def makeWaveletList(Sigma, Theta, imgRad):
    if isinstance(Sigma, int) and (isinstance(Theta, float) or isinstance(Theta,int)):
        makeWavelet(Sigma, Theta, imgRad)
        rList.append(np.matrix.copy(morletReal))
        iList.append(np.matrix.copy(morletImaginary))
    else:
        for sigma in Sigma:
            for theta in Theta:
                makeWavelet(sigma, theta, imgRad)
                rList.append(np.matrix.copy(morletReal))
                iList.append(np.matrix.copy(morletImaginary))
    return

def createGaussian(matrix, sigma):
    [r, c] = matrix.shape
    for i in range(0, r):
        x = (i - imgRad) * 1.
        for j in range(0, c):
            y = (j - imgRad) * 1.
            matrix[i][j] = (1/(2*np.pi*sigma**2)) * math.exp(-(x**2 + y**2)/(2*sigma**2))
    return

def plotWavelets():
    plt.suptitle("Real Wavelets")
    for i in range (0,len(rList)):
        realScaled = rescale(rList[i]) #(255.0 / (rList[i].max() - rList[i].min()) * (rList[i] - rList[i].min())).astype(np.uint8)
        realImage = Image.fromarray(realScaled)
        plt.subplot(3,4, i+1)
        plt.imshow(realImage, cmap='gray')
    plt.show()

    plt.suptitle('Imaginary Wavelets')
    for i in range (0,len(iList)):
        imaginaryScaled = rescale(iList[i]) #(255.0 / (iList[i].max() - iList[i].min()) * (iList[i] - iList[i].min())).astype(np.uint8)
        imaginaryImage = Image.fromarray(imaginaryScaled)
        plt.subplot(3,4, i+1)
        plt.imshow(imaginaryImage, cmap='gray')
    plt.show()
    return

def convolve(picture):
    for i in range (0,len(rList)):
        convolvedR = signal.convolve2d(picture, rList[i])
        rList[i] = np.matrix.copy(convolvedR)
        convolvedI = signal.convolve2d(picture, iList[i])
        iList[i] = np.matrix.copy(convolvedI)
    return

def plotConvolution(picture):
    plt.suptitle("Real Convolutions")
    for i in range(0, len(rList)):
        plt.subplot(3,4, i + 1)
        plt.imshow(rList[i], cmap='gray')
    plt.show()

    plt.suptitle("Imaginary Convolutions")
    for i in range(0, len(rList)):
        plt.subplot(3,4, i + 1)
        plt.imshow(iList[i], cmap='gray')
    plt.show()

    plt.suptitle("Gaussian Convolution")
    convolved = signal.convolve2d(picture, gaussian)
    plt.imshow(convolved, cmap='gray')
    plt.show()
    return

def findMaxMatrix(matrix, mList):
    [r,c] = matrix.shape
    for i in range (0, r):
        for j in range(0, c):
            for n in range (0, len(mList)):
                if matrix[i][j] < abs(mList[n][i][j]):
                    matrix[i][j] = abs(mList[n][i][j])
    return

def plotMatrixHistogram(matrixR, matrixI):
    rescaledR = rescale(matrixR) #((matrixR - matrixR.min()) * (255 / (matrixR.max() - matrixR.min()))).astype(np.uint8)
    rescaledI = rescale(matrixI) #((matrixI - matrixI.min()) * (255 / (matrixI.max() - matrixI.min()))).astype(np.uint8)

    rImg = Image.fromarray(rescaledR)
    iImg = Image.fromarray(rescaledI)
    plt.subplot(1,2,1)
    plt.title("W_max_real")
    plt.imshow(rImg, cmap='gray')
    plt.subplot(1,2,2)
    plt.title("W_max_imaginary")
    plt.imshow(iImg, cmap='gray')
    plt.show()

    plt.subplot(2,1,1)
    plt.title("W_min_real histogram")
    plt.hist(rescaledR)
    plt.subplot(2,1,2)
    plt.hist(rescaledI)
    plt.title("W_max_Imaginary Histogram")
    plt.show()
    return

def ratioEdgeDetection(Wreal, Wimaginary, epsilon):
    ratio = (Wreal + epsilon * 0.001) / (Wimaginary + epsilon)
    D = ratio.max()
    edge = 1 - (ratio/D)
    return edge

def diffEdgeDetection(Wreal, Wimaginary, alpha):
    diff = Wimaginary - Wreal
    edge = np.exp(-(alpha * diff))
    return edge

def plotEdge(ratio, diff):
    rescaledRatio = rescale(ratio)
    rescaledDiff = rescale(diff)
    plt.subplot(1,2,1)
    plt.title("Ratio Edge Detectiion")
    plt.imshow(rescaledRatio, cmap = 'gray')
    plt.subplot(1, 2, 2)
    plt.title("Difference Edge Detection")
    plt.imshow(rescaledDiff, cmap = 'gray')
    plt.show()
    return

#    inputs
Sigma = [1, 3, 6]
Theta = [0, np.pi/4, np.pi/2, np.pi*3/4]
imgRad = 18

#    do not change
morletReal = np.zeros((imgRad*2 + 1, imgRad*2 + 1))
morletImaginary = np.zeros((imgRad*2 + 1, imgRad*2 + 1))
rList = []
iList = []

#    creating wavelets and plotting
makeWaveletList(Sigma, Theta, imgRad)
plotWavelets() ###

#    convoluting and plotting
circle = cv2.imread('circle.jpg',0)
gaussian = np.zeros((imgRad*2 +1, imgRad*2+1))
createGaussian(gaussian, 6)
convolve(circle)
plotConvolution(circle) ###

#    Wreal, Wimaginary
Wreal= np.zeros(rList[0].shape)
Wimaginary = np.zeros(iList[0].shape)
findMaxMatrix(Wreal, rList)
findMaxMatrix(Wimaginary, iList)

#    plot W_max
plotMatrixHistogram(Wreal, Wimaginary) ###

#    edge detection and plotting
epsilon = 1000000
alpha = 0.0000001
edgeRatio = ratioEdgeDetection(Wreal, Wimaginary, epsilon)
edgeDiff =  diffEdgeDetection(Wreal, Wimaginary, alpha)
plotEdge(edgeRatio, edgeDiff)
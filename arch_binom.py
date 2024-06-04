from numpy import *
import numpy as np
from scipy import signal
set_printoptions(precision = 4, suppress = True, linewidth = 300, threshold = 100000)

import ctypes
import sys
import time
import fcntl

# Load the shared library
lib = ctypes.CDLL('./arch_binomial.so')

# Define the function argument and return types

lib.Archetypal.argtypes = [ctypes.POINTER(ctypes.c_int),
                           ctypes.POINTER(ctypes.c_int),
                           ctypes.POINTER(ctypes.c_int),
                           ctypes.POINTER(ctypes.c_double),
                           ctypes.POINTER(ctypes.c_double),
                           ctypes.c_int,
                           ctypes.c_int,
                           ctypes.c_int,
                           ctypes.POINTER(ctypes.c_double)]
lib.Archetypal.restype = ctypes.c_double

lib.binomial_like.argtypes = [ctypes.c_double,
                              ctypes.c_double,
                              ctypes.c_double]
lib.binomial_like.restype = ctypes.c_double



def Arch(datos, noc = 3, analisis = "ut", aic0 = 1e99):
    aic1 = 1e99
    datos = float64(datos)
    maxbin = datos.max(0)
    archivo = "arquetipos_" + analisis + "_%02d.npz"%(noc)
    
    XC = zeros((noc, datos.shape[1])).astype(float64)
    C  = zeros((noc, datos.shape[0])).astype(float64)
    S  = zeros((datos.shape[0], noc)).astype(float64)
    IE = zeros((3,)).astype(float64)

    XC_ptr = XC.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
    A_ptr  = datos.ctypes.data_as (ctypes.POINTER(ctypes.c_int))
    N_ptr  = maxbin.ctypes.data_as (ctypes.POINTER(ctypes.c_int))
    C_ptr  = C.ctypes.data_as (ctypes.POINTER(ctypes.c_double))
    S_ptr  = S.ctypes.data_as (ctypes.POINTER(ctypes.c_double))
    IE_ptr = IE.ctypes.data_as (ctypes.POINTER(ctypes.c_double))
    SSe = lib.Archetypal(A_ptr,       N_ptr,       XC_ptr,         C_ptr,  S_ptr,
                         XC.shape[0], XC.shape[1], datos.shape[0], IE_ptr)

    lk, bic, aic = IE

    if aic < aic0:
        try:
            fl = open(archivo, "rb")
            fcntl.flock(fl, fcntl.LOCK_EX)
            aic1 = load(fl)["ajuste"][2]
            fcntl.flock(fl, fcntl.LOCK_UN)
            fl.close()
        except IOError: aic1 = 1e99
        #print("Obtenido AIC1", IE, aic0, aic1)
        
        open("arquetipos_" + analisis + "_%02d.log"%(noc),"a").write("aic %12.4f aic1 %12.4f\n"%(aic, aic1))
                                                                                     
        if aic < aic1:
            print("Guardando", IE, aic0, aic1)
            fl = open(archivo, "wb")
            fcntl.flock(fl, fcntl.LOCK_EX)
            savez_compressed(fl, XC = XC, S = S, C = C, ajuste = array([SSe, lk, aic, bic]))
            fcntl.flock(fl, fcntl.LOCK_UN)
            fl.close()
    aic = min(aic0, min(aic, aic1))
    #print("Retornando nuevo AIC", aic)

    return lk, aic, bic



def materiales(lista):
    """
    1= Adobe 1
    2= Ladrillo 4 
    3= Chapa  3
    4= Madera  3
    5= Ladrillo y Chapa 4
    6= Ladrillo, Chapa y Madera 4
    7= Chapa y Madera 3
    8= Adobe y Madera 2
    9= Otro -1
    """
    materiales = zeros_like(lista)
    for i in range(len(lista)):
        if lista[i] == 1: materiales[i] = 1
        elif lista[i] == 2: materiales[i] = 4
        elif lista[i] == 3: materiales[i] = 3
        elif lista[i] == 4: materiales[i] = 3
        elif lista[i] == 5: materiales[i] = 4
        elif lista[i] == 6: materiales[i] = 4
        elif lista[i] == 7: materiales[i] = 3
        elif lista[i] == 8: materiales[i] = 2
        else: materiales[i] = -1
    return materiales

def loadmaras():
    #texto = open("lepra_socio.txt",'r').readlines()
    texto = open("/home/caronte/papers_en_redaccion_2/lepra_formosa/lepra_socio_3.txt",'r').readlines()
    matriz = []
    for linea in texto[2:]:
        linea = linea.split(";")
        matriz.append([])
        for celda in linea[1:]:
            try:               celda = int(celda)
            except ValueError: celda = -1
            matriz[-1].append(celda)
    matriz = array(matriz)
    matriz[:,4] = materiales(matriz[:,4])
    matriz[1:,:] -= 1
    return array(matriz).transpose()

def loadARM():
    texto = open("../lepra_formosa/lepra_socio_3.txt",'r').readlines()
    matriz = zeros((len(texto), 17))
    j = 0
    for linea in texto[2:]:
        linea = linea.split(";")
        linea[-4:-1] = linea[-3:]; linea = linea[:-1]
        i = 0
        for celda in linea[1:-1]:
            try:               celda = int(celda)
            except ValueError: celda = -1
            matriz[j, i] = celda
            i += 1
        ocup = int(linea[-1])
        if   ocup == 4: matriz[j, -4]        = 1
        elif ocup >= 1: matriz[j, -4 + ocup] = 1
        else:           matriz[j, -4]        = 1
        j += 1
    matriz[:,1:-6] -= 1
    loc = loadtxt("../lepra_formosa/urbano_rural.txt")
    matriz[:-2,-5]   = loc - 1
    return array(matriz).transpose()[:,:-2].copy()

def loadUT():
    texto = open("../lepra_formosa/lepra_socio_2.txt",'r').readlines()
    ocupv = array([11, 11, 11, 8, 1, 3, 12, 11, 2, 2, 3, 1, 2, 7, 1, 7, 1,
                  2, 8, 7, 1, 5, 10, 7, 11, 8, 8, 5, 8, 8, 8, 4, 1, 8, 7,
                  8, 10, 1, 2, 7, 10, 8, 11, 11, 8, 2, 2, 8, 7, 7, 4, 3,
                  11, 1, 7, 9, 4, 2, 3, 1, 1, 1, 8, 10, 11, 7, 8, 10, 7,
                  11, 6, 1, 3, 2, -1, 1, 8, 1, 8, 1, 1, 3, 2])

    matriz = zeros((len(texto), 17))
    j = 0
    for linea, ocup in zip(texto, ocupv):
        linea = linea.split("\t")
        i = 0
        for celda in linea[1:-1]:
            try:               celda = int(celda)
            except ValueError: celda = -1
            matriz[j, i] = celda
            i += 1
        #ocup = int(linea[-1])
        if ocup >= 1:
            matriz[j, ocup + 4] = 1
        else:
            matriz[j, 5:] = -1
        j += 1
    #matriz[:,-1]   = loc - 1
    matriz[:,1:5] -= 1
    return array(matriz).transpose()

def main(n = -1, distribution = 0, nocmin = 1, nocmax = 6, analisis = "ut"):
    if analisis == "ut":
        M = loadUT()
    else:
        M = loadARM()
    M = int32(M)
    aicm = zeros((10))
    bicm = zeros((10))
    for i in range(nocmin, nocmax):
        aic0 = 1e99
        lk, aic, bic = Arch(M.T.copy(), noc = i + 1, aic0 = aic0, analisis = analisis)
        aic0 = min(aic, aic0)
        aicm[i] = aic
        bicm[i] = bic

if __name__ == '__main__':
    main()
    main(analisis = "arm")

# -*- coding: utf-8 -*-
"""
IMPLEMENTACIÓN DE MATRICES CON COEFICIENTES EN Z/pZ

Consideraciones:
    1) A = [A0, A1, ..., An] donde Ai son vectores columna
"""

from copy import deepcopy
import cuerpoZpZ as ZpZ

def neutro_ad(m, n, p):
    """

    Parámetros
    ----------
    m : ENTERO
        NÚMERO DE FILAS DE LA MATRIZ
    n : ENTERO
        NÚMERO DE COLUMNAS DE LA MATRIZ
    p : ENTERO
        PRIMO

    Devuelve
    --------
    LISTA
        MATRIZ NEUTRO ADITIVO DE TAMAÑO m x n EN Z/pZ

    """
    return [[ZpZ.neutro_ad(p) for i in range(m)] for j in range(n)]

def neutro_mult(n, p):
    """

    Parámetros
    ----------
    n : ENTERO
        NÚMERO DE FILAS Y COLUMNAS DE LA MATRIZ
    p : ENTERO
        PRIMO

    Devuelve
    --------
    LISTA
        MATRIZ NEUTRO MULTIPLICATIVO DE TAMAÑO n x n EN Z/pZ

    """
    return [[ZpZ.neutro_ad(p)]*i + [ZpZ.neutro_mult(p)] + \
    [ZpZ.neutro_ad(p)]*(n-i-1) for i in range(n)]


def suma(A, B, p):
    """

    Parámetros
    ----------
    A : LISTA
        MATRIZ
    B : LISTA
        MATRIZ
    p : ENTERO
        PRIMO

    Devuelve
    --------
    s : LISTA
        SUMA DE LAS MATRICES A Y B EN Z/pZ

    """
    m = len(A[0])
    n = len(A)
    s = neutro_ad(m, n, p)
    for i in range(0, m):
        for j in range(0, n):
            s[i][j] = ZpZ.suma(A[i][j], B[i][j], p)
    return s    

def inv_ad(A, p):
    """

    Parámetros
    ----------
    A : LISTA
        MATRIZ
    p : ENTERO
        PRIMO

    Devuelve
    --------
    inv : LISTA
          MATRIZ INVERSA ADITIVA DE LA MATRIZ A EN Z/pZ

    """
    m = len(A[0])
    n = len(A)
    inv = neutro_ad(m, n, p)
    for i in range(0, m):
        for j in range(0, n):
            inv[i][j] = ZpZ.inv_ad(A[i][j], p)
    return inv  
                      
def mult(A, B, p):
    """

    Parámetros
    ----------
    A : LISTA
        MATRIZ
    B : LISTA
        MATRIZ
    P : ENTERO
        PRIMO

    Devuelve
    --------
    prod : LISTA
           PRODUCTO DE LAS MATRICES A Y B EN Z/pZ

    """
    m1 = len(A[0])
    n1 = len(A)
    m2 = len(B[0])
    n2 = len(B)
    prod = neutro_ad(m1, n2, p)
    for i in range(0, n2):
        for j in range(0, m1):
            prod[i][j] = ZpZ.neutro_ad(p)
            for k in range(0, n1):
                m = ZpZ.mult(A[k][j], B[i][k], p)
                prod[i][j] = ZpZ.suma(prod[i][j], m, p)       
    return prod    

def det(A, p):
    """

    Parámetros
    ----------
    A : LISTA
        MATRIZ
    p : ENTERO
        PRIMO

    Devuelve
    --------
    ENTERO
        DETERMINANTE DE LA MATRIZ A EN Z/pZ

    """
    if len(A) == 1:
        return A[0][0]
    if len(A) == 2:
        return ZpZ.suma(ZpZ.mult(A[0][0], A[1][1], p), \
        ZpZ.inv_ad(ZpZ.mult(A[0][1], A[1][0], p), p), p)
    det_val = ZpZ.neutro_ad(p)
    for j in range(len(A)):
        B = deepcopy(A)
        B.pop(j)
        for k in range(len(B)):
            B[k] = B[k][1:]
        cofactor = det(B, p)
        if j % 2 == 1:
            cofactor = ZpZ.inv_ad(cofactor, p)  
        det_val = ZpZ.suma(det_val, ZpZ.mult(A[j][0], cofactor, p), p)
    return det_val

def eliminacion_gaussiana(A, p):
    """

    Parámetros
    ----------
    A : LISTA
        MATRIZ
    p : ENTERO
        PRIMO

    Devuelve
    --------
    A : LISTA
        MATRIZ RESULTANTE DE APLICAR LA ELIMINACIÓN GAUSSIANA 
        A LA MATRIZ A EN Z/pZ

    """
    
    m = len(A[0])
    n = len(A)
    pivotes = []
    for j in range(0, n):
        i = 0
        while i<m and (A[j][i] == ZpZ.neutro_ad(p) or (i in pivotes)):
            i += 1
        if i<m:
            pivotes.append(i)
            inv = ZpZ.inv_mult(A[j][i], p)
            for k in range(j, n):
                A[k][i] = ZpZ.mult(inv, A[k][i], p)
            for l in range(0, m):
                if l != i:
                    factor = A[j][l]
                    for s in range(j, n):
                        A[s][l] = ZpZ.suma(A[s][l], \
                        ZpZ.mult(ZpZ.inv_ad(factor, p), A[s][i], p), p)
    return A

def rand(m, n, p):
    """

    Parámetros
    ----------
    m : ENTERO
        NÚMERO DE FILAS DE LA MATRIZ
    n : ENTERO
        NÚMERO DE COLUMNAS DE LA MATRIZ
    p : ENTERO
        PRIMO

    Devuelve
    -------
    A : LISTA
        MATRIZ DE COEFICIENTES ALEATORIOS PERTENECIENTES A Z/pZ

    """
    A = neutro_ad(m, n, p)
    for i in range(m):
        for j in range(n):
            A[j][i] = ZpZ.rand(p)
    return A
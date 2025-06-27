# -*- coding: utf-8 -*-
"""
IMPLEMENTACIÓN DE MATRICES CON COEFICIENTES EN Fq

Consideraciones:
    1) A = [A0, A1, ..., An] donde Ai son vectores columna
"""

import random
from copy import deepcopy
import cuerpoFq as Fq


def neutro_ad(m, n, p, f):
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
        MATRIZ NEUTRO ADITIVO DE TAMAÑO m x n EN Fq

    """
    return [[Fq.neutro_ad(p, f) for i in range(m)] for j in range(n)]

def neutro_mult(n, p, f):
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
        MATRIZ NEUTRO MULTIPLICATIVO DE TAMAÑO n x n EN Fq

    """
    
    return [[Fq.neutro_ad(p, f)]*i + [Fq.neutro_mult(p, f)] + \
    [Fq.neutro_ad(p, f)]*(n-i-1) for i in range(n)]

def suma(A, B, p, f):
    """

    Parámetros
    ----------
    A : LISTA
        MATRIZ
    B : LISTA
        MATRIZ
    p : ENTERO
        PRIMO
    f : LISTA
        POLINOMIO

    Devuelve
    --------
    s : LISTA
        SUMA DE LAS MATRICES A Y B EN Fq

    """
    m = len(A[0])
    n = len(A)
    s = neutro_ad(m, n, p, f)
    for i in range(0, m):
        for j in range(0, n):
            s[i][j] = Fq.suma(A[i][j], B[i][j], p, f)
    return s     

def inv_ad(A, p, f):
    """

    Parámetros
    ----------
    A : LISTA
        MATRIZ
    p : ENTERO
        PRIMO
    f : LISTA
        POLINOMIO

    Devuelve
    --------
    inv : LISTA
          MATRIZ INVERSA ADITIVA DE LA MATRIZ A EN Fq

    """
    m = len(A[0])
    n = len(A)
    inv = neutro_ad(m, n, p, f)
    for i in range(0, m):
        for j in range(0, n):
            inv[i][j] = Fq.inv_ad(A[i][j], p, f)
    return inv 
                
def mult(A, B, p, f):
    """

    Parámetros
    ----------
    A : LISTA
        MATRIZ
    B : LISTA
        MATRIZ
    p : ENTERO
        PRIMO
    f : LISTA
        POLINOMIO


    Devuelve
    --------
    prod : LISTA
           PRODUCTO DE LAS MATRICES A Y B EN Fq

    """
    m1 = len(A[0])
    n1 = len(A)
    m2 = len(B[0])
    n2 = len(B)
    prod = neutro_ad(m1, n2, p, f)
    for i in range(0, n2):
        for j in range(0, m1):
            prod[i][j] = Fq.neutro_ad(p, f)
            for k in range(0, n1):
                m = Fq.mult(A[k][j], B[i][k], p, f)
                prod[i][j] = Fq.suma(prod[i][j], m, p, f)       
    return prod 
           
def det(A, p, f):
    """

    Parámetros
    ----------
    A : LISTA
        MATRIZ
    p : ENTERO
        PRIMO
    f : LISTA
        POLINOMIO

    Devuelve
    --------
    ENTERO
        DETERMINANTE DE LA MATRIZ A EN Fq

    """
    if len(A) == 1:
        return A[0][0]
    if len(A) == 2:
        return Fq.suma(Fq.mult(A[0][0], A[1][1], p, f), \
        Fq.inv_ad(Fq.mult(A[0][1], A[1][0], p, f), p, f), p, f)
    det_val = Fq.neutro_ad(p, f)
    for j in range(len(A)):
        B = deepcopy(A)
        B.pop(j)
        for k in range(len(B)):
            B[k] = B[k][1:]
        cofactor = det(B, p, f)
        if j % 2 == 1:
            cofactor = Fq.inv_ad(cofactor, p, f)  
        det_val = Fq.suma(det_val, Fq.mult(A[j][0], cofactor, p, f), p, f)
    return det_val

def eliminacion_gaussiana(A, p, f):
    """

    Parámetros
    ----------
    A : LISTA
        MATRIZ
    p : ENTERO
        PRIMO
    f : LISTA
        POLINOMIO

    Devuelve
    --------
    A : LISTA
        MATRIZ RESULTANTE DE APLICAR LA ELIMINACIÓN GAUSSIANA
        A LA MATRIZ A EN Fq

    """
    m = len(A[0])
    n = len(A)
    pivotes = []
    for j in range(0, n):
        i = 0
        while i<m and (A[j][i] == Fq.neutro_ad(p, f) or (i in pivotes)):
            i += 1
        if i<m:
            pivotes.append(i)
            inv = Fq.inv_mult(A[j][i], p, f)
            for k in range(j, n):
                A[k][i] = Fq.mult(inv, A[k][i], p, f)
            for l in range(0, m):
                if l != i:
                    factor = A[j][l]
                    for s in range(j, n):
                        A[s][l] = Fq.suma(A[s][l], \
                        Fq.mult(Fq.inv_ad(factor, p, f), A[s][i], p, f), \
                        p, f)
    return A

def rand(m, n, p, f):
    """

    Parámetros
    ----------
    m : ENTERO
        NÚMERO DE FILAS DE LA MATRIZ
    n : ENTERO
        NÚMERO DE COLUMNAS DE LA MATRIZ
    p : ENTERO
        PRIMO
    f : LISTA
        POLINOMIO

    Devuelve
    --------
    A : LISTA
        MATRIZ DE COEFICIENTES ALEATORIOS PERTENECIENTES A Fq

    """
    A = neutro_ad(m, n, p)
    for i in range(m):
        for j in range(n):
            deg = random.randint(0, len(f)-2)
            A[j][i] = Fq.rand(deg, p, f)
    return A
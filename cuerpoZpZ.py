# -*- coding: utf-8 -*-
"""
IMPLEMENTACIÓN DEL CUERPO Z/pZ
"""

import random
import auxiliar as aux

def neutro_ad(p):
    """
    
    Parámetros
    ----------
    p : ENTERO
        PRIMO

    Devuelve
    --------
    ENTERO
        NEUTRO ADITIVO DE Z/pZ

    """
    return 0

def neutro_mult(p):
    """

    Parámetros
    ----------
    p : ENTERO
        PRIMO

    Devuelve
    --------
    ENTERO
        NEUTRO MULTIPLICATIVO DE Z/pZ

    """
    return 1

def suma(a, b, p):
    """

    Parámetros
    ----------
    a : ENTERO
        ELEMENTO DE Z/pZ
    b : ENTERO
        ELEMENTO DE Z/pZ
    p : ENTERO
        PRIMO

    Devuelve
    --------
    ENTERO
        SUMA DE LOS ELEMENTOS a Y b EN Z/pZ

    """
    return (a + b)%p

def inv_ad(a, p):
    """

    Parámetros
    ----------
    a : ENTERO
        ELEMENTO DE Z/PZ
    p : ENTERO
        PRIMO

    Devuelve
    --------
    ENTERO
        ELEMENTO INVERSO ADITIVO DE a EN Z/pZ

    """
    return (-a)%p

def mult(a, b, p):
    """

    Parámetros
    ----------
    a : ENTERO
        ELEMENTO DE Z/pZ
    b : ENTERO
        ELEMENTO DE Z/pZ
    p : ENTERO
        PRIMO

    Devuelve
    --------
    ENTERO
        PRODUCTO DE LOS ELEMENTOS a Y b EN Z/pZ

    """
    return(a*b)%p

def inv_mult(a, p):
    """

    Parámetros
    ----------
    a : ENTERO
        ELEMENTO DE Z/pZ
    p : ENTERO
        PRIMO

    Devuelve
    --------
    ENTERO
        ELEMENTO INVERSO MULTIPLICATIVO DE a EN Z/pZ

    """
    g,b,c = aux.gcd(a, p)
    return b%p

def potencia(a, r, p):
    """

    Parámetros
    ----------
    a : ENTERO
        ELEMENTO DE Z/pZ
    r : ENTERO
        EXPONENTE AL QUE SE QUIERE ELEVAR EL ELEMENTO a
    p : ENTERO
        PRIMO

    Devuelve
    --------
    x : ENTERO
        ELEMENTO a^r PERTENECIENTE A Z/pZ

    """
    if r<0:
        a = inv_mult(a, p)
        r = -r
    if r == 0:
        x = neutro_mult(p)
    elif r % 2 == 0:
        x = potencia(a, r//2, p)
        x = mult(x, x, p)
    else:
        x = potencia(a, (r-1)//2, p)
        x = mult(x, x, p)
        x = mult(a, x, p)
    return x 

def rand(p):
    """

    Parámetros
    ----------
    p : ENTERO
        PRIMO

    Devuelve
    --------
    ENTERO
        ELEMENTO ALEATORIO DE Z/pZ

    """
    return random.randint(0, p-1)
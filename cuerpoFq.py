# -*- coding: utf-8 -*-
"""
IMPLEMENTACIÓN DEL CUERPO Fq:= (Z/pZ)[x]/<f(x)>

Consideraciones:
    1) f debe ser mónico e irreducible
    2) q = p^n con n = deg(f)
""" 

import anilloZpZx as ZpZx

def neutro_ad(p, f):
    """

    Parámetros
    ----------
    p : ENTERO
        PRIMO
    f : LISTA
        POLINOMIO

    Devuelve
    --------
    LISTA
        NEUTRO ADITIVO DE (Z/pZ)[x]/<f(x)>

    """
    return ZpZx.neutro_ad(p)

def neutro_mult(p, f):
    """

    Parámetros
    ----------
    p : ENTERO
        PRIMO
    f : LISTA
        POLINOMIO

    Devuelve
    --------
    LISTA
        NEUTRO MULTIPLICATIVO DE (Z/pZ)[x]/<f(x)>

    """
    return ZpZx.neutro_mult(p)

def suma(g, h, p, f):
    """

    Parámetros
    ----------
    g : LISTA
        POLINOMIO
    h : LISTA
        POLINOMIO
    p : ENTERO
        PRIMO
    f : LISTA
        POLINOMIO

    Devuelve
    --------
    LISTA
        SUMA DE LOS ELEMENTOS g Y h EN (Z/pZ)[x]/<f>

    """
    return ZpZx.suma(g, h, p)

def inv_ad(g, p, f):
    """

    Parámetros
    ----------
    g : LISTA
        POLINOMIO
    p : ENTERO
        PRIMO
    f : LISTA
        POLINOMIO

    Devuelve
    --------
    LISTA
        ELEMENTO INVERSO ADITIVO DE g EN (Z/pZ)[x]/<f>

    """
    return ZpZx.inv_ad(g, p)

def mult(g, h, p, f):
    """

    Parámetros
    ----------
    g : LISTA
        POLINOMIO
    h : LISTA
        POLINOMIO.
    p : ENTERO
        PRIMO
    f : LISTA
        POLINOMIO

    Devuelve
    --------
    LISTA
        PRODUCTO DE LOS ELEMENTOS g Y h EN Z/pZ[x]

    """
    return ZpZx.div(ZpZx.mult(g, h, p), f, p)[1]

def inv_mult(g, p, f):
    """

    Parámetros
    ----------
    g : LISTA
        POLINOMIO
    p : ENTERO
        PRIMO
    f : LISTA
        POLINOMIO

    Devuelve
    --------
    LISTA
        ELEMENTO INVERSO MULTIPLICATIVO DE g EN (Z/pZ)[x]/<f>

    """
    return ZpZx.div(ZpZx.gcd_ext(f, g, p)[2], f, p)[1]

def potencia(g, r, p, f):
    """

    Parámetros
    ----------
    g : LISTA
        ELEMENTO DE (Z/pZ)[x]/<f>
    r : ENTERO
        EXPONENTE AL QUE SE QUIERE ELEVAR EL ELEMENTO g
    p : ENTERO
        PRIMO
    f : LISTA
        POLINOMIO

    Devuelve
    --------
    x : ENTERO
        ELEMENTO g^r PERTENECIENTE A (Z/pZ)[x]/<f>

    """
    if r<0:
        g = inv_mult(g, p, f)
        r = -r
    if r == 0:
        x = neutro_mult(p, f)
    elif r % 2 == 0:
        x = potencia(g, r//2, p, f)
        x = mult(x, x, p, f)
    else:
        x = potencia(g, (r-1)//2, p, f)
        x = mult(x, x, p, f)
        x = mult(g, x, p, f)
    return x 

def rand(n, p, f):
    """

    Parámetros
    ----------
    n : ENTERO
        GRADO DE UN POLINOMIO (n<deg(f))
    p : ENTERO
        PRIMO
    f : LISTA
        POLINOMIO

    Devuelve
    --------
    LISTA
        POLINOMIO ALEATORIO DE GRADO n DE (Z/pZ)[x]/<f>

    """
    return ZpZx.rand(n, p)
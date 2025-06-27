# -*- coding: utf-8 -*-
"""
IMPLEMENTACIÓN DEL ANILLO Fq[x] = ((Z/pZ)[z]/<f(z)>)[x]

Consideraciones:
    1) f(z) debe ser mónico e irreducible
    2) q = p^n con n = deg(f(z))
""" 

import random
from copy import deepcopy
import auxiliar as aux
import cuerpoFq as Fq

def neutro_ad(p, f):
    """

    Parámetros
    ----------
    p : ENTERO
        PRIMO
    f : LISTA
        POLINOMIO

    Devuelve
    -------
    LISTA
        NEUTRO ADITIVO DE ((Z/pZ)[x]/<f(x)>)[x]

    """
    return []

def neutro_mult(p, f):
    """

    Parámetros
    ----------
    p : ENTERO
        PRIMO
    f : LISTA
        POLINOMIO

    Devuelve
    -------
    LISTA
        NEUTRO MULTIPLICATIVO DE ((Z/pZ)[x]/<f(x)>)[x]
    
    """
    return [Fq.neutro_mult(p, f)]

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
        SUMA DE LOS ELEMENTOS g Y h EN ((Z/pZ)[x]/<f>)[x]

    """
    m = max(len(g), len(h))
    s = [neutro_ad(p, f)]*m
    g = g + [neutro_ad(p, f)]*(m-len(g))
    h = h + [neutro_ad(p, f)]*(m-len(h))
    for i in range(0, m):
        s[i] = Fq.suma(g[i], h[i], p, f)
    return reduce(s, p, f)

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
        ELEMENTO INVERSO ADITIVO DE g EN ((Z/pZ)[x]/<f>)[x]

    """
    inv = [neutro_ad(p, f)]*len(g)
    for i in range(0, len(g)):
        inv[i] = Fq.inv_ad(g[i], p, f)
    return reduce(inv, p, f)

def mult(g, h, p, f):
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
        PRODUCTO DE LOS ELEMENTOS g Y h EN ((Z/pZ)[x]/<f>)[x]

    """
    deg = (len(g)-1)+(len(h)-1)
    m = [neutro_ad(p, f)]*(deg+1)
    for i in range(0, len(g)):
        for j in range(0, len(h)):
            m[i+j] = Fq.suma(Fq.mult(g[i], h[j], p, f), m[i+j], p, f)
    return reduce(m, p, f) 
                        
def div(g, h, p, f):
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
    cociente : LISTA
        COCIENTE DE LA DIVISIÓN DEL POLINOMIO g 
        ENTRE EL POLINOMIO h EN ((Z/pZ)[x]/<f>)[x]
    resto : LISTA
        RESTO DE LA DIVISIÓN DEL POLINOMIO g 
        ENTRE EL POLINOMIO h EN ((Z/pZ)[x]/<f>)[x]

    """
    if h == neutro_ad(p, f):
        raise Exception("No se puede dividir por cero")
    else:
        grado_pol1 = len(g)-1
        grado_pol2 = len(h)-1
        pol1 = deepcopy(g)
        cociente = [Fq.neutro_ad(p, f)]*(grado_pol1-grado_pol2+1)
        while grado_pol1 >= grado_pol2 and pol1 != neutro_ad(p, f):
            pol2 = deepcopy(h)
            a = pol1[len(pol1)-1]
            b = pol2[len(pol2)-1]
            c = Fq.mult(a, Fq.inv_mult(b, p, f), p, f)
            dif = grado_pol1 - grado_pol2
            cociente[dif] = c
            for i in range(0, len(pol2)):
                pol2[i] = Fq.mult(Fq.inv_ad(c, p, f), pol2[i], p, f)
            pol2 = [Fq.neutro_ad(p, f)]*dif + pol2
            pol1 = suma(pol1, pol2, p, f)
            grado_pol1 = len(pol1)-1
        resto = pol1
        return cociente, resto

def gcd(g, h, p, f):
    """

    Parameters
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
        GCD DE LOS POLINOMIOS g Y h en ((Z/pZ)[x]/<f>)[x]

    """
    while h != neutro_ad(p, f):
        cociente, resto = div(g, h, p, f)
        g, h = h, resto
    inv = [Fq.inv_mult(g[len(g)-1], p, f)]
    return mult(inv, g, p, f)

def gcd_ext(g, h, p, f):
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
    gcd : LISTA
        GCD DE LOS POLINOMIOS g Y h en ((Z/pZ)[x]/<f>)[x]
    r : LISTA
    s : LISTA
        POLINOMIOS TALES QUE gcd = g*r + h*s

    """
    r0, s0 = neutro_mult(p, f), neutro_ad(p, f)
    r1, s1 = neutro_ad(p, f), neutro_mult(p, f)
    while h != neutro_ad(p, f):
        cociente, resto = div(g, h, p, f)
        g, h = h, resto
        aux1 = r1 
        r1 = suma(r0, inv_ad(mult(cociente, r1, p, f), p, f), p, f)
        r0 = aux1
        aux2 = s1
        s1 = suma(s0, inv_ad(mult(cociente, s1, p, f), p, f), p, f)
        s0 = aux2
    inv = [Fq.inv_mult(g[len(g)-1], p, f)]
    gcd = mult(inv, g, p, f)
    r = mult(inv, r0, p, f)
    s = mult(inv, s0, p, f)
    return gcd, r, s

def reduce(g, p, f):
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
        POLINOMIO g REPRESENTADO CON UNA LISTA 
        CON ÚNICAMENTE LOS COEFICIENTES NECESARIOS

    """
    if g == neutro_ad(p, f):
        return g
    else:
        if (g[len(g)-1] == Fq.neutro_ad(p, f)):
            g.pop(len(g)-1)
            return reduce(g, p, f)
        else:
            return g

def interpola_lagrange(k, a, b, p, f):
    """

    Parámetros
    ----------
    k : ENTERO
        NÚMERO DE PUNTOS A INTERPOLAR
    a : LISTA
        CONTIENE LAS PRIMERAS COORDENADAS 
        DE LOS PUNTOS DE INTERPOLACIÓN
    b : LISTA
        CONTIENE LAS SEGUNDAS COORDENADAS 
        DE LOS PUNTOS DE INTERPOLACIÓN
    p : ENTERO
        PRIMO

    Devuelve
    --------
    g : LISTA
        POLINOMIO DE INTERPOLACIÓN DE LAGRANGE 
        DE LOS PUNTOS {(a_i, b_i)}_{i}

    """
    g = neutro_ad(p, f)
    for i in range(0, k):
        L_i = neutro_mult(p, f)
        for j in range(0, k):
            if j != i:
                p1 = \
                Fq.inv_mult(Fq.suma(a[i], \
                Fq.inv_ad(a[j], p, f), p, f), p, f)
                p0 = Fq.mult(Fq.inv_ad(a[j], p, f), p1, p, f)
                pol = [p0, p1]
                L_i = mult(L_i, pol, p, f)
        g = suma(g, mult([b[i]], L_i, p, f), p, f)
    return g

def evalua(g, a, p, f):
    """

    Parámetros
    ----------
    g : LISTA
        POLINOMIO
    a : LISTA
        ELEMENTO DE (Z/pZ)[x] DONDE SE QUIERE EVALUAR EL POLINOMIO g
    p : ENTERO
        PRIMO
    f : LISTA
        POLINOMIO

    Devuelve
    --------
    g_a : LISTA
        RESULTADO DE EVALUAR EL POLINOMIO g EN EL ELEMENTO a

    """
    g_a = g[0]
    potencia = Fq.neutro_mult(p, f)
    for i in range(1, len(g)):
        potencia = Fq.mult(potencia, a, p, f)
        g_a = Fq.suma(g_a, Fq.mult(g[i], potencia, p, f), p, f)
    return g_a

def potencia_modulo(g, r, h, p, f):
    """

    Parámetros
    ----------
    g : LISTA
        POLINOMIO
    r : ENTERO, r>=0
        POTENCIA A LA QUE SE QUIERE ELEVAR EL POLINOMIO g
    h : LISTA
        POLINOMIO
    p : ENTERO
        PRIMO
    f : LISTA
        POLINOMIO

    Devuelve
    --------
    a : LISTA
        POLINOMIO a = g^r (mod h) EN ((Z/pZ)[x]/<f>)[x]

    """
    a = neutro_mult(p, f)
    b = g
    while r>0:
        if (r%2) == 1:
            a = div(mult(a, b, p, f), h, p, f)[1]
        b = div(mult(b, b, p, f), h, p, f)[1]
        r = r//2
    return a

def irreducible(h, p, f):
    """

    Parámetros
    ----------
    h : LISTA
        POLINOMIO A COMPROBAR SI ES IRREDUCIBLE
    p : ENTERO
        PRIMO
    f : LISTA
        POLINOMIO

    Devuelve
    --------
    irred : BOOLEANO
        TRUE/FALSE PARA INDICAR SI h ES IRREDUCIBLE SEGÚN 
        EL TEST DE RABIN

    """
    irred = True
    n = len(f)-1
    m = len(h)-1
    ps = aux.calcula_ps(m)
    i = 0
    x = [Fq.neutro_ad(p, f), Fq.neutro_mult(p, f)]
    x_potencia = potencia_modulo(x, (p**n)**m, h, p, f)
    pol = suma(x_potencia,  inv_ad(x, p, f), p, f)
    if div(pol, h, p, f)[1] != neutro_ad(p, f):
        irred = False
    while irred and i<len(ps):
        pi = ps[i]
        aux1 = [Fq.neutro_ad(p, f), Fq.neutro_mult(p, f)]
        for j in range(pi):
            aux1 = potencia_modulo(aux1, p, h, p, f)
        aux2 = \
        suma(aux1, \
        [Fq.neutro_ad(p, f), Fq.inv_ad(Fq.neutro_mult(p, f), p, f)], p, f)
        g = gcd(h, aux2, p, f)
        if g != neutro_mult(p, f):
            irred = False
        i += 1
    return irred

def rand(n, p, f):
    """

    Parameters
    ----------
    n : ENTERO
        GRADO DE UN POLINOMIO
    p : ENTERO
        PRIMO
    f : LISTA
        POLINOMIO

    Devuelve
    --------
    pol : LISTA
        POLINOMIO ALEATORIO DE GRADO n DE ((Z/pZ)[x]/<f>)[x]

    """
    pol = []
    for i in range(n+1):
        deg = random.randint(0, len(f)-2)
        pol.append(Fq.rand(deg, p, f))
    return pol
# -*- coding: utf-8 -*-
"""
IMPLEMENTACIÓN DEL ANILLO (Z/pZ)[x]

Consideraciones: 
    1) El polinomio c0 + c1x + ... + cn*x^n se representa 
    con la lista [c0, c1, c2, ..., cn]
    2) El polinomio nulo se representa con [] 
"""

from copy import deepcopy
import auxiliar as aux
import cuerpoZpZ as ZpZ

def neutro_ad(p):
    """

    Parámetros
    ----------
    p : ENTERO
        PRIMO

    Devuelve
    --------
    LISTA
        NEUTRO ADITIVO DE (Z/pZ)[x]

    """
    return []

def neutro_mult(p):
    """

    Parámetros
    ----------
    p : ENTERO
        PRIMO

    Devuelve
    --------
    LISTA
        NEUTRO MULTIPLICATIVO DE (Z/pZ)[x]

    """
    return [ZpZ.neutro_mult(p)]

def suma(g, h, p):
    """

    Parámetros
    ----------
    g : LIST
        POLINOMIO
    h : LIST
        POLINOMIO
    p : ENTERO
        PRIMO

    Devuelve
    --------
    LISTA
        SUMA DE LOS ELEMENTOS g Y h EN (Z/pZ)[x]

    """
    m = max(len(g), len(h))
    s = [ZpZ.neutro_ad(p)]*m
    g = g + [ZpZ.neutro_ad(p)]*(m-len(g))
    h = h + [ZpZ.neutro_ad(p)]*(m-len(h))
    for i in range(0, m):
        s[i] = ZpZ.suma(g[i], h[i], p)
    return reduce(s, p)

def inv_ad(g, p):
    """

    Parámetros
    ----------
    g : LIST
        POLINOMIO
    p : ENTERO
        PRIMO

    Devuelve
    --------
    inv : LISTA
        ELEMENTO INVERSO ADITIVO DE g EN (Z/pZ)[x]

    """
    inv = [ZpZ.neutro_ad(p)]*len(g)
    for i in range(0, len(g)):
        inv[i] = ZpZ.inv_ad(g[i], p)
    return inv

def mult(g, h, p):
    """

    Parámetros
    ----------
    g : LISTA
        POLINOMIO
    h : LISTA
        POLINOMIO
    p : ENTERO
        PRIMO

    Devuelve
    --------
    LISTA
        PRODUCTO DE LOS ELEMENTOS g Y h EN (Z/pZ)[x]

    """
    deg = (len(g)-1)+(len(h)-1)
    m = [ZpZ.neutro_ad(p)]*(deg+1)
    for i in range(0, len(g)):
        for j in range(0, len(h)):
            m[i+j] = ZpZ.suma(ZpZ.mult(g[i], h[j], p), m[i+j], p)
    return reduce(m, p)

def div(g, h, p):
    """

    Parámetros
    ----------
    g : LISTA
        POLINOMIO
    h : LISTA
        POLINOMIO
    p : ENTERO
        PRIMO

    Devuelve
    --------
    cociente : LISTA
        COCIENTE DE LA DIVISIÓN DEL POLINOMIO g 
        ENTRE EL POLINOMIO h EN (Z/pZ)[x]
    resto : LISTA
        RESTO DE LA DIVISIÓN DEL POLINOMIO g 
        ENTRE EL POLINOMIO h EN (Z/pZ)[x]

    """
    if h == neutro_ad(p):
        raise ValueError("No se puede dividir por cero")
    else:
        grado_pol1 = len(g)-1
        grado_pol2 = len(h)-1
        pol1 = deepcopy(g)
        cociente = [ZpZ.neutro_ad(p)]*(grado_pol1-grado_pol2+1)
        while grado_pol1 >= grado_pol2:
            pol2 = deepcopy(h)
            a = pol1[len(pol1)-1]
            b = pol2[len(pol2)-1]
            c = ZpZ.mult(a, ZpZ.inv_mult(b, p), p)
            dif = grado_pol1 - grado_pol2
            cociente[dif] = c
            for i in range(0, len(pol2)):
                pol2[i] = ZpZ.mult(ZpZ.inv_ad(c, p), pol2[i], p)
            pol2 = [ZpZ.neutro_ad(p)]*dif + pol2
            pol1 = suma(pol1, pol2, p)
            grado_pol1 = len(pol1)-1
        resto = pol1
        return cociente, resto
    
def gcd(g, h, p):
    """

    Parámetros
    ----------
    g : LISTA
        POLINOMIO
    h : LISTA
        POLINOMIO
    p : ENTERO
        PRIMO

    Devuelve
    --------
    LISTA
        GCD DE LOS POLINOMIOS g Y h en (Z/pZ)[x]

    """
    while len(h) >= 1:
        cociente, resto = div(g, h, p)
        g, h = h, resto
    inv = [ZpZ.inv_mult(g[len(g)-1], p)]
    return mult(inv, g, p)

def gcd_ext(g, h, p):
    """

    Parámetros
    ----------
    g : LISTA
        POLINOMIO
    h : LISTA
        POLINOMIO
    p : ENTERO
        PRIMO

    Devuelve
    --------
    gcd : LISTA
        GCD DE LOS POLINOMIOS g Y h en (Z/pZ)[x]
    r : LISTA
    s : LISTA
        POLINOMIOS TALES QUE gcd = g*r + h*s

    """
    r0, s0 = neutro_mult(p), neutro_ad(p)
    r1, s1 = neutro_ad(p), neutro_mult(p)
    while len(h) >= 1:
        cociente, resto = div(g, h, p)
        g, h = h, resto
        aux1 = r1 
        r1 = suma(r0, inv_ad(mult(cociente, r1, p), p), p)
        r0 = aux1
        aux2 = s1
        s1 = suma(s0, inv_ad(mult(cociente, s1, p), p), p)
        s0 = aux2
    inv = [ZpZ.inv_mult(g[len(g)-1], p)]
    gcd = mult(inv, g, p)
    r = mult(inv, r0, p)
    s = mult(inv, s0, p)
    return gcd, r, s

def reduce(g, p):
    """

    Parámetros
    ----------
    g : LISTA
        POLINOMIO
    p : ENTERO
        PRIMO

    Devuelve
    --------
    LISTA
        POLINOMIO g REPRESENTADO CON UNA LISTA 
        CON ÚNICAMENTE LOS COEFICIENTES NECESARIOS

    """
    if g == neutro_ad(p):
        return g
    else:
        if (g[len(g)-1] == ZpZ.neutro_ad(p)):
            g.pop(len(g)-1)
            return reduce(g, p)
        else:
            return g

def interpola_lagrange(k, a, b, p):
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
    g = neutro_ad(p)
    for i in range(0, k):
        l_i = neutro_mult(p)
        for j in range(0, k):
            if j != i:
                p1 = \
                ZpZ.inv_mult(ZpZ.suma(a[i], ZpZ.inv_ad(a[j], p), p), p)
                p0 = ZpZ.mult(ZpZ.inv_ad(a[j], p), p1, p)
                pol = [p0, p1]
                l_i = mult(l_i, pol, p)
        g = suma(g, mult([b[i]], l_i, p), p)
    return g

def evalua(g, a, p):
    """

    Parámetros
    ----------
    g : LISTA
        POLINOMIO
    a : ENTERO
        ELEMENTO DE Z/pZ DONDE SE QUIERE EVALUAR EL POLINOMIO g
    p : ENTERO
        PRIMO

    Devuelve
    --------
    g_a : ENTERO
        RESULTADO DE EVALUAR EL POLINOMIO g EN EL ELEMENTO a

    """
    g_a = g[0]
    potencia = ZpZ.neutro_mult(p)
    for i in range(1, len(g)):
        potencia = ZpZ.mult(potencia, a, p)
        g_a = ZpZ.suma(g_a, ZpZ.mult(g[i], potencia, p), p)
    return g_a

def potencia_modulo(g, r, h, p):
    """

    Parámetros
    ----------
    g : LISTA
        POLINOMIO
    r : ENTERO, r>=0
        POTENCIA A LA QUE SE QUIERE ELEVAR EL POLINOMIO g
    h : LISTA
        MÓDULO QUE SE VA A TOMAR PARA CALCULAR LA POTENCIA
    p : ENTERO
        PRIMO

    Devuelve
    --------
    a : LISTA
        POLINOMIO a = g^r (mod h) EN (Z/pZ)[x]

    """
    a = neutro_mult(p)
    b = g
    while r>0:
        if (r%2) == 1:
            a = div(mult(a, b, p), h, p)[1]
        b = div(mult(b, b, p), h, p)[1]
        r = r//2
    return a

def irreducible(h, p):
    """

    Parameters
    ----------
    h : LISTA
        POLINOMIO A COMPROBAR SI ES IRREDUCIBLE
    p : ENTERO
        PRIMO

    Devuelve
    --------
    irred : BOOLEANO
        TRUE/FALSE PARA INDICAR SI h ES IRREDUCIBLE SEGÚN
        EL TEST DE RABIN

    """
    irred = True
    n = len(h)-1
    ps = aux.calcula_ps(n)
    i = 0
    x = [ZpZ.neutro_ad(p), ZpZ.neutro_mult(p)]
    x_potencia = potencia_modulo(x, p**n, h, p)
    pol = suma(x_potencia,  inv_ad(x, p), p)
    if div(pol, h, p)[1] != neutro_ad(p):
        irred = False
    while irred and i<len(ps):
        pi = ps[i]
        aux1 = [ZpZ.neutro_ad(p), ZpZ.neutro_mult(p)]
        for j in range(pi):
            aux1 = potencia_modulo(aux1, p, h, p)
        aux2 = \
        suma(aux1, \
        [ZpZ.neutro_ad(p), ZpZ.inv_ad(ZpZ.neutro_mult(p), p)], p)
        g = gcd(h, aux2, p)
        if g != neutro_mult(p):
            irred = False
        i += 1
    return irred

def rand(n, p):
    """

    Parámetros
    ----------
    n : ENTERO
        GRADO DE UN POLINOMIO
    p : ENTERO
        PRIMO

    Devuelve
    --------
    pol : LISTA
        POLINOMIO ALEATORIO DE GRADO n DE (Z/pZ)[x]

    """
    pol = []
    for i in range(n+1):
        pol.append(ZpZ.rand(p))
    return pol
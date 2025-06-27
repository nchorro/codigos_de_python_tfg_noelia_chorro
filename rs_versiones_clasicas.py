# -*- coding: utf-8 -*-
"""
CÓDIGOS DE REED-SOLOMON - VERSIONES CLÁSICAS
"""

import cuerpoFq as Fq
import anilloFqx as Fqx
import matricesFq as mat
from math import ceil 

def rs_algebra_lineal(n, k, a, y, p, f):
    """

    Parámetros
    ----------
    n : ENTERO
        NÚMERO DE PUNTOS DE EVALUACIÓN a_i
    k : ENTERO, k<= n-2
        COTA SUPERIOR DEL GRADO DEL POLINOMIO s(x) 
        CUYOS COEFICIENTES CONSTITUYEN EL MENSAJE ORIGINAL
    a : LISTA
        CONTIENE LOS PUNTOS DE EVALUACIÓN a_1, a_2, ..., a_n
    y : LISTA
        CONTIENE LOS VALORES RECIBIDOS A TRAVÉS DEL 
        CANAL CON RUIDO DE LAS EVALUACIONES
        DE LOS a_i EN s(x)
    p : ENTERO
        PRIMO

    Devuelve
    --------
    s : LISTA
        POLINOMIO s(x) QUE CONSTITUYE LA PALABRA CLAVE 
        QUE SE QUIERE ENVIAR

    """
    
    # Número máximo de errores
    t = (n-k)//2
    
    # Grados de los polinomios e(x) y h(x) buscados
    deg_e = t
    deg_h = t-1 
    
    # Calculamos el polinomio de interpolación de Lagrange
    g = Fqx.interpola_lagrange(n, a, y, p, f)
    deg_g = len(g)-1  
    
    # Calculamos el polinomio w(x)
    w = calcula_w(a, p, f)
    deg_w = len(w)-1
    
    # Construimos la matriz del sistema lineal homogéneo
    filas = max(deg_g + deg_e, deg_w + deg_h)+1
    deg_u = k + ((n-k+1)//2)-1 
    ec = filas-(deg_u+1)
    sist = mat.neutro_ad(ec, deg_e + deg_h + 2, p, f)
    for i in range(0, deg_e + 1):
        col = [Fqx.neutro_ad(p, f)]*i+ g
        while len(col)<filas:
            col.append(Fq.neutro_ad(p, f))
        while len(col)>filas:
            col.pop(len(col)-1)
        sist[i] = col[-ec:]
    for j in range(deg_e + 1, deg_e + deg_h + 2):
        col = [Fqx.neutro_ad(p, f)]*(j - deg_e - 1) + w
        while len(col)<filas:
            col.append(Fq.neutro_ad(p, f))
        while len(col)>filas:
            col.pop(len(col)-1)
        sist[j] = col[-ec:]
        
    # Aplicamos la eliminación Gaussiana de la matriz del sistema
    A = mat.eliminacion_gaussiana(sist, p, f)
    
    # Obtenemos una solución del sistema
    coefs = Fqx.neutro_mult(p, f)*(deg_e + deg_h + 2)
    i = 0
    while i<ec:
        c = Fq.neutro_ad(p, f)
        encontrado = False
        pos = 0
        j = 0
        while j<(deg_e + deg_h + 2):
            if encontrado:
                c = Fq.suma(c,Fq.inv_ad(A[j][i], p, f), p, f)
            if A[j][i] == Fq.neutro_mult(p, f) and not encontrado:
                encontrado = True
                pos = j
            j += 1
        if encontrado:
            coefs[pos] = c
        i += 1
    e = Fqx.reduce(coefs[:deg_e+1], p, f)
    h = Fqx.reduce(coefs[-deg_h-1:], p, f)
    
    # Con los polinomios e(x) y h(x), calculamos u(x) y recuperamos s(x)
    wh = Fqx.mult(w, h, p, f)
    ge = Fqx.mult(g, e, p, f)
    u = Fqx.suma(ge, wh, p, f)
    s = Fqx.div(u, e, p, f)[0]
    
    return s

def rs_gcd_ext_euclideo(n, k, a, y, p, f):
    """

    Parámetros
    ----------
    n : ENTERO
        NÚMERO DE PUNTOS DE EVALUACIÓN a_i
    k : ENTERO
        COTA SUPERIOR DEL GRADO DEL POLINOMIO s(x) 
        CUYOS COEFICIENTES CONSTITUYEN EL MENSAJE ORIGINAL
    a : LISTA
        CONTIENE LOS PUNTOS DE EVALUACIÓN a_1, a_2, ..., a_n
    y : LISTA
        CONTIENE LOS VALORES RECIBIDOS A TRAVÉS DEL 
        CANAL CON RUIDO DE LAS EVALUACIONES
        DE LOS a_i EN s(x)
    p : ENTERO
        PRIMO

    Devuelve
    --------
    s : LISTA
        POLINOMIO s(x) QUE CONSTITUYE LA PALABRA CLAVE 
        QUE SE QUIERE ENVIAR

    """
    
    # Grado del polinomio u(x) buscado
    deg_u = k + ((n-k+1)//2)-1
    
    # Calculamos el polinomio de interpolación de Lagrange
    g = Fqx.interpola_lagrange(n, a, y, p, f)
    
    # Calculamos el polinomio w(x)
    w = calcula_w(a, p, f)
    
    # Inicialización del algoritmo
    h0 = Fqx.neutro_mult(p, f)
    e0 = Fqx.neutro_ad(p, f)
    h1 = Fqx.neutro_ad(p, f)
    e1 = Fqx.neutro_mult(p, f)
    u0 = w
    u1 = g
    
    # Bucle principal
    while (len(u1)-1) > deg_u:
        
        cociente, resto = Fqx.div(u0, u1, p, f)
        
        u0 = u1
        u1 = resto 
        
        aux_h = \
        Fqx.suma(h0, Fqx.inv_ad(Fqx.mult(cociente, h1, p, f), p, f), p, f)
        h0 = h1
        h1 = aux_h
        
        aux_e = \
        Fqx.suma(e0, Fqx.inv_ad(Fqx.mult(cociente, e1, p, f), p, f), p, f)
        e0 = e1
        e1 = aux_e 
    
    # Cuando termina el bucle, se recupera s(x)    
    s = Fqx.div(u1, e1, p, f)[0] 
    
    return s

def calcula_w(a, p, f):
    """

    Parámetros
    ----------
    a : LISTA
        CONTIENE LOS PUNTOS DE EVALUACIÓN a_1, a_2, ..., a_n
    p : ENTERO
        PRIMO
    f : LISTA
        POLINOMIO

    Devuelve
    --------
    w : LISTA
        POLINOMIO PRODUCTO DE TODOS LOS TÉRMINOS DE LA FORMA (x-a_i)

    """
    w = Fqx.neutro_mult(p, f)
    
    for i in range(len(a)):
        p_i = [Fq.inv_ad(a[i], p, f), Fq.neutro_mult(p, f)]
        w = Fqx.mult(w, p_i, p, f)
        
    return w  
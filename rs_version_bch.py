# -*- coding: utf-8 -*-
"""
CÓDIGOS DE REED-SOLOMON - VERSIÓN BCH
"""

import cuerpoFq as Fq
import anilloFqx as Fqx
import matricesFq as mat
from copy import deepcopy

def calcula_Mt_vt(k, n, alpha, i, r, t, p, f):
    """

    Parámetros
    ----------
    k : ENTERO
        DIMENSIÓN DEL CÓDIGO BCH
    n : ENTERO
        LONGITUD DEL CÓDIGO BCH
    alpha : LISTA
        GENERADOR DE Fq*
    i : ENTERO
        POTENCIA
    r : LISTA
        POLINOMIO RECIBIDO POR EL CANAL CON RUIDO
    t : ENTERO
        NÚMERO MÁXIMO DE ERRORES DE TRANSMISIÓN
    p : ENTERO
        PRIMO
    f : LISTA
        POLINOMIO

    Devuelve
    --------
    Mt : LISTA
        MATRIZ DEL SISTEMA PARA HALLAR g
    vt : LISTA
        VECTOR DE TÉRMINOS INDEPENDIENTES DEL SISTEMA PARA HALLAR g
    t : ENTERO
        NÚMERO DE ERRORES OCURRIDOS

    """
    
    # Calculamos la potencia alpha^i
    alpha_i = Fq.potencia(alpha, i, p, f)
    
    # Calculamos simultáneamente las potencias de alpha 
    # y las evaluamos en el polinomio r(x), lo cual es 
    # equivalente a evaluar en el polinomio e(x)
    e_alphas = [Fqx.evalua(r, alpha_i, p, f)]
    for j in range(1, 2*t):
        alpha_i = Fq.mult(alpha_i, alpha, p, f)
        e_alphas.append(Fqx.evalua(r, alpha_i, p, f))
        
    # Construimos la matriz M_t y el vector v_t
    Mt = mat.neutro_ad(t, t, p, f)
    for j in range(0, t):
        Mt[j] = e_alphas[j:t+j]
    vt = e_alphas[-t:]
    
    # Si la matriz tiene determinante nulo, bajamos t 
    # una unidad hasta encontrar el adecuado
    if mat.det(Mt, p, f) == Fq.neutro_ad(p, f) and t!= 0:
        t = t-1
        Mt, vt, t = calcula_Mt_vt(k, n, alpha, i, r, t, p, f)
        
    return Mt, vt, t

def rs_bch(k, n, alpha, i, g, r, p, f):
    """

    Parámetros
    ----------
    k : ENTERO
        DIMENSIÓN DEL CÓDIGO BCH
    n : ENTERO
        LONGITUD DEL CÓDIGO BCH
    alpha : ENTERO
        GENERADOR DE Fq*
    i : ENTERO
        POTENCIA
    r : LISTA
        POLINOMIO RECIBIDO POR EL CANAL CON RUIDO
    p : ENTERO
        PRIMO
    f : LISTA
        POLINOMIO

    Devuelve
    -------
    s : LISTA
        POLINOMIO s(x) QUE CONSTITUYE LA PALABRA CLAVE 
        QUE SE QUIERE ENVIAR

    """
    
    # Valor de t inicial, tomamos el máximo
    t0 = (n-k)//2
    Mt, vt, t = calcula_Mt_vt(k, n, alpha, i, r, t0, p, f)
    
    # Aplicamos la eliminación gaussiana
    A = mat.eliminacion_gaussiana(Mt+[vt], p, f)
    new_Mt = A[:t]
    new_vt = A[len(A)-1]
    
    # Calculamos la solución del sistema que son
    # los coeficientes [w0, w1, ..., w_{t-1}]
    coefs = deepcopy(new_vt)
    
    for j in range(0, len(new_Mt)):
        for s in range(0, len(new_Mt[0])):
            if new_Mt[j][s] == Fq.neutro_mult(p, f):
                coefs[s] = new_vt[j]
                
    # Montamos w añadiendo el coeficiente principal: 
    # w(x) = x^{t} -w_{t-1}x^{t-1} - ... - w_0
    w = Fqx.inv_ad(coefs, p, f) + Fqx.neutro_mult(p, f)
    
    # Hallamos las raíces de w(x) que son los beta_{i}
    betas = []
    if (Fqx.evalua(w, Fq.neutro_mult(p, f), p, f) == Fqx.neutro_ad(p, f)):
        betas.append(0)
    beta = 1
    candidato = Fq.neutro_mult(p, f)
    while beta<n:
        candidato = Fq.mult(candidato, alpha, p, f)
        if (Fqx.evalua(w, candidato, p, f) == Fqx.neutro_ad(p, f)):
            betas.append(beta)
        beta += 1
    
    # Hallamos el vector del sistema para calcular el polinomio error
    alpha_i = Fq.potencia(alpha, i, p, f)
    aux = alpha_i
    r_alpha = [Fqx.evalua(r, alpha_i, p, f)]
    for j in range(t-1):
        aux = Fq.mult(aux, alpha, p, f)
        r_alpha.append(Fqx.evalua(r, aux, p, f))

    # Hallamos la matriz del sistema para calcular el polinomio error
    aux2 = []
    for b in betas:
        e1 = Fq.potencia(alpha_i, b, p, f)
        col = [e1]
        for j in range(1, t):
            e1 = Fq.mult(e1, Fq.potencia(alpha, b, p, f), p, f)
            col.append(e1)
        aux2.append(col)
    
    # Hacemos la eliminación gaussiana
    B = mat.eliminacion_gaussiana(aux2 + [r_alpha], p, f)
    new_B = B[:t]
    new_r = B[len(B)-1]
    
    # Obtenemos el polinomio e
    e_coefs = deepcopy(new_r)
    for j in range(0, len(new_B)):
        for s in range(0, len(new_B[0])):
            if new_B[j][s] == Fq.neutro_mult(p, f):
                e_coefs[s] = new_r[j]
    if len(betas) == 0:
        e = Fqx.neutro_ad(p, f)
    else:
        e = [Fq.neutro_ad(p, f)]*(betas[len(betas)-1]+1)
        for k in range(len(betas)):
            e[betas[k]] = e_coefs[k]
    
    #Recuperamos c y, con él, recuperamos s
    c = Fqx.suma(r, Fqx.inv_ad(e, p, f), p, f)
    s = Fqx.div(c, g, p, f)[0]
    
    return s

def calcula_g(k, n, alpha, i, p, f):
    """

    Parámetros
    ----------
    k : ENTERO
        DIMENSIÓN DEL CÓDIGO BCH
    n : ENTERO
        LONGITUD DEL CÓDIGO BCH
    alpha : LISTA
        GENERADOR DE Fq*
    i : ENTERO
        POTENCIA
    p : ENTERO
        PRIMO
    f : LISTA
        POLINOMIO

    Returns
    -------
    g : LISTA
        CALCULA EL POLINOMIO GENERADOR DEL CÓDIGO

    """
    
    alpha_ij = [Fq.potencia(alpha, i, p, f)]
    
    g = \
    Fqx.suma([Fq.neutro_ad(p, f), Fq.neutro_mult(p, f)], \
    Fqx.inv_ad(alpha_ij, p, f), p, f)
    for j in range(1, n-k):
        alpha_ij = Fqx.mult(alpha_ij, [alpha], p, f)
        a = Fqx.suma([Fq.neutro_ad(p, f), Fq.neutro_mult(p, f)], \
        Fqx.inv_ad(alpha_ij, p, f), p, f)
        g = Fqx.mult(g, a, p, f)
    
    return g
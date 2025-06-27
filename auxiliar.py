# -*- coding: utf-8 -*-
"""
AUXILIAR
"""

def gcd(x,y):
    """

    Parámetros
    ----------
    x : ENTERO
    y : ENTERO

    Devuelve
    --------
    m : ENTERO
    a : ENTERO
    b : ENTERO
    
    m = gcd(x, y)
    m = x*a + y*b

    """
    xespar = x%2 == 0
    yespar = y%2 == 0
    if x < 0:
        m,a,b = gcd(-x,y)
        a = -a
    elif y < 0:
        m,a,b = gcd(x,-y)
        b = -b
    elif x == 0:               
        m,a,b = y,0,1
    elif y == 0:             
        m,a,b = x,1,0
    elif xespar and yespar:
        m,a,b = gcd(x//2, y//2)
        m *= 2
    elif xespar:
        m,a,b = gcd(x//2, y)
        if a % 2 == 0:
           a //= 2
        else:
           a = (a+y)//2
           b = b - x//2
    elif yespar:
        m,a,b = gcd(x, y//2)
        if b % 2 == 0:
           b //= 2
        else:
           b = (b+x)//2
           a = a - y//2
    elif x > y:
        m,a,b = gcd(y, x-y)
        a,b = b,a-b
    else:
        m,a,b = gcd(x, y-x)
        a -= b
    return m,a,b

def es_primo(n):
    """
    
    Parámetros
    ----------
    n : ENTERO

    Devuelve
    --------
    TRUE: SI n ES PRIMO
    FALSE: SI n ES COMPUESTO
    
    """
    
    if n==2 or n==3:
        return True
    
    elif n%2==0:
        return False
    
    else:
        i=1
        while i*i<=n:
            i+=2
            if n%i==0:
                return False
    return True

def calcula_ps(m):
    """

    Parámetros
    ----------
    m : ENTERO

    Devuelve
    --------
    ps : LISTA
        CONTIENE TODOS LOS DIVISORES PRIMOS DE m 

    """
    ps = [1]
    for i in range(2, m):
        if es_primo(i) and (m%i == 0):
            ps.append(m//i)
    return ps




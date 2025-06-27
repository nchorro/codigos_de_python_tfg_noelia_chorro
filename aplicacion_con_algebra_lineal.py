# -*- coding: utf-8 -*-
"""
GRABACIÓN Y LECTURA DE DATOS EN UN CD USANDO CÓDIGOS DE REED-SOLOMON
CON ALGORITMO DE DECODIFICACIÓN DE ÁLGEBRA LINEAL
"""

import anilloFqx as Fqx
import rs_versiones_clasicas as rs
from copy import deepcopy
import random


# FUNCIONES PARA EL ENTRELAZADO Y DESENTRELAZADO DE DATOS

def entrelazado(bloques):
    e = []
    for k in range(len(bloques[0])):
        for j in range(len(bloques)):
            e.append(bloques[j][k])
    ent = []
    for t in range(0, len(e), len(bloques[0])):
        ent.append(e[t: t+len(bloques[0])])
    return ent

def desentrelazado(bloques):
    e = []
    for k in range(len(bloques)):
        for j in range(len(bloques[0])):
            e.append(bloques[k][j])
    ent = []
    for t in range(0, len(bloques)):
        x = []
        for s in range(0, len(bloques[0])):
            x.append(e[t + s*len(bloques)])
        ent.append(x)
    return ent
    

# CONFIGURACIÓN: PARÁMETROS

# Cuerpo finito
p = 2

# Polinomio irreducible
f = [1, 0, 1, 1, 1, 0, 0, 0, 1]

# Capa 1: Código de Reed-Solomon 1
n1 = 32
k1 = 28

# Capa 2: Código de Reed-Solomon 2
n2 = 28
k2 = 24


# Elementos aleatorios de Fq para la evaluación
a1 = []
while len(a1) < n1:
    elem = [random.randint(0,1) for i in range(8)]
    if elem not in a1:
        a1.append(elem)
        
a2 = []
while len(a2) < n2:
    elem = [random.randint(0,1) for i in range(8)]
    if elem not in a2:
        a2.append(elem)


# Queremos grabar en un CD el siguiente mensaje
# mensaje = b"Codigos de Reed-Solomon." * 109  # 24 bytes
r = 13
mensaje = b"Codigos de Reed-Solomon." * r  # 24 bytes


# MENSAJE A GRABAR EN EL CD

with open("intro.txt", "r", encoding="utf-8") as file:
    texto = file.read()

mensaje = texto.encode(encoding="utf-8")[:k2*109]

# Primero convertimos el mensaje a formato binario
binario = []
for b in mensaje:
    bits = []
    for w in range(8):
        bits.append(b%2)
        b //= 2
    binario.append(bits)

# Lo transformamos en 109 listas de 24 bytes
s2 = [binario[j:j+k2] for j in range(0,len(binario),k2)]


# PROCESO DE GRABACIÓN DE DATOS EN EL CD

# Codificamos con RS en la capa 2 (codificamos cada bloque de 24 bytes)
y2 = []
for j in range(len(s2)):
    y2_j = []
    for k in range(len(a2)):
        y2_j.append(Fqx.evalua(s2[j], a2[k], p, f))
    y2.append(y2_j)

    
# Aplicamos el entrelazado
s1 = entrelazado(y2)

# Codificamos con RS en la capa 1
y1 = []
for j in range(len(s1)):
    y1_j = []
    for k in range(len(a1)):
        y1_j.append(Fqx.evalua(s1[j], a1[k], p, f))
    y1.append(y1_j)

# Se graba c1 en el CD Y ocurren los errores
r1 = deepcopy(y1)
r1[101] = [[random.randint(0,1) for x in range(8)] for y in range(32)]
r1[102] = [[random.randint(0,1) for x in range(8)] for y in range(32)]
r1[103] = [[random.randint(0,1) for x in range(8)] for y in range(32)]
r1[104] = [[random.randint(0,1) for x in range(8)] for y in range(32)]
r1[105] = [[random.randint(0,1) for x in range(8)] for y in range(32)]
r1[106] = [[random.randint(0,1) for x in range(8)] for y in range(32)]
r1[107] = [[random.randint(0,1) for x in range(8)] for y in range(32)]


# PROCESO DE LECTURA DE DATOS DEL CD

# Decodificar con RS en la capa 1
s1_rs = []
for j in range(len(r1)):
    s1_rs_j = rs.rs_gcd_ext_euclideo(n1, k1, a1, r1[j], p, f)
    if len(s1_rs_j) != k1:
        s1_rs_j = s1_rs_j + [[]]*(k1-len(s1_rs_j))
    s1_rs.append(s1_rs_j)

# Desentrelazado
r2 = desentrelazado(s1_rs)

# Decodificar con RS en la capa 2
s2_rs = []
for j in range(len(r2)):
    s2_rs_j = rs.rs_gcd_ext_euclideo(n2, k2, a2, r2[j], p, f)
    if len(s2_rs_j) != k2:
        s2_rs_j = s2_rs_j + [[]]*(k2-len(s2_rs_j))
    s2_rs.append(s2_rs_j)

# Añadimos los bits = 0, eliminados por redundancia
for j in range(109):
    for k in range(k2):
        s2_rs[j][k] += [0] * (8-len(s2_rs[j][k]))

# Convertimos el mensaje en formato binario en caracteres de texto
msg_ints = [0] * (109*k2)
for j in range(109):
    for k in range(24):
        for w in range(8):
            msg_ints[j*k2+k] += (2 ** w) * s2_rs[j][k][w]
msg_bytes = bytes(msg_ints)

# Mensaje leido en el CD
mensaje_leido = msg_bytes.decode(encoding="utf-8")
print('El mensaje leido en el CD es : ', mensaje_leido)
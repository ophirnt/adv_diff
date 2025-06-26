#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 01:56:05 2019

@author: ophir
"""

import numpy as np
import matplotlib.pyplot as plt

def elem(epsilon, beta, h, zeta):
    ''' Constrói matriz do elemento '''
    
    return epsilon/h * np.array([[1,-1],[-1,1]]) + \
           beta/2 * np.array([[-1, 1],[-1,1]]) + \
           beta*zeta/2 * np.array([[1,-1],[-1,1]])
           
def matriz_elementos(elem, n):
    ''' Constrói a matriz com todos os elementos da malha '''
    
    a = []
    
    for i in range(n):
        a.append(elem)
    
    a = np.array(a)
    
    return a


# Conectividade
ien = np.array([[1,2],[2,3],[3,4],[4,5],[5,6],[6,7],[7,8],[8,9],[9,10],[10,11]]) \
     - np.array([[1,1],[1,1],[1,1],[1,1],[1,1],[1,1],[1,1],[1,1],[1,1],[1,1]])

# Matrizes globais

A = np.zeros([11,11])
A_upwind = np.zeros([11,11])
A_exato = np.zeros([11,11])

# Valores do problema

h = 1/10
epsilon = 0.5
beta = 1#10
Peclet = beta*h/(2*epsilon)


# Construção da matriz dos elementos
a = matriz_elementos(elem(epsilon, beta, h, 0),10)
a_upwind = matriz_elementos(elem(epsilon, beta, h, 1),10)
a_exato = matriz_elementos(elem(epsilon, beta, h, 1/np.tanh(Peclet) - 1/Peclet),10)

# Montagem das matrizes globais
for k in range(10):
    for I in range(2):
        for J in range(2):
            
            A[ien[k,I], ien[k,J]] = A[ien[k,I], ien[k,J]] + a[k,I,J]
            A_upwind[ien[k,I], ien[k,J]] = A_upwind[ien[k,I], ien[k,J]] + a_upwind[k,I,J]
            A_exato[ien[k,I], ien[k,J]] = A_exato[ien[k,I], ien[k,J]] + a_exato[k,I,J]


# Construção do vetor F
            
F = np.zeros(9)

# Aplicação das condições de Dirichlet

A = A[1:10]
A_upwind = A_upwind[1:10]
A_exato = A_exato[1:10]

F_inst = F -A[:,10]*1 -A[:,0]*0
F_upwind = F -A_upwind[:,10]*1 -A_upwind[:,0]*0
F_exato = F -A_exato[:,10]*1 -A_exato[:,0]*0

A = A[:,1:10]
A_upwind = A_upwind[:,1:10]
A_exato = A_exato[:,1:10]

# Solução exata

pos = np.linspace(0, 1, 100)
sol_exata = lambda x: (np.e**(beta/epsilon * x) -1)/(np.e**(beta/epsilon) -1)
solucao = sol_exata(pos)

# Resolução dos problemas discretos

u_inst = np.linalg.solve(A,F_inst)
u_upwind = np.linalg.solve(A_upwind,F_upwind)
u_exato = np.linalg.solve(A_exato,F_exato)

# Adicionando os valores conhecidos aos vetores solução
u_inst = np.append(u_inst, np.array([1.]))
u_inst= np.append(np.array([0.]), u_inst)

u_upwind = np.append(u_upwind, np.array([1.]))
u_upwind = np.append(np.array([0.]), u_upwind)

u_exato = np.append(u_exato, np.array([1.]))
u_exato = np.append(np.array([0.]), u_exato)

# Gráficos

print("WAH")
# print(sol_exata(np.array([0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0])))
print(sol_exata(np.array([0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.0])))

nodes = np.linspace(0, 1, 11)

plt.subplot(2,2,1)

plt.plot(pos,solucao, color = 'black')
plt.plot(nodes, u_inst, linestyle = '--', marker='x', color='black')
plt.plot(nodes, u_upwind, linestyle = '-.', marker='^', color='black')
plt.plot(nodes, u_exato, linestyle = ':', marker='o', color='black')
plt.grid(True)
plt.legend( ('Solução analítica', 'Não-estabilizada', 'Upwind', 'Nodalmente exata') )
plt.title('Todas as soluções')
plt.xlabel('Posição em x')
plt.ylabel('Quantidade da propriedade u')

plt.subplot(2,2,2)

plt.plot(pos,solucao, color = 'black')
plt.plot(nodes, u_inst, linestyle = '--', marker='x', color='black')
plt.grid(True)
plt.legend( ('Solução analítica', 'Não-estabilizada') )
plt.title('Resultado sem estabilização')
plt.xlabel('Posição em x')
plt.ylabel('Quantidade da propriedade u')

plt.subplot(2,2,3)

plt.plot(pos,solucao, color = 'black')
plt.plot(nodes, u_upwind, linestyle = '-.', marker='^', color='black')
plt.grid(True)
plt.legend( ('Solução analítica', 'Upwind') )
plt.title('Resultado com upwind')
plt.xlabel('Posição em x')
plt.ylabel('Quantidade da propriedade u')

plt.subplot(2,2,4)

plt.plot(pos,solucao, color = 'black')
plt.plot(nodes, u_exato, linestyle = ':', marker='o', color='black')
plt.grid(True)
plt.legend( ('Solução analítica', 'Nodalmente exata') )
plt.title('Resultado nodalmente exato')
plt.xlabel('Posição em x')
plt.ylabel('Quantidade da propriedade u')

plt.tight_layout()
plt.show()

            
  
            



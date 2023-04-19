# On s'intérèsse au problème de N corps. C'est un problèmes de mécanique céleste
# qui consiste à calucler la trajectoires de N corps s'attirant mutuellement.

# Dans ce projet nous nous intéresserons unique à l'orbite d'une comète tournant autour
# du soleil. Le soleil est supposé fixe. Les trajectoires de la comète se trouvent sur un plan

import numpy as np
import matplotlib.pyplot as plt

from schemas_numerique import *

def schema_euler_implicite(Y_init, t_init, T, N, f):
    """On cherche une approximation d'une équation sous la forme : y' = f(t,y) 
    en se servant du schéma d'euler implicite

    Args:
        y_init (tuple): position initial
        t_init (int)  : temps à partir duquel on s'intérèsse à la solution
        T (int)       : période
        N (int)       : nombre de points de la subdivision de l'interlave [t_init, t_init + T]
        f (function)  : la fonction que l'on souhaite approximer

    Returns:
        function : solution approximative de l'équation y' = f(t,y)
    """
    
    # Initialisation
    Y = [[Y_init, (1, 1)]] # Y = ((q1, q2), (q1', q2')) = (Y1, Y2) système d'ordre 1
    delta_t = T / N # pas d'une subdivision uniforme
    subdivision = np.arange(t_init, t_init + T, delta_t) # subdivision de l'intervalle [t_init, t_init + T]
    
    # Itération : Y_{n+1} = Y_n + delta_t * f(t, Y_N), ie celle du schema d'euler implicite
    for t in subdivision:
        f_Yn = f(t, Y[-1][0], Y[-1][1])
        iteration = [[Y[-1][i][j] + delta_t * f_Yn[i][j] for i in [0, 1]] for j in [0, 1]] # dans R4
        Y.append(iteration) 
    
    return ([elmt[0] for elmt in Y], [elmt[1] for elmt in Y])

if __name__ == "__main__":
    # Constantes
    m = 1
    G = 1
    
    # Système des accélérations des coordonnées q1 et q2 dans le plan orbital
    deuxieme_loi_Newton = lambda t, q1, q2: (-m * G * q1 / (np.sqrt(q1**2 + q2**2)), -m * G * q2 / (np.sqrt(q1 ** 2 + q2 ** 2)))
    
    # On considère l'équation Y'=f(t, Y1, Y2) où Y1, Y2 \in \R^2 et Y = (y, y') et on cherche à trouver y \in \R2
    f = lambda t, Y1, Y2: (Y2, deuxieme_loi_Newton(t, Y1[0], Y1[1]))
    y, y_prime = schema_euler_implicite((1, 1), 0, 60, 100, f)
    q1 = [elmt[0] for elmt in y]
    q2 = [elmt[1] for elmt in y]
    
    # Graphique de la trajctoires de la comète
    plt.plot(q1, q2)
    plt.title("Trajectoires de la comète")
    plt.tight_layout()
    plt.show()
    
    # Équation de l'énergie de la comète
    E = lambda t, q1, q2: (m/2) * (q1(t+1) ** 2 + q2(t+1) ** 2) - G * m ** 2 / np.sqrt((q1(t) ** 2 + q2(t) ** 2))
    
    # Il faudra penser a justifier le choix du plan orbital et non de R3
# On s'intérèsse au problème de N corps. C'est un problèmes de mécanique céleste
# qui consiste à calucler la trajectoires de N corps s'attirant mutuellement.

# Dans ce projet nous nous intéresserons unique à l'orbite d'une comète tournant autour
# du soleil. Le soleil est supposé fixe. Les trajectoires de la comète se trouvent sur un plan

import numpy as np
import matplotlib.pyplot as plt

def euler_explicite(T, N, y_init, f):
    # Initialisation
    y = np.zeros((len(y_init), N+1))
    y[:,0] = y_init
    h = T / N
    t = np.linspace(0, T, num=N+1)
    
    print(y)
    
    # Itération
    for n in range(N):
        y[:, n+1] = (y[:,n] + h * f(t[n], y[:,n]))
        
    return t, y

if __name__ == "__main__":
    # Constantes
    m = 1 # masse du soleil
    G = 1 # Constante gravitationnelle
    T = 1
    N = 1
    y_init = np.array([3, -4, -2, -1])
    
    # On considère l'équation Y'=f(t, Y1, Y2) où Y1, Y2 \in \R^2 et Y = (y, y') et on cherche à trouver y \in \R^2
    f = lambda t, q: np.array([q[2], q[3], 
        -m * G * q[0] / ((q[0]**2 + q[1]**2)**(3/2)), 
        -m * G * q[1] / ((q[0]**2 + q[1]**2)**(3/2))])
    subdivision, approx = euler_explicite(T, N, y_init, f)
    
    # Graphique de la trajctoires de la comète
    plt.plot(approx[0], approx[1])
    plt.title("Trajectoires de la comète")
    plt.tight_layout()
    plt.show()
    
    # Équation de l'énergie de la comète
    E = lambda t, q1, q2, q3, q4: (m/2) * (q3 ** 2 + q4 ** 2) - G * m ** 2 / np.sqrt((q1 ** 2 + q2 ** 2))
    plt.plot(subdivision, [E(None, approx[0][t], approx[1][t], approx[2][t], approx[3][t]) for t in range(len(subdivision))])
    plt.show()
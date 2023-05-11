# On s'intérèsse au problème de N corps. C'est un problèmes de mécanique céleste
# qui consiste à calucler la trajectoires de N corps s'attirant mutuellement.

# Dans ce projet nous nous intéresserons unique à l'orbite d'une comète tournant autour
# du soleil. Le soleil est supposé fixe. Les trajectoires de la comète se trouvent sur un plan

import numpy as np
import matplotlib.pyplot as plt

def schema_euler_explicite(T, N, y_init, f):
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

def schema_semi_implicite(T, N, y_init, f):
    # Initialisation
    y = np.zeros((len(y_init), N+1))
    y[:,0] = y_init
    h = T / N
    t = np.linspace(0, T, num=N+1)
    
    # Itération
    for n in range(N):
        y[:, n+1] = y[:,n] / (1 -  h * f(t[n], y[:,n]))
        
    return t, y

def erreur(exacte, approx):
    return 

if __name__ == "__main__":
    # Constantes
    m = 1 # masse du soleil
    G = 1 # Constante gravitationnelle
    T = 25
    N = 2500
    y_init = np.array([1., 1., 0., 1.])
    
    # On considère l'équation Y'=f(t, Y1, Y2) où Y1, Y2 \in \R^2 et Y = (y, y') et on cherche à trouver y \in \R^2
    f = lambda t, q: np.array([q[2], q[3], 
        -m * G * (q[0] / ((q[0]**2 + q[1]**2)**(3/2))), 
        -m * G * (q[1] / ((q[0]**2 + q[1]**2)**(3/2)))])
    subdivision, approx_euler = schema_euler_explicite(T, N, y_init, f)
    _, approx_semi_implicite = schema_semi_implicite(T, N, y_init, f)
    
    # Représentations graphiques
    fig, ax = plt.subplots(3)
    
    # Graphique de la trajctoires de la comète
    ax[0].plot(approx_euler[0], approx_euler[1], label="euler explicite")
    # ax[0].plot(approx_semi_implicite[0], approx_semi_implicite[1], label="semi_implicite")
    ax[0].set_title("Trajectoires de la comète")
    ax[0].set_xlabel("position sur l'axe q1")
    ax[0].set_ylabel("position sur l'axe q2")
    ax[0].legend()
    
    # Graphique de la vitesse de la comète en fonction du temps
    ax[1].plot(subdivision, approx_euler[2], label="vitesse q1 euler")
    ax[1].plot(subdivision, approx_euler[3], label="vitesse q2 euler")
    ax[1].set_title("Vitesse de la comète en fonction du temps")
    ax[1].set_xlabel("temps")
    ax[1].set_ylabel("vitesse")
    ax[1].legend()
    
    # Équation de l'énergie de la comète
    E = lambda t, q1, q2, q3, q4: (m/2) * (q3 ** 2 + q4 ** 2) - G * m ** 2 / np.sqrt((q1 ** 2 + q2 ** 2))
    ax[2].plot(subdivision, [E(None, approx_euler[0][t], approx_euler[1][t], approx_euler[2][t], approx_euler[3][t]) for t in range(len(subdivision))])
    ax[2].set_title("Énergie de la comète en fonction du temps")
    ax[2].set_xlabel("temps")
    ax[2].set_ylabel("énergie")
    ax[2].legend()
    
    filepath = './Results/trajectoires_cometeT={}'.format(T)
    
    plt.tight_layout()
    plt.savefig(filepath)
    print("Le graphique a été sauvegardé dans ", filepath)
    plt.show()
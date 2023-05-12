# On s'intérèsse au problème de N corps. C'est un problèmes de mécanique céleste
# qui consiste à calucler la trajectoires de N corps s'attirant mutuellement.

# Dans ce projet nous nous intéresserons unique à l'orbite d'une comète tournant autour
# du soleil. Le soleil est supposé fixe. Les trajectoires de la comète se trouvent sur un plan

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

def schema_euler_explicite(T, N, y_init, f):
    # Initialisation
    y = np.zeros((len(y_init), N+1))
    y[:,0] = y_init
    h = T / N
    t = np.linspace(0, T, num=N+1)
    
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
    Id = np.diag([1] * len(y_init))
    
    # Itération
    for n in range(N):
        # y[:, n+1] = y[:,n] * np.linalg.inv(Id -  h * f(t[n], y[:,n]))
        pass
        
    return t, y

def erreur(exacte, approx):
    return np.max(np.abs(exacte - approx))

if __name__ == "__main__":
    # Constantes
    m = 1 # masse du soleil
    G = 1 # Constante gravitationnelle
    T = 100
    N = np.array([100, 1000, 5000, 50000])
    y_init = np.array([1., 1., 0., 1.])
    
    # On considère l'équation Y'=f(t, Y1, Y2) où Y1, Y2 \in \R^2 et Y = (y, y') et on cherche à trouver y \in \R^2
    f = lambda t, q: np.array([q[2], q[3], 
        -m * G * (q[0] / ((q[0]**2 + q[1]**2)**(3/2))), 
        -m * G * (q[1] / ((q[0]**2 + q[1]**2)**(3/2)))])
    
    
    subdivision, approx_euler = schema_euler_explicite(T, N[-1], y_init, f)
    _, approx_semi_implicite = schema_semi_implicite(T, N[-1], y_init, f)
    
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
    
    # Équation de l'énergie de la comète; q \in \R^4
    E = lambda t, q: (m/2) * (q[2] ** 2 + q[3] ** 2) - G * m ** 2 / ((q[0] ** 2 + q[1]** 2)**(1/2))
    # ax[2].plot(subdivision, [E(None, approx_euler[0][t], approx_euler[1][t], approx_euler[2][t], approx_euler[3][t]) for t in range(len(subdivision))])
    ax[2].plot(subdivision, E(subdivision, approx_euler), label="Euler")
    ax[2].set_title("Énergie de la comète en fonction du temps")
    ax[2].set_xlabel("temps")
    ax[2].set_ylabel("énergie")
    ax[2].legend()
    
    filepath = './Results/trajectoires_cometeT={}N={}'.format(T, N[-1])
    
    plt.tight_layout()
    plt.savefig(filepath)
    print("Le graphique a été sauvegardé dans ", filepath)
    plt.show()
    
    # Calcul de l'erreur
    erreur_euler = []
    
    for n in N:
        subdivision, approx_euler = schema_euler_explicite(T, n, y_init, f)
        # plt.plot(subdivision, E(subdivision, approx_euler), label="N={}".format(n))
        # energie_exacte = np.array([-0.2 for t in subdivision])
        energie_exacte = np.array([(4 * G * m * t) ** (1/2) for t in subdivision])
        erreur_euler.append(erreur(energie_exacte, E(subdivision, approx_euler)))
    
    # print("Ordre de convergence : ", stats.linregress(np.log(T/n), np.log(erreur_euler)))
    
    # plt.plot(0,0)
    plt.plot(np.log(N), np.log(erreur_euler), label="Erreur euler")
    plt.title("Erreur")
    plt.xlabel("log(N)")
    plt.ylabel("log(erreur)")
    plt.legend()
    plt.show()
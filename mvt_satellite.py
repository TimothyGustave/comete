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
        y[:, n+1] = np.linalg.inv(Id -  h * f(t[n], y[:,n])) @ y[:,n] 
        
    return t, y

def erreur(exacte, approx):
    return np.max(np.abs(exacte - approx))

def A(q):
    # J'ai pas encore trouvé de manière élégante de faire ça
    A = np.zeros((4, 4))
    nr = (q[0]**2 + q[1]**2)**(3/2)
    
    # if nr == 0:
    #     print(q)
    
    A[0, 2] = 1
    A[1, 3] = 1
    A[2, 0] = -m * G / nr
    A[3, 1] = -m * G / nr
    
    return A

if __name__ == "__main__":
    # Constantes
    m = 1 # masse du soleil
    G = 1 # Constante gravitationnelle
    T = 100
    N = np.array([100, 1000, 5000, 50000])
    ys_init = np.array([[1., 1., 0., 1.],
                       [1., 1., 0.5, 1.],
                       [1., 1., -1., 1.]])
    
        
    # On considère l'équation Y'=f(t, Y1, Y2) où Y1, Y2 \in \R^2 et Y = (y, y') et on cherche à trouver y \in \R^2
    f = lambda t, q: np.array([q[2], q[3], 
        -m * G * (q[0] / ((q[0]**2 + q[1]**2)**(3/2))), 
        -m * G * (q[1] / ((q[0]**2 + q[1]**2)**(3/2)))])
    
    g = lambda t, q: A(q) # C'est moche il faut essayer de trouver quelques choses de mieux
    
    E = lambda t, q: (m/2) * (q[2] ** 2 + q[3] ** 2) - G * m ** 2 / ((q[0] ** 2 + q[1]** 2)**(1/2))
    
    
    # Représentations graphiques
    fig, ax = plt.subplots(1, 3)
    titre = ["Ellipse", "Parabole", "Hyperbole"]
    
    for i in range(3):
        
        y_init = ys_init[i]
        
        subdivision, approx_euler = schema_euler_explicite(T, N[-1], y_init, f)
        _, approx_semi_implicite = schema_semi_implicite(T, N[-1], y_init, g)
    
        # Graphique de la trajctoires de la comète
        # On va par la suite tracer toutes les situations possibles des coniques.
        ax[i].plot(approx_euler[0], approx_euler[1], label="euler explicite")
        ax[i].plot(approx_semi_implicite[0], approx_semi_implicite[1], label="semi_implicite")
        ax[i].set_title(titre[i])
        ax[i].set_xlabel("position sur l'axe q1")
        ax[i].set_ylabel("position sur l'axe q2")
        ax[i].legend()
        
    print(approx_semi_implicite)
        
    filepath_trajectoires = './Results/trajectoires_comete_coniquesT={}N={}'.format(T, N[-1])
    
    plt.tight_layout()
    plt.savefig(filepath_trajectoires)
    print("Le graphique a été sauvegardé dans ", filepath_trajectoires)
    plt.show()
    
    # Calcul de l'erreur
    fig, ax = plt.subplots(3)
    
    filepath_energie = "./Results/energie"
    
    erreur_euler = []
    erreur_semi_implicite = []
    
    for n in N:
        subdivision, approx_euler = schema_euler_explicite(T, n, ys_init[0], f)
        _, approx_semi_implicite = schema_semi_implicite(T, n, ys_init[0], g)
        
        print("La longueur: ", len(y_init * len(subdivision)), len(subdivision))
        
        # energie_exacte = E(subdivision, y_init * (len(subdivision) / 4)) # l'énergie est constante, et donc égale au conditions initiales
        energie_exacte = [E(t, y_init) for t in subdivision]
        
        print("La 2eme longueur: ", len(energie_exacte))
        print(energie_exacte)
        
        energie_approx_euler = E(subdivision, approx_euler)
        energie_approx_semi_implicite = E(subdivision, approx_semi_implicite)
        
        erreur_euler.append(erreur(energie_exacte, energie_approx_euler))
        erreur_semi_implicite.append(erreur(energie_exacte, energie_approx_semi_implicite))
        
        ax[1].plot(subdivision, E(subdivision, approx_euler), label="N={}".format(n))
        if n == 50000:
            ax[2].plot(subdivision, E(subdivision, approx_semi_implicite), label="N={}".format(n))
    
    # Tracer de l'énérgie en fonction du temps
    ax[0].plot(approx_euler[0], approx_euler[1], label="euler explicite")
    ax[0].plot(approx_semi_implicite[0], approx_semi_implicite[1], label="semi_implicite")
    ax[0].set_title("Trajectoires de la comète")
    ax[0].set_xlabel("q1")
    ax[0].set_ylabel("q2")
    ax[0].legend()
    ax[1].set_title("Énergie pour le schéma d'Euler Explicite")
    ax[1].set_xlabel("temps")
    ax[1].set_ylabel("énergie")
    ax[1].legend()
    ax[2].set_title("Énergie pour le schéma semi-implicite")
    ax[2].set_xlabel("temps")
    ax[2].set_ylabel("énergie")
    ax[2].legend()
    
    plt.tight_layout()
    plt.savefig(filepath_energie)
    print("La figure a été sauvé a l'emplacement : ", filepath_energie)
    plt.show()
    
    print(erreur_semi_implicite)
    
    # Tracer de l'erreur et de l'ordre de convergence
    regression_lineaire_euler = stats.linregress(np.log(N), np.log(erreur_euler))
    regression_lineaire_semi_implicite = stats.linregress(np.log(N), np.log(erreur_semi_implicite))
    ordre_convergence_euler = regression_lineaire_euler.slope
    ordre_convergence_semi_implicite = regression_lineaire_semi_implicite.slope
    # ordonne_origine_euler = regression_lineaire_euler
    print(regression_lineaire_euler)
    print("Ordre de convergence du schéma d'euler : ", ordre_convergence_euler)
    print("Ordre de convergence du schéma semi implicite : ", ordre_convergence_semi_implicite)
    
    plt.plot(np.log(N), np.log(erreur_euler), label="Erreur euler")
    plt.plot(np.log(N), np.log(erreur_semi_implicite), label="Erreur semi implicite")
    plt.title("Erreur")
    plt.xlabel("log(N)")
    plt.ylabel("log(erreur)")
    plt.legend()
    plt.savefig("./Results/erreur.png")
    print("La figure a été sauvée à l'emplacement : ", "./Results/erreur.png")
    plt.show()
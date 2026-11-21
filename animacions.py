import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import PillowWriter

from Euler_explicit import x, matriu_1, matriu_2, matriu_3, dt_1, dt_2, dt_3

#LES ANIMACIONS ES PODEN VEURE A LA CARPETA "ANIMAVIONS" DEL REPOSTORI

def crear_animacio(matriu, dt, nom_fitxer):

    fig, ax = plt.subplots()

    indices = np.linspace(0, len(matriu)-1, min(200, len(matriu)), dtype=int)

    linia, = ax.plot(x, matriu[0])

    ax.set_xlim(0, 1)
    ax.set_ylim(np.min(matriu), np.max(matriu))

    ax.set_xlabel('Posició x')
    ax.set_ylabel('Temperatura T(x,t)')
    ax.set_title('Evolució temporal de la temperatura')

    def actualitzar(frame):
        linia.set_ydata(matriu[frame])
        temps = frame * dt
        ax.set_title(f'Evolució de la temperatura: t = {temps:.5f}')
        return linia,                      # <-- coma: retorna una tupla

    animacio = animation.FuncAnimation(fig, actualitzar, frames=indices,
                                       interval=40, blit=False)   # <-- blit=False

    animacio.save(nom_fitxer, writer=PillowWriter(fps=25))

    return animacio

#ANIMACIONS EULER EXPLÍCIT
crear_animacio(matriu_1, dt_1, 'Euler_explicit_051.gif')
crear_animacio(matriu_2, dt_2, 'Euler_explicit_049.gif')
crear_animacio(matriu_3, dt_3, 'Euler_explicit_025.gif')

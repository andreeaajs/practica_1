import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import PillowWriter

from Euler_explicit import x, matriu_1, matriu_2, matriu_3, dt_1, dt_2, dt_3

# Configurem l'estil i el format dels gràfics.

plt.rcParams.update({
    "font.family": "serif", # font de la lletra
    "mathtext.fontset": "cm", # les formules es fan amb "computer modern" igual que LaTeX
    "font.size": 10, # mida de la lletra
    "axes.labelsize": 10, # mida de les etiquetes
    "legend.fontsize": 8, # mida del text de la llegenda
    "xtick.labelsize": 9, # mida dels números de les marques dels eixos
    "ytick.labelsize": 9, 
    "xtick.direction": "in", "ytick.direction": "in", # Les marques cap a "dins"
    "xtick.top": True, "ytick.right": True, # dibuxa també marques als costats dret i superior
    "xtick.minor.visible": True, "ytick.minor.visible": True, # Afegeix marques menors
    "axes.linewidth": 0.8, # Gruix de la caixa dels eixos
    "lines.linewidth": 1.2, # gruix de les corbes de plot
    "legend.frameon": False, # Treu la quadricula
    "figure.dpi": 120, # marca la resolució
    "savefig.bbox": "tight", # si es guarda la figura, retalla l'espai en blanc
})


# LES ANIMACIONS ES PODEN VEURE A LA CARPETA "ANIMACIONS" DEL REPOSTORI.

# Definim una funció per generar i guardar l'animació de la temperatura.
def crear_animacio(matriu, dt, nom_fitxer):

    fig, ax = plt.subplots()

    indices = np.linspace(0, len(matriu)-1, min(200, len(matriu)), dtype=int)   # Seleccionem fins a 200 fotogrames distribuïts uniformement al llarg de tota la simulació.

    # Comencem la gràfica amb la temperatura inicial i establim els límits dels eixos.
    linia, = ax.plot(x, matriu[0])
    ax.set_xlim(0, 1)
    ax.set_ylim(np.min(matriu), np.max(matriu))
    ax.set_xlabel('Posició x')
    ax.set_ylabel('Temperatura T(x,t)')
    ax.set_title('Evolució temporal de la temperatura')

    # Actualitzem la temperatura i el temps corresponent en cada fotograma de l'animació.
    def actualitzar(frame):
        linia.set_ydata(matriu[frame])
        temps = frame * dt
        ax.set_title(f'Evolució de la temperatura: t = {temps:.5f}')
        return linia,   # Retornem la línia gràfica actualitzada dins d'una tupla d'un sol element.
    
    # Creem l'animació a partir dels índexs de temps seleccionats.
    # interval=40 indica el temps entre fotogrames (en ms) i blit=False
    # fa que es redibuixque la figura completa en cada actualització.
    animacio = animation.FuncAnimation(fig, actualitzar, frames=indices,
                                       interval=40, blit=False)
    
    # Guardem l'animació en format GIF amb una velocitat de 25 fotogrames per segon.
    animacio.save(nom_fitxer, writer=PillowWriter(fps=25))

    plt.show()
    plt.close(fig)

    return animacio

# ANIMACIONS EULER EXPLÍCIT
crear_animacio(matriu_1, dt_1, 'Euler_explicit_051.gif')
crear_animacio(matriu_2, dt_2, 'Euler_explicit_049.gif')
crear_animacio(matriu_3, dt_3, 'Euler_explicit_025.gif')
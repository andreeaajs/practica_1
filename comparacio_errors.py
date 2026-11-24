import numpy as np
import matplotlib.pyplot as plt
from Euler_explicit import x, err_1 as e1, err_2 as e2, err_3 as e3
from Euler_implicit import error_4 as e4, error_5 as e5
from Crank_Nicolson import error_4 as e6, err_5 as e7

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


# Definim els colors de cada mètode.
c_exp, c_imp, c_cn = "#0072B2", "#D55E00", "#009E73"

# Representem els errors numèrics dels diferents mètodes.
plt.figure()

plt.plot(x,e2, "-", label='explícit $0.49 \Delta \hat{x}^2$', color = c_exp)
plt.plot(x,e3, "--", label='explícit $0.25 \Delta \hat{x}^2$', color=c_exp)
plt.plot(x,e4, "-", label='implícit $ \Delta \hat{x}^2$', color=c_imp)
plt.plot(x,e5, "--", label='implícit $0.5 \Delta \hat{x}^2$', color=c_imp)
plt.plot(x,e6, "-", label='C-N $ \Delta \hat{x}^2$', color=c_cn)
plt.plot(x,e7, "--", label='C-N $0.5 \Delta \hat{x}^2$', color=c_cn)

plt.ylim(0,1.8*10**(-5))
plt.xlabel('$\hat{x}$')
plt.ylabel('Error numèric')

plt.legend(ncol=2, loc="upper center", fontsize=9, columnspacing=1.0)
plt.tight_layout()
plt.show()

# Valor màxim dels errors en cada cas.
print('Error màxim Euler explícit amb gamma 0.51:', max(e1))
print('Error màxim Euler explícit amb gamma 0.49:', max(e2))
print('Error màxim Euler explícit amb gamma 0.25:', max(e3))
print('Error màxim Euler implícit amb gamma 1:', max(e4))
print('Error màxim Euler implícit amb gamma 0.5:', max(e5))
print('Error màxim Cranck-Nicolson amb gamma 1:', max(e6))
print('Error màxim Cranck-Nicolson amb gamma 0.5:', max(e7))
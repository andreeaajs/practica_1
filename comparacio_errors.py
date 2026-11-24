import numpy as np
import matplotlib.pyplot as plt
from Euler_explicit import x, err_1 as e1, err_2 as e2, err_3 as e3
from Euler_implicit import error_4 as e4, error_5 as e5
from Crank_Nicolson import error_4 as e6, err_5 as e7

# Configurem l'estil i el format dels gràfics.

plt.rcParams.update({
    "font.family": "serif",
    "mathtext.fontset": "cm", 
    "font.size": 10,
    "axes.labelsize": 10,
    "legend.fontsize": 8,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
    "xtick.minor.visible": True, "ytick.minor.visible": True,
    "axes.linewidth": 0.8,
    "lines.linewidth": 1.2,
    "legend.frameon": False,
    "figure.dpi": 120,
    "savefig.bbox": "tight",
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
plt.title('Comparació dels errors numèrics dels diferents mètodes')
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
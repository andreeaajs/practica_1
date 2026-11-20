import numpy as np
import _euler_explicit_ 
import matplotlib.pyplot as plt
from solucio_analitica import T

# Perque els grafics quedin millor------------------------------------------------------------

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

#---------------------------------------------------------------------------------------------


#PARÀMETRES

sigma, V_pic, L = 0.472, 40.0, 0.02      
kappa = 0.56                             
Tc = 36.5                                

P_ext = sigma*V_pic**2/(2*L**2)        
dTc   = P_ext*L**2/kappa                 

beta = (Tc + 273.15)/dTc               
beta_C = Tc/dTc                          
T_50 = (50 - Tc)/dTc                     
T_80 = (80 - Tc)/dTc 

print(beta,T_50,T_80)

N = 101
t_a = 0.025
T_a = 0.459
x = np.linspace(0, 1, N)





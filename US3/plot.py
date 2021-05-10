import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.constants import codata  
import scipy.constants as const
import numpy as np
import uncertainties as unc
import uncertainties.unumpy as unp
from uncertainties.unumpy import (nominal_values as noms, std_devs as stds)
import math


k15 = np.array([ 85  , 159 , 344 , 342  , 464  , 586 ])
k30 = np.array([-134 ,-281 ,-439 ,-659  ,-903  ,-1172])
k60 = np.array([ 232 , 537 , 903 , 1318 , 1893 , 2271])
m15 = np.array([ 49  , 61  , 98  , 146  , 195  , 0   ])
m30 = np.array([-61  ,-98  ,-146 ,-244  ,-342  , 0   ])
m60 = np.array([ 110 , 220 , 366 , 525  , 537  , 0   ])
g15 = np.array([ 0   , 49  , 61  , 85   , 122  , 134 ])
g30 = np.array([ 0   ,-61  ,-85  ,-146  ,-208  ,-269 ])
g60 = np.array([ 61  , 98  , 159 , 232  , 330  , 427 ])

#Strömungsgeschwindigkeiten berechnen
c = 1800
f = 2000000
def F(Frequenzverschiebung, Dopplerwinkel):
    return (Frequenzverschiebung * c) / (np.cos(Dopplerwinkel * (np.pi / 180)) * 2 * f)

vk15 = F(k15, 80.06)
vk30 = F(k30, 70.57)
vk60 = F(k60, 54.74)
vm15 = F(m15, 80.06)
vm30 = F(m30, 70.57)
vm60 = F(m60, 54.74)
vg15 = F(g15, 80.06)
vg30 = F(g30, 70.57)
vg60 = F(g60, 54.74)

Geschw = np.array([vk15, vk30, vk60, vm15, vm30, vm60, vg15, vg30, vg60])

np.savetxt('data/Geschwindigkeiten.dat', (Geschw), header='v')

#plot
def A(x):
    return (x * 2 * f) / c

a1 = A(vk15)
a2 = A(vm15)
a3 = A(vg15)

plt.errorbar(vk15, a1, fmt='r.', label=r'7mm / 15$^{\circ}$')
plt.errorbar(vm15, a2, fmt='b.', label=r'10mm / 15$^{\circ}$')
plt.errorbar(vg15, a3, fmt='g.', label=r'16mm / 15$^{\circ}$')
plt.xlabel(r'Strömungsgeschwindigkeit [$\frac{m}{s}$]')
plt.ylabel(r'Δ$v / cos(\alpha)')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('plot1.pdf')
plt.close('all')

#Strömungsprofil

f70 = np.array([42, 305, 404, 317, 256, 122]) #Frequenzverschiebung
f45 = np.array([ 0, 134, 171, 159, 134, 110])
i70 = np.array([19, 61 , 60 , 69 , 93 , 49 ])/100 #Intensität
i45 = np.array([ 0, 53 , 69 , 80 , 74 , 54 ])/100
Tiefe = np.array([13,14,15,16,17,18])


v70 = F(f70, 80.06) #Geschwindigkeit
v45 = F(f45, 80.06)

plt.errorbar(Tiefe, v70, fmt='r.', label=r'Geschwindigkeit')
plt.errorbar(Tiefe, i70, fmt='b.', label=r'Intensität')
plt.xlabel(r'Messtiefe [$\mu s$]')
plt.ylabel(r'I [$\frac{10MV^2}{s}$] / $v$ [$\frac{m}{s}]')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('plot2.pdf')
plt.close('all')


plt.errorbar(Tiefe, v45, fmt='r.', label=r'Geschwindigkeit')
plt.errorbar(Tiefe, i45, fmt='b.', label=r'Intensität')
plt.xlabel(r'Messtiefe [$\mu s$]')
plt.ylabel(r'Strömungsgeschwindigkeit[$\frac{m}{s}]')
plt.ylabel(r'I [$\frac{10MV^2}{s}$] / $v$ [$\frac{m}{s}]')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('plot3.pdf')
plt.close('all')

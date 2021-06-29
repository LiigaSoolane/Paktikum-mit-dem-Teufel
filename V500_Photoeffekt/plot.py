import matplotlib.pyplot as plt
import numpy as np
from tabulate import tabulate

def pl(x,y,a,c, alpha, omega):
    y=np.sqrt(y)
    params, cov = np.polyfit(x[alpha:omega], y[alpha:omega], deg=1, cov=True)
    errors = np.sqrt(np.diag(cov))
    #for name, value, error in zip('ab', params, errors):
    #    print(f'{name} = {value:.3f} ± {error:.3f}')
    x0=np.linspace(x[alpha], x[omega], 100)
    y0=params[1]+x0*params[0]
    c+1
    plt.figure(c)
    plt.plot(x, y, 'r*', label='Messwerte')
    plt.plot(x0, y0, label='Messwerte')
    plt.xlabel(r'U [V]')
    plt.ylabel(r'$\sqrt{I}$')
    plt.legend(loc='best')
    plt.savefig(a)
    plt.close()
    ug=params[1]/params[0]
    return c, ug

redU, redI = np.genfromtxt('drot.txt',unpack=True)
greenU, greenI = np.genfromtxt('dgrün.txt',unpack=True)
blueU, blueI = np.genfromtxt('dblau.txt',unpack=True)
yelU, yelI = np.genfromtxt('dgelb.txt',unpack=True)
yellU, yellI = np.genfromtxt('dgelb1.txt',unpack=True)

c=-1
c, redug =pl(redU, redI, 'build/plot.pdf', c, 6, 10)
c, greenug =pl(greenU, greenI, 'build/fuckyou.pdf', c, 8, 12)
c, blueug =pl(blueU, blueI, 'build/hurensohn.pdf', c, 8, 13)
c, yelug =pl(yelU, yelI, 'build/verdammtescheiße.pdf', c, 18, 28)
c, yellug =pl(yellU, yellI, 'build/fresse.pdf', c, 8, 11)

lamda = np.array([614, 577.579, 577.579, 546, 434])
lamda = 299792458*10**9/(lamda)
ug = np.array([redug, yelug, yellug, greenug, blueug])
#table ={'1': lamda, '3': ug}
#print(tabulate (table, tablefmt="latex"))
params, cov = np.polyfit(lamda, ug, deg=1, cov=True)
errors = np.sqrt(np.diag(cov))

for name, value, error in zip('ab', params, errors):
    print(f'{name} = {value} ± {error}')


x0=np.linspace(lamda[0], lamda[4], 100)
y0=params[1]+x0*params[0]
a=params[0]

c+1
plt.figure(c)
plt.plot(lamda, ug, 'r*', label=r'für $U_G$ bestimmte Werte')
plt.plot(x0, y0, label='lineare Regression')
plt.xlabel('f [Hz]')
plt.ylabel(r'$U_G$ [V]')
plt.legend(loc='best')
plt.savefig('build/torte.pdf')

theo=4.136*10**(-15)
print((theo-a)/theo)


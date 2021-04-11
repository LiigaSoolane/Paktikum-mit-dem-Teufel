import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from uncertainties import ufloat

t, T1, p1, T2, p2, N = np.genfromtxt('Daten.dat', unpack=True)
T1 = T1+273.15
T2 = T2+273.15
for i in range(35):
    UT1=ufloat(T1[i], 0.1)
    UT2=ufloat(T2[i], 0.1)
p1 = p1 + 1
p2 = p2 + 1
t = 60*t
x = np.linspace(t[0], t[35], 80)
q = [5*60,15*60,25*60,35*60]
m1 = 4
cw = 4184
mc = 750
kappa=1.14
rho=5.51*10**(-3)
T=273.15




#a)
plt.plot(t, T1, 'r*', label='T1')
plt.plot(t, T2, 'g*', label='T2')
plt.xlabel(r'Zeit [sec]')
plt.ylabel(r'Temperatur [K]')
plt.legend(loc='best')
plt.savefig('build/plot1.pdf')


#b)
def func1(x, a, b, c):
    return a*x**2+b*x+c

params, cov = curve_fit(func1, t, T1)
a, b, c = params[0], params[1], params[2]
yfit1 = a*x**2+b*x+c
errors = np.sqrt(np.diag(cov))
print('a,b,c = ', params)
print('errors = ', errors)
a, b, c = ufloat(params[0], errors[0]), ufloat(params[1], errors[1]), ufloat(params[2], errors[2])


params, cov = curve_fit(func1, t, T2)
errors1 = np.sqrt(np.diag(cov))
a2, b2, c2 = params[0], params[1], params[2]
yfit2 = a2*x**2+b2*x+c2
print('a,b,c = ', params)
print('errors = ', errors1)
a2, b2, c2 = ufloat(params[0], errors1[0]), ufloat(params[1], errors1[1]), ufloat(params[2], errors1[2])

plt.plot(x, yfit1, 'r-', label= 'Ausgleichsgerade T1')
plt.plot(x, yfit2, 'g-', label= 'Ausgleichsgerade T2')
plt.savefig('build/plot2.pdf')

#c)
def func2(x,a,b):
    return 2*a*x+b
for i in q:
    v=func2(i,a,b)
    vi=func2(i,a2,b2)
    print('dT1/dt = ', v)
    print('dT2/dt = ', vi)

#d)
k=5
for i in q:
    v=(m1*cw+mc)*func2(i,a,b)/120
    va = T1[k]/(T1[k]-T2[k])
    k = k+10
    print('v = ', v)
    print('vid = ', va)

#e)
v=[0,0,0,0]
k=0
for i in q:
    v[k]=(m1*cw+mc)*func2(i,a2,b2)
    k=k+1
print('dQ2/dt = ', v)

x_plot=1/T1
y_plot=np.log(p1)
xfit = np.linspace(x_plot[0],x_plot[35],80)
params, cov = np.polyfit(x_plot, y_plot, deg=1, cov=True)
print('a,b = ',params)
errors2 = np.sqrt(np.diag(cov))
print('aerr,berr=',errors2)
L = ufloat(-params[0]*8.31446261815324,errors2[0]*8.31446261815324)
error=errors2[0]*8.31446261815324
print('Lerr=',error)
print('L = ', L)

plt.figure(0)
plt.plot(x_plot, y_plot, 'r*', label='Messdaten des Reservoir 1')
plt.plot(xfit, params[0]*xfit+params[1], 'g-', label='fit')
plt.ylabel(r'ln($\dfrac{p}{p_0}$)')
plt.xlabel(r'1/Temperatur [1/K]')
plt.legend(loc='best')
plt.savefig('build/plot3.pdf')

mt=v
k=0
for i in v:
    mt[k]=i/L
    k=k+1
print('m/t = ',mt)
k=0
for i in v:
    mt[k]=mt[k]*120.9
    k=k+1
print('m/t = ',mt)

#f)
k=5
h=0
nmech=([0,0,0,0])
for i in q:
    ro=1/(T2[k]/(rho*T*p2[k]))
    nmech[h]=(1/(kappa-1))*((p1[k]*((p2[k]/p1[k])**(1/kappa)))-p2[k])*mt[h]/(ro * 1000)
    h=h+1
    k=k+10
    print('rho = ',ro)
print('Nmech = ', nmech)
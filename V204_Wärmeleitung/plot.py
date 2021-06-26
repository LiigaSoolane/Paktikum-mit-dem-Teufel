import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks
import uncertainties.unumpy as unp
from uncertainties import ufloat

def amplitudes(temp):
    maxs, _ = find_peaks(temp, distance = 15)
    mins, _ = find_peaks(-temp, distance = 30)
    #print(mins_ind)

    #print(temp[maxs_ind], temp[mins_ind])
    #amps = np.empty(len(maxs_ind))
    #dt = np.empty(len(maxs_ind))
#
    #dt[0] = t[mins_ind[0]]
    #amps[0] = temp[maxs_ind[0]] - temp[0] 
#
    #j = 1
    #while j < len(maxs_ind):
    #    amps[j] = temp[maxs_ind[j]] - temp[mins_ind[j-1]]
    #    dt[j] = t[mins_ind[j]] - t[mins_ind[j-1]]
#
    #    j += 1

    print(maxs)
    print(mins)
    amps = 0
    
    return amps
    


###### Statische Messung

n_stat, stat_1, stat_2, stat_3, stat_4, stat_5, stat_6, stat_7, stat_8 = np.genfromtxt('statisch.txt', unpack=True)

t_stat = n_stat/5

stat_1 += 273.15
stat_2 += 273.15
stat_3 += 273.15
stat_4 += 273.15
stat_5 += 273.15
stat_6 += 273.15
stat_7 += 273.15
stat_8 += 273.15

### Plots
plt.figure()
plt.plot(t_stat, stat_1, label = 'Messing, breit')
plt.plot(t_stat, stat_4, label = 'Messing, schmal')

plt.xlabel('Zeit [s]')
plt.ylabel('Temperatur [°C]')
plt.legend()

plt.savefig('verlauf_mess.pdf')
plt.clf()

plt.figure()
plt.plot(t_stat, stat_5, label = 'Aluminium')
plt.plot(t_stat, stat_8, label = 'Edelstahl')


plt.xlabel('Zeit [s]')
plt.ylabel('Temperatur [°C]')
plt.legend()

plt.savefig('verlauf_alu_edel.pdf')
plt.clf()

plt.figure()
plt.plot(t_stat, stat_7 - stat_8, label = 'Temperaturdifferenz Edelstahl')
plt.plot(t_stat, stat_2 - stat_1, label = 'Temperaturdifferenz Messing, breit')

plt.xlabel('Zeit [s]')
plt.ylabel('Temperatur [°C]')

plt.legend()

plt.savefig('differenz_stat.pdf')
plt.clf()


### Wärmestrom
names = np.array(['Messing, breit', 'Messing, schmal', 'Aluminium', 'Edelstahl'])
temps = np.array([stat_1, stat_4, stat_5, stat_8])
diffs = np.array([stat_1 - stat_2, stat_3 - stat_4, stat_5 - stat_6, stat_7 - stat_8])

dx = 3 * 1e-2

A = (1.2 * 1e-2) * (0.4 * 1e-2)
A_schmal = (0.7 * 1e-2) * (0.4 * 1e-2)

rho_mess = 8520
rho_alu = 2800
rho_edel = 8000

c_mess = 385
c_alu = 830
c_edel = 400


### Wärmeleitfähigkeit

kappa_mess = 120
kappa_alu = 235
kappa_edel = 21

n = np.array([5, 300 * 5, 600 * 5, 900 * 5, 5403])
cs = np.array([c_mess, c_mess, c_alu, c_edel])
kappas = np.array([kappa_mess, kappa_mess, kappa_alu, kappa_edel])
As = np.array([A, A_schmal, A, A])


### despues de setecientos segundos

k = 700 * 5
print(f'Temp. nach 700: 1:{stat_1[k]:.3f}, schm:{stat_4[k]:.3f}, Alu:{stat_5[k]:.3f}, Edel:{stat_8[k]:.3f}')


### stuff en cinco tiempos differentes

j = 0

for name, c, kappa, A, temp, diff in zip(names, cs, kappas, As, temps, diffs):
    print(f'{name}->')
    j = 0
    while j < 5:
        Qdt = - kappa * A * diff[n[j]]/dx
        print(f'{name}_Q: {Qdt:.3f}, diff: {diff[n[j]]:.3f}')
        j += 1


###### dynamische Messung, 80s
n_dyn1, dyn1_1, dyn1_2, dyn1_3, dyn1_4, dyn1_5, dyn1_6, dyn1_7, dyn1_8 = np.genfromtxt('dynamisch1.txt', unpack=True)
t_dyn1 = n_dyn1 * 2

dyn1_1 += 273.15
dyn1_2 += 273.15
dyn1_3 += 273.15
dyn1_4 += 273.15
dyn1_5 += 273.15
dyn1_6 += 273.15
dyn1_7 += 273.15
dyn1_8 += 273.15

plt.figure()
plt.plot(t_dyn1, dyn1_1, label = 'außen')
plt.plot(t_dyn1, dyn1_2, label = 'innen')

plt.xlabel('Zeit [s]')
plt.ylabel('Temperatur [°C]')
plt.legend()

#plt.grid()
plt.savefig('dyn_1_mess.pdf')
plt.clf()

plt.figure()
plt.plot(t_dyn1, dyn1_4, label = 'außen')
plt.plot(t_dyn1, dyn1_5, label = 'innen')

plt.xlabel('Zeit [s]')
plt.ylabel('Temperatur [°C]')
plt.legend()

#plt.grid()
plt.savefig('dyn_1_alu.pdf')
plt.clf()


### get amplitudes and phase shift

A1 = amplitudes(dyn1_1)
A2 = amplitudes(dyn1_2)

A6 = amplitudes(dyn1_6)
A5 = amplitudes(dyn1_5)

A1_max = np.array([ 149,  297,  447,  567,  684,  841,  999, 1158, 1317, 1475, 1637])
A1_min = np.array([ 0, 174,  341,  504,  587,  746,  907, 1068, 1228, 1389, 1548])
A1 = dyn1_1[A1_max] - dyn1_1[A1_min]

A2_max = np.array([ 93,  252,  412,  535,  651,  810,  970, 1130, 1290, 1450, 1610])
A2_min = np.array([ 0, 166,  327,  487,  575,  728,  888, 1049, 1208, 1368, 1528])
A2 = dyn1_2[A2_max] - dyn1_2[A2_min]

dt_12_ = (t_dyn1[A1_max] -  t_dyn1[A2_max])/4
dt_12 = ufloat(np.mean(dt_12_), np.std(dt_12_, ddof=1))

A1_m = ufloat(np.mean(A1), np.std(A1, ddof = 1))
A2_m = ufloat(np.mean(A2), np.std(A2, ddof = 1))

print(f'Phasendifferenz Messing {dt_12}, Amp 1 {A1_m}, 2 {A2_m}')

A6_max = np.array([  89,  248,  409,  529,  647,  807,  967, 1127, 1286, 1446, 1606])
A6_min = np.array([ 0, 162,  324,  485,  572,  725,  884, 1045, 1205, 1365, 1524])
A6 = dyn1_6[A6_max] - dyn1_6[A6_min]

A5_max = np.array([ 117,  269, 428,  547,  665,  823,  983, 1144, 1303, 1462, 1622])
A5_min = np.array([   1,  171,  334,  497,  582,  738,  898, 1058, 1218, 1378, 1538])

A5 = (dyn1_5[A5_max] - dyn1_5[A5_min])
print('a5', A5)

A6_m = ufloat(np.mean(A6), np.std(A6, ddof=1))
A5_m = ufloat(np.mean(A5), np.std(A5, ddof=1))

print('a6', A6)

dt_56_ = (t_dyn1[A5_max] -  t_dyn1[A6_max])/4
dt_56 = ufloat(np.mean(dt_56_), np.std(dt_56_, ddof=1))
print(f'Phasendifferenz Aluminium {dt_56}, Amp 5 {A5_m}, 6 {A6_m}')


kappa_mess_ex = (rho_mess * c_mess * dx**2)/(2 * dt_12 * unp.log(A2_m/A1_m))
kappa_alu_ex = (rho_alu * c_alu * dx**2)/(2 * dt_56 * unp.log(A6_m/A5_m))

abw_me = 100 * (kappa_mess_ex - kappa_mess)/kappa_mess
abw_al = 100 * (kappa_alu_ex - kappa_alu)/kappa_alu

print(f'kappa Messing {kappa_mess_ex}, Abweichung {abw_me}')
print(f'kappa Aluminium {kappa_alu_ex}, Abweichung {abw_al}')


###### dynamische Messung, 200s

n_dyn2, dyn2_1, dyn2_2, dyn2_3, dyn2_4, dyn2_5, dyn2_6, dyn2_7, dyn2_8 = np.genfromtxt('dynamisch2.txt', unpack = True)

t_dyn2 = n_dyn2 * 2

dyn2_1 += 273.15
dyn2_2 += 273.15
dyn2_3 += 273.15
dyn2_4 += 273.15
dyn2_5 += 273.15
dyn2_6 += 273.15
dyn2_7 += 273.15
dyn2_8 += 273.15

plt.figure()
plt.plot(t_dyn2, dyn2_7, label = 'außen')
plt.plot(t_dyn2, dyn2_8, label = 'innen')

plt.xlabel('Zeit [s]')
plt.ylabel('Temperatur [°C]')
plt.legend()
#plt.grid()
plt.savefig('dyn_2.pdf')
plt.clf()

plt.figure()
plt.plot(t_dyn2, dyn2_8, label = 'innen')
plt.xlabel('Zeit [s]')
plt.ylabel('Temperatur [°C]')
plt.legend()
#plt.grid()
plt.savefig('dyn_2_close.pdf')
plt.clf()


### get amplitudes and phase shift

A7 = amplitudes(dyn2_7)

A7_max = np.array([ 133,  281,  439,  557,  676,  835,  995, 1154, 1314, 1474, 1633])
A7_min = np.array([ 0, 179,  345,  506,  587,  748,  907, 1068, 1227, 1388, 1547])
A7 = dyn2_7[A7_max] - dyn2_7[A7_min]

A7_m = ufloat(np.mean(A7), np.std(A7, ddof=1))

#A8 = amplitudes(dyn2_8)
#A8_max = np.array([  19,   57, 1418])
#A8_min = np.array([  20,   53, 1414])
#A8 = dyn2_8[A8_max] - dyn2_8[A8_min]

A8_1 = np.array([29.5, 34, 37.1, 38.3, 40, 43.9])
A8_2 = np.array([31.7, 37.7, 38, 39.1, 41,  44.9])

A8 = A8_2 - A8_1

A8_m = ufloat(np.mean(A8), np.std(A8, ddof=1))
print('a8', A8_m, '; a7', A7_m)

dt_78_ = np.array([80, 100, 77, 89, 91])

dt_78 = ufloat(np.mean(dt_78_), np.std(dt_78_))

kappa_edel_ex = (rho_edel * c_edel * dx**2)/(2 * 80 * unp.log(A7_m/A8_m))

abw_ed = 100 * (kappa_edel_ex - kappa_edel)/kappa_edel
print(f'delta t {dt_78}')
print(f'kappa Edelstahl {kappa_edel_ex}, Abweichung {abw_ed}')


### Wellenlänge und Frequenz

lam_1 = unp.sqrt((4 * np.pi * kappa_mess_ex * 80)/(rho_mess * c_mess))
lam_5 = unp.sqrt((4 * np.pi * kappa_alu_ex * 80)/(rho_alu * c_alu))
lam_7 = unp.sqrt((4 * np.pi * kappa_edel_ex * 200)/(rho_edel * c_edel))

lam_1_t = unp.sqrt((4 * np.pi * kappa_mess * 80)/(rho_mess * c_mess))
lam_5_t = unp.sqrt((4 * np.pi * kappa_alu* 80)/(rho_alu * c_alu))
lam_7_t = unp.sqrt((4 * np.pi * kappa_edel * 200)/(rho_edel * c_edel))

abw_1 = 100 * (lam_1 - lam_1_t)/lam_1_t
abw_2 = 100 * (lam_5 - lam_5_t)/lam_5_t
abw_3 = 100 * (lam_7 - lam_7_t)/lam_7_t

print(f'lambda: Messing {lam_1}, Alu {lam_5}, Edelstahl {lam_7}')
print(f'lambda: Messing {lam_1_t}, Alu {lam_5_t}, Edelstahl {lam_7_t}')
print(f'Abweichung: Messing {abw_1}, Alu {abw_2}, Edelstahl {abw_3}')
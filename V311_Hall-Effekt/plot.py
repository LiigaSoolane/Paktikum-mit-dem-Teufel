import matplotlib.pyplot as plt
import numpy as np

# konstantes Uh
B1, I1, B2, I2 = np.genfromtxt('hallwerte.txt', unpack=True)

#fitte B1 und B2
params, cov = np.polyfit(I1, B1, deg=1, cov=True)
errors = np.sqrt(np.diag(cov))

x = np.linspace(0, 5, num = 1000)


plt.plot(I1, B1, 'rx', label="Messwerte")
plt.plot(x, params[0] * x + params[1], label="Lineare Regression")
plt.legend(loc="best")
plt.xlabel('Stromstärke I [A]')
plt.ylabel('Magnetische Flussdichte B [mT]')
plt.savefig('unswitched.pdf')

poroms, cav = np.polyfit(I2, B2, deg=1, cov=True)
errars = np.sqrt(np.diag(cav))

plt.plot(I2, B2, 'rx', label="Messwerte")
plt.plot(x, poroms[0] * x + poroms[1], label= "Lineare Regression")
plt.legend(loc="best")
plt.xlabel('Stromstärke I [A]')
plt.ylabel('Magnetische Flussdichte B [mT]')
plt.savefig('switched.pdf')
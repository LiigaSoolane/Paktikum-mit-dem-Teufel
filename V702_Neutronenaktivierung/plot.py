import matplotlib.pyplot as plt
import numpy as np

tv, Nv = np.genfromtxt('Vanadium.dat', unpack=True)
tr, Nr = np.genfromtxt('Rhodium.dat', unpack=True)


Nu = np.array(129, 143, 144, 136, 126, 158)
tu = 300


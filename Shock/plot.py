import matplotlib.pyplot as plt
import numpy as np
import sys

t = 0

for arg in sys.argv:
    t = arg

datos = np.loadtxt("./estado_"+str(t)+".dat")

x = datos[:,0]
v = datos[:,1]
p = datos[:,2]
rho = datos[:,3]

figura = plt.figure(figsize=(12.0, 17.0))

ax = figura.add_subplot(3, 1, 1)
ax.plot(x, v)
ax.set_title("$\mathrm{Velocidad}$")
ax.set_xlabel("$x (m)$")
ax.set_ylabel("$Vel (m/s)$")

ax = figura.add_subplot(3, 1, 2)
ax.plot(x, p)
ax.set_title("$\mathrm{Presion}$")
ax.set_xlabel("$x (m)$")
ax.set_ylabel("$Presion (kN/m^2)$")

ax = figura.add_subplot(3, 1, 3)
ax.plot(x, rho)
ax.set_title("$\mathrm{Densidad}$")
ax.set_xlabel("$x (m)$")
ax.set_ylabel("$Densidad (kg/m^3)$")

plt.savefig("estado_"+t+".pdf")
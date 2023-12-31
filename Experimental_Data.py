import numpy as np
import matplotlib.pyplot as plt

Rm = ['[0.01]', '[0.001]', '[0.1]']
K = ["{707.10678119, -873.4060921, 344.75201939, -145.05769211}", "{2236.0679775,-2088.64485917,963.92838115,-348.10668719}", "{223.607,-440.240,134.153,-72.80}"]
data = []
L = [0]
T_s = 0.00667
cut = [900,0,600] # even number

for i in range(len(Rm)):
    with open(f"directory/data ({i+1}).txt", 'r') as file:
        for line in file:
            parts = line.split(',')
            for part in parts:
                try:
                    number = float(part)
                    data.append(number)
                except ValueError:
                    continue
        L.append(len(data))


for i in range(len(Rm)):
    plt.figure(figsize=(12, 6))

    t = np.array([i for i in range(len(data[L[i]+cut[i]:L[i]+cut[i]+18000:2]))])*T_s
    pos = np.array(data[L[i]+cut[i]:L[i]+cut[i]+18000:2])-(0.89/2)
    ang = data[L[i]+1+cut[i]:L[i]+cut[i]+18000:2]

    plt.subplot(1, 2, 1)
    plt.plot(t, pos)
    plt.title("Cart Position")
    plt.ylabel("$x_c [m]$")
    plt.xlabel("t [s]")
    plt.grid()

    plt.subplot(1, 2, 2)
    plt.plot(t, ang)
    plt.title("Angle")
    plt.ylabel(r"$\theta [rad]$")
    plt.xlabel("t [s]")
    plt.grid()
    plt.show()

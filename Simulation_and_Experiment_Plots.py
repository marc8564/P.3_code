import numpy as np
import matplotlib.pyplot as plt
import control as ctrl
from matplotlib import gridspec

a_c = 3.2 # friction coeficient
a_p = 4.6e-3 # friction coeficient
m_c = 6.28 #mass of cart
m_p = 0.250 #mass of pendulum
g = 9.82 #gravitational acceleration
l = 0.3325 #length of pendulum
K_m = 93.4/1000 #motor constant
R_p = 0.028 # pulley radius

# System matrices
A_0 = np.array([[0.0, 1.0, 0.0, 0.0],
                [0.0, -a_c/m_c, m_p*g/m_c, -a_p/(l*m_c)],
                [0.0, 0.0, 0.0, 1.0],
                [0.0, -a_c/(l*m_c), g/l + m_p*g/(l*m_c), -a_p/(m_p*l**2) - a_p/(l**2*m_c)]])

B = np.array([[0.0], [1/m_c], [0.0], [1/(l*m_c)]])

#_______________________________Simulation________________________#

def ode1(S,K):
    u = -np.dot(K, S) 

    z1 = S[0]
    z2 = S[1]
    z3 = S[2]
    z4 = S[3]

    dz1dt = z2
    dz2dt = (-a_c*z2 + u - m_p * l * np.sin(z3) * z4**2 + m_p * g * np.cos(z3) * np.sin(z3) - a_p * z4 / l) * np.cos(z3) / (m_c + m_p - m_p * np.cos(z3)**2)
    dz3dt = z4
    dz4dt = (g/l) * np.sin(z3) + (-a_p * z4) / (m_p * l**2) + (np.cos(z3) * (-a_c * z2 + u - m_p * l * np.sin(z3) * z4**2 + m_p * g * np.cos(z3) * np.sin(z3) + (-a_p * z4) / l) / (l * (m_c + m_p - m_p * np.cos(z3)**2)))

    return np.array([dz1dt, dz2dt, dz3dt, dz4dt])

#RK4
def sim(IC, time, K):
    z_history = np.array([])
    I = np.array([])
    u = np.array([])
    z = IC

    t_step = 0.00667

    for i in range(len(time)):

        k1 = t_step * ode1(z, K)
        k2 = t_step * ode1(z + k1 / 2, K)
        k3 = t_step * ode1(z + k2 / 2, K)
        k4 = t_step * ode1(z + k3, K)

        I = np.append(I, np.dot(K,z)/((K_m)/R_p))
        u = np.append(u, np.dot(K,z))

        z = z + (k1 + 2 * k2 + 2 * k3 + k4) / 6
        z_history = np.append(z_history, z)

    A = z_history[::4]
    B = z_history[1::4]
    C = z_history[2::4]
    D = z_history[3::4]
    return [A,B,C,D,u,I]

t = np.linspace(0, 0.00667*1000, 1000)
S1 = np.array([0.0, 0.0, 0.1, 0.0])
S2 = np.array([0.0, 0.0, 0.5, 0.0])
S3 = np.array([0.0, 0.0, 1, 0.0])


Q = np.array([[5000, 0, 0, 0],
              [0, 1, 0, 0],
              [0, 0, 1, 0],
              [0, 0, 0, 1]])

R = np.array([0.1])

K, S_1, E = ctrl.lqr(A_0, B, Q, R)

KK = [np.array([-707.10678119, -344.75201939 ,873.4060921, 145.05769211]),np.array([-2236.0679775, -963.92838115, 2088.64485917, 348.10668719]),np.array([-223.607,-134.153,440.240,72.80])]

#________________________PLOTS__________________________________________#


Rm = ['[0.01]', '[0.001]', '[0.1]']
K_ = ["{707.10678119, -873.4060921, 344.75201939, -145.05769211}", "{2236.0679775,-2088.64485917,963.92838115,-348.10668719}", "{223.607,-440.240,134.153,-72.80}"]
data = []
L = [0]
T_s = 0.00667
cut = [900,0,600] #even number

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
    t = np.array([i for i in range(len(data[L[i]+cut[i]:L[i+1]:2]))])*T_s
    if len(t)>9000:
        t = t[:9000:]


    plt.subplot(1, 2, 1)
    plt.axhline(xmax=30, xmin=0, y=0.445, c="black")
    plt.axhline(xmax=30, xmin=0, y=-0.445, c="black")
    plt.plot(t, sim(S1, t, KK[i])[0], label="Simulation")
    plt.plot(np.array([i for i in range(len(data[L[i]+cut[i]:L[i+1]:2]))])*T_s, np.array(data[L[i]+cut[i]:L[i+1]:2])-(0.89/2), label="Experiment")

    plt.xlabel("t [s]")
    plt.ylabel("$x_c$ [m]")
    plt.title("Cart Position")
    plt.legend(loc = "upper right")
    plt.grid()

    plt.subplot(1, 2, 2)
    t = np.array([i for i in range(len(data[L[i]+1+cut[i]:L[i+1]:2]))])*T_s
    plt.plot(t, sim(S1, t, KK[i])[2], label="Simulation")
    plt.plot(np.array([i for i in range(len(data[L[i]+1+cut[i]:L[i+1]:2]))])*T_s, data[L[i]+1+cut[i]:L[i+1]:2], label = "Experiment")

    plt.xlabel("t [s]")
    plt.ylabel(r"$\theta$ [rad]")
    plt.title("Angle")
    plt.legend(loc="upper right")
    plt.grid()

    plt.suptitle("", fontsize=16)
    plt.tight_layout()
    plt.show()
    #plt.savefig(f"directory_/{i+1}", dpi=500)

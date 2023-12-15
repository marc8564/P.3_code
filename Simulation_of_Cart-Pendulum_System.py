import numpy as np
import matplotlib.pyplot as plt
import control as ctrl
from matplotlib import gridspec

a_c = 3.2
a_p = 4.6e-3
m_c = 6.28
m_p = 0.250
g = 9.82
l = 0.3325
K_m = 93.4/1000
R_m = 0.028

# System matrices
A_0 = np.array([[0.0, 1.0, 0.0, 0.0],
                [0.0, -a_c/m_c, m_p*g/m_c, -a_p/(l*m_c)],
                [0.0, 0.0, 0.0, 1.0],
                [0.0, -a_c/(l*m_c), g/l + m_p*g/(l*m_c), -a_p/(m_p*l**2) - a_p/(l**2*m_c)]])

B = np.array([[0.0], [1/m_c], [0.0], [1/(l*m_c)]])

#_______________________________Simulering________________________#

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

    return np.array([dz1dt, dz2dt[0], dz3dt, dz4dt[0]])

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

        I = np.append(I, np.dot(K,z)/((K_m)/R_m))
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

#print(f"K-matrix: {K}")
#print(f"Eigenvalues: {np.linalg.eigvals(A_0-np.dot(B,K))}")

#________________________PLOTS__________________________________________#

plt.figure(figsize=(12, 6))

gs = gridspec.GridSpec(3,3)

fig = plt.figure(1)
gridspec.GridSpec(6,4)

plt.subplot2grid((6,4), (0,0), colspan=2, rowspan=3)
plt.locator_params(axis='x', nbins=5)
plt.locator_params(axis='y', nbins=5)

plt.axhline(xmax=30, xmin=0, y=0.445, c="black")
plt.axhline(xmax=30, xmin=0, y=-0.445, c="black")
plt.plot(t, sim(S1, t, K)[0], label=r"$\theta_0 = 0.1$")
plt.plot(t, sim(S2, t, K)[0], label=r"$\theta_0 = 0.5$")
plt.plot(t, sim(S3, t, K)[0], label=r"$\theta_0 = 1.0$")

plt.xlabel("t [s]")
plt.ylabel("$x_c$ [m]")
plt.title("Cart Position")
plt.legend(loc="upper right")
plt.grid()

plt.subplot2grid((6,4), (3,0), colspan=2, rowspan=3)
plt.locator_params(axis='x', nbins=5)
plt.locator_params(axis='y', nbins=5)

plt.axhline(xmax=30, xmin=0, y=78.9079229, c="black")
plt.axhline(xmax=30, xmin=0, y=-78.9079229, c="black")
plt.plot(t, sim(S1, t, K)[5], label=r"$\theta_0 = 0.1$")
plt.plot(t, sim(S2, t, K)[5], label=r"$\theta_0 = 0.5$")
plt.plot(t, sim(S3, t, K)[5], label=r"$\theta_0 = 1.0$")

plt.xlabel('t [s]')
plt.ylabel("i [A]")
plt.title("Current")
plt.legend(loc="upper right")
plt.grid()

plt.subplot2grid((6,4), (0,2), colspan=2, rowspan=2)
plt.locator_params(axis='x', nbins=5)
plt.locator_params(axis='y', nbins=5)

plt.plot(t, sim(S1, t, K)[1], label=r"$\theta_0 = 0.1$")
plt.plot(t, sim(S2, t, K)[1], label=r"$\theta_0 = 0.5$")
plt.plot(t, sim(S3, t, K)[1], label=r"$\theta_0 = 1.0$")

plt.xlabel("t [s]")
plt.ylabel("$\dot{x}_c$ [m/s]")
plt.title("Cart Velocity")
plt.legend(loc = "upper right")
plt.grid()

plt.subplot2grid((6,4), (2,2), colspan=2, rowspan=2)
plt.locator_params(axis='x', nbins=5)
plt.locator_params(axis='y', nbins=5)

plt.plot(t, sim(S1, t, K)[2], label=r"$\theta_0 = 0.1$")
plt.plot(t, sim(S2, t, K)[2], label=r"$\theta_0 = 0.5$")
plt.plot(t, sim(S3, t, K)[2], label=r"$\theta_0 = 1.0$")

plt.xlabel("t [s]")
plt.ylabel(r"$\theta$ [rad]")
plt.title("Angle")
plt.legend(loc="upper right")
plt.grid()

plt.subplot2grid((6,4), (4,2), colspan=2, rowspan=2)
plt.locator_params(axis='x', nbins=5)
plt.locator_params(axis='y', nbins=5)

plt.plot(t, sim(S1, t, K)[3], label=r"$\theta_0 = 0.1$")
plt.plot(t, sim(S2, t, K)[3], label=r"$\theta_0 = 0.5$")
plt.plot(t, sim(S3, t, K)[3], label=r"$\theta_0 = 1.0$")

plt.xlabel("t [s]")
plt.ylabel(r" $\dot{\theta}$ [rad/s]")
plt.title("Angular Velocity")
plt.legend(loc="upper right")
plt.grid()


plt.suptitle("", fontsize=16)
plt.tight_layout()
plt.show()
#plt.savefig(f"directory", dpi=500)


C = np.eye(4)
def observability_matrix(A, C):
    n = A.shape[0]
    O = C

    for i in range(n-1):
        O = np.vstack((O, np.dot(C, np.linalg.matrix_power(A,n+1))))

    return O

print("Observability Matrix:")
print(observability_matrix(A_0, C))

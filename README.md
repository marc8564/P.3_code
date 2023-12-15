# P.3_code
The python code features 3 scripts. The first, titled "Simulation_of_Cart-Pendulum_System.py", features a simulation of a cart pendulum system, using RK4. A controller is found in the script with LQR and is used to stabilize the simulation. The result is then plotted and saved. The script, called "Experimental_Data.py", initializes by importing data collected from the real experiment into one list. The list is then manipulated to plot the desired data.
The third script "Simulation_and_Experiment_Plots.py" combines the previous two scripts and plots the simulation along with the experiment.
The file "data (1).txt" contains the experimental data obtained using the K-matrix:
    K = [707.10678119, -873.4060921, 344.75201939, -145.05769211]
The file "data (2).txt" contains the experimental data obtained using the K-matrix:
    K = [-2236.0679775, -963.92838115, 2088.64485917, 348.10668719]
The file "data (3).txt" contains the experimental data obtained using the K-matrix:
    K = [-223.607,-134.153,440.240,72.80]

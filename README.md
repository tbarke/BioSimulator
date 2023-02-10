# BioSimulator
Biological Organism Simulator for a one-dimensional environment that emphasizes in Information metrics (mutual information, subjective information)

Main Run Script is testbench.py that will run default settings for a simulation and create figures for the outputted data

for a generic explanation of the simulation and cell strategy, see below

The main simulation settings, cell strategy types, enviornments, stress rates





-----------------------------------
Cell Simulation:
The simulation consists of cells that have the ability move left or right in their two dimensional enviornment. Their "goal" or objective is what the cell's
are judged on. In this case the cell's growth rate is the objective. This growth rate is estimated by the expression sum over t = [0, T) of (log_2\frac{p_(t+1)}{p(t)}), 
where T is the length of the simulation in time p(t) is the 

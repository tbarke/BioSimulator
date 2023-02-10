# BioSimulator
Biological Organism Simulator for a one-dimensional environment that emphasizes in Information metrics (mutual information, subjective information)

Main Run Script is testbench.py that will run default settings for a simulation and create figures for the outputted data

for a generic explanation of the simulation and cell strategy, see below

The main simulation settings, cell strategy types, enviornments, stress rates are set in testbench.py. It will simulate it for a certain time and record and save the cell locations, concentrations, bound receptors, cell counts, and other parameters in the Data folder according to the data and cell run name. It will then take in cell signal information and compute the informaiton rate estimates and cell growth estimates. Finally it will output the results in the Results folder into .pdfs for the simulation parameters.

Other functions possible create a .gif file of the simulation concentrations and cells over simulated time, create a heat map of internal molecule counts, etc.

-----------------------------------
Cell Simulation:
The simulation consists of cells that have the ability move left or right in their two dimensional periodic enviornment. Their "goal" or objective is what the cell's
are judged on. In this case the cell's growth rate is the objective. This growth rate is estimated by the expression sum over t = [0, T] of (log_2(p_(t+1))(p(t))), 
where T is the length of the simulation in time p(t) is the population count of the cells.

The cells are able to move according to their strategy and bound receptors defined by a bionomal distribution, absorb molecules from their enviornment, divide once they reach a particular internal molecule threshold, and die if either of their internal molecule counts reach 0 according to some survival cost. 

Their are two molecule types in the enviornment that the cell can respond to and both molecules are required for cell propogation and division.

What the simulation is attempting to measure is if the cell strategy that has a larger growth rate has a lower pure information rate from its enviornment.

Information Rate:
The information rate as seen by the cells expresses how much the signals from the enviornemnt has control over the cell's response. It is estimated in a few ways considering how the signal is defined. All estimations use mutual information or MI. This is a mathematical measure that defines the correlation and information rate between two sets of data. One measure is estimated through the MI between the external cell concentrations and their bound receptors. Another measure is estimated through the MI between the external cell concentrations and the cell response (i.e. the movement). Another type is similar to the above estimates, but only considers the gradient of the external concentrations or the signal differences corresponding to either side of the cell.

Cell Stress:
The cell stress is defined as the survival cost of the cell per time divided by the absorption rate of the cell. The higher the cell stress is the harder it is for them to survive.

Different Cell Enviornments:
The cells operate in different enviornments. One is defined by the von Mises distribution (close to a normal distribution but with a periodic variance) Another is the constant uphill downhill enviornment which applies a constant gradient to the cell at any given time. Another is the constant conentration which applies a constant concentration everywhere in the enviornment. The last is the Mana From Heaven model which is created by randomly dropping concentraitons of molecule types that diffuse and degrade at a particular rate.

Cell Strategies:


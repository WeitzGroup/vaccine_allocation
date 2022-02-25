# vaccine_allocation
A1. Figures - Main Text 

Run “plot_all.m” for generating all figures, if mat files are present. To generate mat files, run "main_all.m"

Figure 1: 
Dependence of the optimal fraction of vaccines to be donated to κ (coupling coefficient) and the vaccination rate. 


Figure 2: 
Variation of fatalities per 10^7 over 1 year in countries A and B for different values of µ. The value of µ which minimizes the fatalities in country A is termed the ‘optimal policy’. The simulation is run for low, medium and high coupling coefficients, κ (with κ ∈ {10^−6 , 10^−4, 10^−2}) and the daily vaccination rate is 0.28% of the total population

Figure 3: 
Plots of the fatalities per 107 over 1 year for the no-sharing, optimal and hybrid policies in countries A and B. The analysis is done for different coupling coefficients (κ ∈ [10^−8 , 10^−1]) and three vaccination rates (50%, 75% and 100% of the population vaccinated in 1 year, assuming vaccine stock lasts). 

A2. Figures - Supplementary Information

Figure S1: 
Fatalities in countries A and B when the (a) no-sharing policy, (b) optimal policy and (c) hybrid policy is implemented, over different κ ∈[10^−8,10^−1] and vaccination rates (from 0.14% to 0.42% of the population daily).

Figure S2:
Change in fatalities in countries A and B when comparing the hybrid policy with the optimal policy, over different κ ∈ [10−8,10−1] and vaccination rates (from 0.14% to 0.42% of the population daily).

B1. Description of key scripts - Optimal Control 

1. 	dLdx.m: implementation of dL/dx 
2. 	dLdu.m: implementation of dL/du 
3. 	dJdu.m: implementation of dJ/du
4. 	dfdx.m: implementation of df/dx
5.	dfdu.m: implementation of df/du
6. 	normdJdu: implementation of |dJ/du|
6. 	costate_f.m: costate system
7. 	costate_solver.m: numerical integration of costate system based on Euler method (backward direction)
8. 	state_f.m: state dynamical system
9. 	state_solver.m: numerical integration of state system based on Euler method (forward direction)
10. 	cost.m: cost functional
11. 	proj_operator: piecewise constant function approximation of continuous control input
12: 	L2innerProduct: inner product (in L2 space) of two signals 


B2. Description of key scripts - Feedback Control 

1.	bang_bang_control.m: Takes input of system state and provides the contact rate for susceptibles based on line params.
2.	compute_man_days_fraction.m: Based on contact rates and states, computes person-hours worked during lockdown
3.	ga_heuristic.m: Generates mat file with optimal line params for I-R plane to maximize cost function using inbuilt genetic algorithm function
4.	heuristic_applied_grid.m: Computes state trajectories based on optimal line parameters 
5.	heuristic_cost_new: Cost function to be maximized using GA
6.	plot_ga_full_grid: outputs 11 figures which make Fig. S9
7.	plot_ga_reduced_grid: outputs 3 figures which make Fig. 4

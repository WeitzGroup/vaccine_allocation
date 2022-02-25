# vaccine_allocation
A1. Figures - Main Text 

Run “plot_all.m” for generating all figures, if mat files are present. To generate mat files, run "main.m"

Figure 1: 
Dependence of the optimal fraction of vaccines to be donated to κ (coupling coefficient) and the vaccination rate. 


Figure 2: 
Variation of fatalities per 10^7 over 1 year in countries A and B for different values of µ. The value of µ which minimizes the fatalities in country A is termed the ‘optimal policy’. The simulation is run for low, medium and high coupling coefficients, κ (with κ ∈ {10^−6 , 10^−4, 10^−2}) and the daily vaccination rate is 0.28% of the total population. For plotting, change line 50 in plot_all.m with appropriate value of κ, as needed ({10^−6 , 10^−4, 10^−2})

Figure 3: 
Plots of the fatalities per 107 over 1 year for the no-sharing, optimal and hybrid policies in countries A and B. The analysis is done for different coupling coefficients (κ ∈ [10^−8 , 10^−1]) and three vaccination rates (50%, 75% and 100% of the population vaccinated in 1 year, assuming vaccine stock lasts). 

A2. Figures - Supplementary Information

Figure S1: 
Fatalities in countries A and B when the (a) no-sharing policy, (b) optimal policy and (c) hybrid policy is implemented, over different κ ∈[10^−8,10^−1] and vaccination rates (from 0.14% to 0.42% of the population daily).

Figure S2:
Change in fatalities in countries A and B when comparing the hybrid policy with the optimal policy, over different κ ∈ [10−8,10−1] and vaccination rates (from 0.14% to 0.42% of the population daily).

B. Description of key scripts

1. 	force_of_infection.m: implementation of force of infection for 2 coupled systems 
2. 	state_f.m: state dynamical system 
3. 	state_solver.m: numerical integration of state system based on Euler method (forward direction)
4. 	main.m: code to find optimal sharing fraction for different vaccination rates and coupling factor; creates data files (.mat) for plotting
5. 	plot_all.m: plot all graphs using saved data

%% Parameters
% model parameters
parsM.Tinc = 4; % length of incubation period
parsM.Tinf = 6; % duration patient is infectious
parsM.etaI = 0.1; % *true* transmission effectivenesss, note transmission rate beta ~ etaI*contactRate
parsM.mu = 1e-2; % case fatality ratio
parsM.c_baseline = 5; % baseline contact rate
parsM.Ntot = 1e7; % total number of population (neglect death)
parsM.numSVar = 10; % number of state variables
parsM.numCVar = 1; % number of control variables

% contact rates - c_S, c_E, c_I, c_R, c_V
parsM.cA = [parsM.c_baseline; parsM.c_baseline; parsM.c_baseline/2; ...
    parsM.c_baseline; parsM.c_baseline];
parsM.cB = [parsM.c_baseline; parsM.c_baseline; parsM.c_baseline/2; ...
    parsM.c_baseline; parsM.c_baseline];

parsM.kappa = 0;

parsM.total_vaccines = 0.7*parsM.Ntot;

% simulation params
parsS.idx = 1;
parsS.step = 0.02;

usa_vac_rate = 0.5/(6*30); % USA ~6 months for 50% fully vaccinated
parsS.vaccination_rate_baseline = 1 * usa_vac_rate * parsM.Ntot;                           
% numerical solver parameters;
parsT.dt = 1e-1;

%% Simulation of outbreak
% outbreak simulation with true model
parsT.t0 = 0;
parsT.tf = 12*30;   % 12 months

ini_infected_1_base = 500;
ini_infected_2_base = 500;

%ini_infected_multiplier_vector = [0.1, 0.3, 0.5, 1,  2,  5, 10];
ini_infected_multiplier_vector = 1;

% costate boundary condition
lambda_tf = zeros(12,1);
lambda_tf(6) = 1;       % grad of D_A

initial_state.A = [parsM.Ntot - ini_infected_1_base, 0, ...
    ini_infected_1_base, 0, 0, 0];
initial_state.B = [parsM.Ntot - ini_infected_2_base, 0, ...
    ini_infected_2_base, 0, 0, 0];
%% kappa vector
kappa_iter =1;
kappa_vec = 0;
for i = 8:-1:2
    for j = 1:9
        kappa_vec(kappa_iter) = j * 10^(-i);
        kappa_iter = kappa_iter +1;
    end    
end
kappa_vec(kappa_iter) = 10^(-1);
lambda = usa_vac_rate * parsM.Ntot;
vac_donated_coarse_opt = zeros(length(0.5*lambda:0.05*lambda:1.5*lambda),length(kappa_vec));

%% coarse search for start point
mu_vec = 0:0.03:1;
vac_rate_idx = 1;
if ~isfile('coarse_start_pt_10^7_pop.mat')
    for vaccination_rate = 0.5*lambda:0.05*lambda:1.5*lambda 
    parsS.vaccination_rate_baseline = vaccination_rate;
    vac_rate_idx
    kappa_idx = 1;
        for kappa = kappa_vec
            parsM.kappa = kappa;
            for i =1:length(mu_vec)
                parsS.VA = parsM.total_vaccines * (1-(mu_vec(i)));
                parsS.VB = parsM.total_vaccines * (mu_vec(i));
                state_sol_test= state_solver(parsM, parsT, parsS,initial_state);
                deaths_A(i) = state_sol_test.A(end,end);
                deaths_B(i) = state_sol_test.B(end,end);
            end
            [min_death_A, idx] = min(deaths_A);
            mu_optimal_val = mu_vec(idx);
            vac_donated_coarse_opt(vac_rate_idx, kappa_idx) = mu_optimal_val;
            kappa_idx = kappa_idx + 1;
        end
        save('coarse_start_pt_10^7_pop', 'vac_donated_coarse_opt')
        vac_rate_idx = vac_rate_idx + 1;
    end
    save('coarse_start_pt_10^7_pop', 'vac_donated_coarse_opt')
end
%% gradient desc for optimal
load('coarse_start_pt_10^7_pop.mat')
vac_rate_idx = 1;
for vaccination_rate = 0.5*lambda:0.05*lambda:1.5*lambda 
    parsS.vaccination_rate_baseline = vaccination_rate;
    vac_rate_idx
    
kappa_idx = 1;
%kappa_vec = [10^(-8), 10^(-7),10^(-6),10^(-5), 10^(-4), 10^(-3),...
%    10^(-2), 10^(-1)];


for kappa = kappa_vec
parsM.kappa = kappa;
kappa
dv = 0.02;
dm = 0.004;

vac_donated = vac_donated_coarse_opt(vac_rate_idx, kappa_idx); % use coarse opt as start pt for gradient desc
grad_prev = 10^5;
grad = 10^5;
parsS.step = 0.02;
while abs(grad) > 0.1  
    parsS.VA = parsM.total_vaccines * (1-vac_donated);
    parsS.VB = parsM.total_vaccines * vac_donated;
    
    state_sol = state_solver(parsM, parsT, parsS, initial_state);
    J1 = state_sol.A(end, end);
  
    parsS.VA = parsM.total_vaccines * (1-(vac_donated + dv));
    parsS.VB = parsM.total_vaccines * (vac_donated + dv);
    state_sol_test = state_solver(parsM, parsT, parsS, initial_state);
    J2 = state_sol_test.A(end,end);
    grad_prev = grad;
    grad = (J2-J1)/dv;
    
    if grad*grad_prev < 0
        parsS.step = parsS.step/2;
    end
    
    if abs(grad) > 1
        grad = grad/abs(grad);
    end    
    vac_donated = vac_donated - parsS.step*grad;
   
    if vac_donated >= 0.5
        vac_donated = 0.5;
        break
    end
    if vac_donated <= 0
        vac_donated = 0;
        break
    end
end  
    
    vac_don_save(vac_rate_idx, kappa_idx) =  vac_donated;
    kappa_idx = kappa_idx + 1;
        
end 
    save('vac_don_save_70_perc_10^7_pop', 'vac_don_save')
    vac_rate_idx = vac_rate_idx + 1;
end    
%%
figure(1)
for i = 1:vac_rate_idx-1
semilogx(kappa_vec, vac_don_save(i, :),'Linewidth', 3);
hold on
%val = ini_frac_min + (i-1)*ini_frac_delta;
%title("ini fraction = " + val )
%ylim([0, 0.4])
end
xlabel('Kappa', 'FontName', 'Times New Roman','FontSize',20,'Interpreter','latex'); 
ylabel('Optimal fraction', 'FontName', 'Times New Roman','FontSize',20, 'Interpreter','latex');
%%
save('vac_don_save_70_perc_10^7_pop.mat', 'vac_don_save')



%% comparison of baseline, 1/3, 1/2 and optimal strats
fatalities_baseline.A = 0 * vac_don_save;   %  no share strat
fatalities_baseline.B = fatalities_baseline.A;

fatalities_optimal = fatalities_baseline;   % optimal share strat

fatalities_third_share = fatalities_baseline;   % mu = 1/3 strat

fatalities_half_share = fatalities_baseline; % mu = 0.5 strat

vac_rate_vec = 0.5*lambda:0.05*lambda:1.5*lambda;

for i = 1:length(vac_rate_vec)
    for j = 1:length(kappa_vec)
        parsM.kappa = kappa_vec(j);
        parsS.vaccination_rate_baseline = vac_rate_vec(i);
        vac_donated = [0, 1/3, 1/2, vac_don_save(i, j)];
        for k = 1:4    
            % setup for no share
            parsS.VA = parsM.total_vaccines * (1-vac_donated(k));
            parsS.VB = parsM.total_vaccines * vac_donated(k);
            state_sol = state_solver(parsM, parsT, parsS, initial_state);
            if k == 1       % no share
                fatalities_baseline.A(i, j) = state_sol.A(end, end);
                fatalities_baseline.B(i, j) = state_sol.B(end, end);
            elseif k == 2   % 1/3 share    
                fatalities_third_share.A(i, j) = state_sol.A(end, end);
                fatalities_third_share.B(i, j) = state_sol.B(end, end);
            elseif k ==3    % 1/2 share
                fatalities_half_share.A(i, j) = state_sol.A(end, end);
                fatalities_half_share.B(i, j) = state_sol.B(end, end);
            elseif k == 4   % optimal share
                fatalities_optimal.A(i, j) = state_sol.A(end, end);
                fatalities_optimal.B(i, j) = state_sol.B(end, end);
            end        
        end    
    end
end  


%%
fatalities_all.baseline = fatalities_baseline;
fatalities_all.optimal = fatalities_optimal;
fatalities_all.third = fatalities_third_share;
fatalities_all.half = fatalities_half_share;

save('strategy_compare_10^7_pop.mat', 'fatalities_all');

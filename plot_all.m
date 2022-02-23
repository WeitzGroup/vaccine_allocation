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

parsM.vaccine_reduction = 0.1;         %90 perc efficacy
parsM.total_vaccines = 0.7*parsM.Ntot;

% optimization algorithm parameters (backtracking step size)
parsC.alpha = 0.1; % alpha is *very* important 
parsC.beta = 0.5;

% simulation params
parsS.idx = 1;
parsS.step = 0.02;
                          

% numerical solver parameters;
parsT.dt = 1e-1;
parsT.t0 = 0;
parsT.tf = 12*30;   % 12 months


ini_infected_1_base = 500;
ini_infected_2_base = 500;
initial_state.A = [parsM.Ntot - ini_infected_1_base, 0, ...
    ini_infected_1_base, 0, 0, 0];
initial_state.B = [parsM.Ntot - ini_infected_2_base, 0, ...
    ini_infected_2_base, 0, 0, 0];

usa_vac_rate = 0.5/(6*30);
lambda = usa_vac_rate * parsM.Ntot;

%% reduction in fatalities for different mu
mu_vec = 0:0.01:1;
parsM.kappa = 10^(-6);
parsS.vaccination_rate_baseline = 1*lambda;

for i =1:length(mu_vec)
    parsS.VA = parsM.total_vaccines * (1-(mu_vec(i)));
    parsS.VB = parsM.total_vaccines * (mu_vec(i));
    state_sol_test= state_solver(parsM, parsC, parsT, parsS,initial_state);
    deaths_A(i) = state_sol_test.A(end,end);
    deaths_B(i) = state_sol_test.B(end,end);
end
%%
[min_death_A, idx] = min(deaths_A);
mu_optimal_val = mu_vec(idx);

figure;
plot(mu_vec, deaths_A/(parsM.Ntot/10^7), 'b','Linewidth', 3)
hold on
plot(mu_vec, deaths_B/(parsM.Ntot/10^7), 'r','Linewidth', 3)
%plot(mu_vec, deaths_A + deaths_B, '--', 'Linewidth', 3)
yline(deaths_A(1)/(parsM.Ntot/10^7), '--b', 'Linewidth', 3)
yline(deaths_B(1)/(parsM.Ntot/10^7), '--r', 'Linewidth', 3)
plot(mu_optimal_val ,deaths_A(idx)/(parsM.Ntot/10^7), 'o', 'MarkerSize',10,...
    'MarkerEdgeColor','blue', 'Linewidth', 2)
plot(mu_optimal_val ,deaths_B(idx)/(parsM.Ntot/10^7), 'o', 'MarkerSize',10,...
    'MarkerEdgeColor','red', 'Linewidth', 2)

plot(0.33 ,deaths_A(34)/(parsM.Ntot/10^7), 'd', 'MarkerSize',10,...
    'MarkerEdgeColor','blue', 'Linewidth', 2)
plot(0.33 ,deaths_B(34)/(parsM.Ntot/10^7), 'd', 'MarkerSize',10,...
    'MarkerEdgeColor','red', 'Linewidth', 2)

reduction_A_optimal = (deaths_A(1) - deaths_A(idx))/deaths_A(1)*100;
reduction_B_optimal = (deaths_B(1) - deaths_B(idx))/deaths_B(1)*100;

reduction_A_33 = (deaths_A(1) - deaths_A(34))/deaths_A(1)*100;
reduction_B_33 = (deaths_B(1) - deaths_B(34))/deaths_B(1)*100;

%str = {'$\uparrow $Optimal strategy'};
%text(mu_optimal_val, deaths_A(idx),'Optimal fraction$\downarrow $', 'FontSize', 20, Interpreter='latex')

fatalities_A = compose('%.1f',reduction_A_optimal);
str1 = {fatalities_A{1},'\% reduction in A'};
str1 = strjoin(str1, {'\b'});

fatalities_B = compose('%.1f',reduction_B_optimal);
str2 = {fatalities_B{1},'\% reduction in B'};
str2 = strjoin(str2, {'\b'});
% print fatality reduction
%text(0, 0, str1, 'FontSize', 20, Interpreter='latex')
%text(0, 10000, str2, 'FontSize', 20, Interpreter='latex')


xlim([0, 0.5])

axis square

legend('Fatalities in A ','Fatalities in B','Fatalities in A (no sharing)',...
    'Fatalities in B (no sharing)','Fatalities in A (optimal)',...
    'Fatalities in B (optimal)','Fatalities in A (33\% sharing)',...
    'Fatalities in B (33\% sharing)','Interpreter','latex')
xlabel('$\mu$ (Fraction donated)','Interpreter','latex')
ylabel('Fatalities per $10^7$ over 1 year','Interpreter','latex')
set(gca,'FontSize',20);
legend boxoff
title('$\kappa = 10^{-2}$')



%% graphs for optimal strats
% heatmap showing optimal mu for different kappa and vaccination rates
%load('vac_don_save_70_perc_sym_diff_extended_final.mat')  %for 10^6
load('vac_don_save_70_perc_10^7_pop.mat')

kappa_iter =1;
kappa_vec = 0;
for i = 8:-1:2
    for j = 1:9
        kappa_vec(kappa_iter) = j * 10^(-i);
        kappa_iter = kappa_iter +1;
    end    
end
kappa_vec(kappa_iter) = 10^(-1);
figure;
%xvals = 0.5*lambda:0.1*lambda:1.5*lambda;
%xvals = {'0.5\lambda_0', '0.6\lambda_0','0.7\lambda_0','0.8\lambda_0','0.9\lambda_0'...
%    ,'\lambda_0','1.1\lambda_0','1.2\lambda_0','1.3\lambda_0','1.4\lambda_0','1.5\lambda_0'};
%yvals = {10^(-8), 10^(-7),10^(-6),10^(-5), 10^(-4), 10^(-3),...
%    10^(-2), 10^(-1)};
%yvals = {'R_0 = 40','R_0 = 60','R_0 = 80','R_0 = 100'};
xvals = 0.5:0.05:1.5;
xvals = xvals * 100/360;
yvals = kappa_vec;
h = heatmap(xvals, yvals, vac_don_save','CellLabelColor','none');
%h = heatmap(vac_don_save','CellLabelColor','none');
h.Colormap = parula;

xlabel('Percentage vaccinated daily ($\%$)')
ylabel('Coupling coefficient, $\kappa$')
set(gca,'FontSize',20);
h.XDisplayLabels = compose('%.2f',str2double(h.XDisplayLabels));
for i = 1:3:63
    yvals_cells{1, i} = yvals(i);
    yvals_cells{1, i + 1} = "";
    yvals_cells{1, i + 2} = "";
end    
for i = 2:2:21
    xvals_cells{1, i-1} = xvals(i-1);
    xvals_cells{1, i} = "";
end
xvals_cells{1, 21} = xvals(21);

yvals_cells{1, 64} = yvals(64);
h.YDisplayLabels = yvals_cells;
h.NodeChildren(3).XAxis.Label.Interpreter = 'latex';
h.NodeChildren(3).YAxis.Label.Interpreter = 'latex';
h.NodeChildren(3).Title.Interpreter = 'latex';
h.NodeChildren(3).YDir='normal'; 

figure;
contour(vac_don_save', [0.33, 0.33],'--k', 'LineWidth',3)
hold on
contour(vac_don_save', [3*10^(-3),3*10^(-3)],'--r', 'LineWidth',3)

xticks(1:2:21)
xticklabels(xvals(1:2:21)*360)
title('No sharing', Interpreter='latex')
title('$33\%$ sharing', Interpreter='latex')
xlabel('Percentage vaccinated yearly ($\%$)',Interpreter='latex')
ylabel('Coupling coefficient, $\kappa$')
set(gca,'FontSize',20);

% 
% figure;
% yyaxis right
% plot(1:7, -20:20:100)
% yticks(-20:20:100)
% yticklabels({'-20\%','0\%','20\%','40\%','60\%','80\%','100\%'})
% set(groot,'defaultAxesTickLabelInterpreter','latex');  
% set(gca,'FontSize',20);
% ax = gca;
% ax.YAxis(2).Color = 'k';
%% compute optimal v/s no share - change in fatalities in A and B (create matrix)
%load('strategy_compare.mat') % for 10^6 pop
load('strategy_compare_10^7_pop.mat')

% fatality change optimal 
fatality_reduction_optimal.A = 100*(fatalities_all.baseline.A - fatalities_all.optimal.A)./fatalities_all.baseline.A;
fatality_reduction_optimal.B = 100*(fatalities_all.baseline.B - fatalities_all.optimal.B)./fatalities_all.baseline.B;

fatality_reduction_third.A = 100*(fatalities_all.baseline.A - fatalities_all.third.A)./fatalities_all.baseline.A;
fatality_reduction_third.B = 100*(fatalities_all.baseline.B - fatalities_all.third.B)./fatalities_all.baseline.B;

fatality_reduction_half.A = 100*(fatalities_all.baseline.A - fatalities_all.half.A)./fatalities_all.baseline.A;
fatality_reduction_half.B = 100*(fatalities_all.baseline.B - fatalities_all.half.B)./fatalities_all.baseline.B;

%%
% heatmaps for optimal
plot_heatmap(xvals, yvals, fatality_reduction_optimal.A, -20, 100)
plot_heatmap(xvals, yvals, fatality_reduction_optimal.B, -20, 100);

% heatmaps for 1/3 reduction
plot_heatmap(xvals, yvals, fatality_reduction_third.A,-20, 100);
plot_heatmap(xvals, yvals, fatality_reduction_third.B,-20, 100);

% heatmaps for 1/2 reduction
plot_heatmap(xvals, yvals, fatality_reduction_half.A,-20, 100);
plot_heatmap(xvals, yvals, fatality_reduction_half.B,-20, 100);

%% hybrid strat
% heatmap comparing hybrid strat (max(1/3, optimal)) with optimal

fatality_hybrid = fatalities_all.third;
fatality_reduction_hybrid_compare_optimal = fatality_reduction_third;
fatality_reduction_hybrid = fatality_reduction_third;     

for i =1:length(vac_don_save(:, 1))
    for j = 1:length(vac_don_save(1,:))
        if vac_don_save(i, j) > 1/3
            fatality_hybrid.A(i,j) = fatalities_all.optimal.A(i,j);
            fatality_hybrid.B(i,j) = fatalities_all.optimal.B(i,j);

            fatality_reduction_hybrid_compare_optimal.A(i,j) = 0;
            fatality_reduction_hybrid_compare_optimal.B(i,j) = 0;

            fatality_reduction_hybrid.A(i,j) = fatality_reduction_optimal.A(i,j);
            fatality_reduction_hybrid.B(i,j) = fatality_reduction_optimal.B(i,j);

        else
            fatality_reduction_hybrid_compare_optimal.A(i,j) = (fatality_hybrid.A(i,j)-fatalities_all.optimal.A(i,j))/fatalities_all.optimal.A(i,j);
            fatality_reduction_hybrid_compare_optimal.B(i,j) = (fatality_hybrid.B(i,j)-fatalities_all.optimal.B(i,j))/fatalities_all.optimal.B(i,j);
        end
    end    
end    

% compare hybrid to no-share
plot_heatmap(xvals, yvals, fatality_reduction_hybrid.A,-20, 100);
plot_heatmap(xvals, yvals, fatality_reduction_hybrid.B,-20, 100);

% compare hybrid to optimal
plot_heatmap(xvals, yvals, 100*fatality_reduction_hybrid_compare_optimal.A, 0, 100);
plot_heatmap(xvals, yvals, -100*fatality_reduction_hybrid_compare_optimal.B, 0, 100);

%% actual fatalities in no-share, optimal and hybrid
% no-share
plot_heatmap(xvals, yvals, fatalities_all.baseline.A, 0, 8*10^4);
plot_heatmap(xvals, yvals, fatalities_all.baseline.B, 0, 8*10^4);

% optimal
plot_heatmap(xvals, yvals, fatalities_all.optimal.A, 0, 8*10^4);
plot_heatmap(xvals, yvals, fatalities_all.optimal.B, 0, 8*10^4);

% hybrid
plot_heatmap(xvals, yvals, fatality_hybrid.A, 0, 8*10^4);
plot_heatmap(xvals, yvals, fatality_hybrid.B, 0, 8*10^4);

%%
figure;
yyaxis right
plot(0:8, 0:10^4:8*10^4)
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
title('Fatalities in A', Interpreter='latex')
set(gca,'FontSize',20);

figure;
title('Fatalities in B', Interpreter='latex')
set(gca,'FontSize',20);

%% graph of optimal reduction and hybrid reduction for a and b for vac rate of 100% per year
for i = 1:5:11
figure;
semilogx(kappa_vec, fatalities_all.baseline.A(i,:), ':b', 'Linewidth', 3)
hold on
semilogx(kappa_vec, fatalities_all.baseline.B(i,:), ':r', 'Linewidth', 3)
semilogx(kappa_vec, fatalities_all.optimal.A(i, :),'b','Linewidth', 2)
semilogx(kappa_vec, fatalities_all.optimal.B(i, :),'r','Linewidth', 2)

plot(kappa_vec, fatality_hybrid.A(i, :), '--b','Linewidth', 3)
plot(kappa_vec, fatality_hybrid.B(i, :), '--r','Linewidth', 3)
xline(10^(-4),'-.k','Linewidth', 2)
xticks([10^(-8),10^(-7),10^(-6),10^(-5),10^(-4),10^(-3),10^(-2),10^(-1)])               
%yticks(-20:10:100)
legend('Fatalities in A (No-sharing)','Fatalities in B (No-sharing)','Fatalities in A (Optimal)','Fatalities in B (Optimal)',...
    'Fatalities in A (Hybrid)','Fatalities in B (Hybrid)','Interpreter','latex')
legend boxoff
ax = gca;
ax.XAxisLocation = 'origin';    
ylabel('Fatalities per $10^7$ in 1 year', interpreter='latex')
%ylabel({'Excess';'deaths(\%)'}, interpreter='latex')
xlabel('Coupling coefficient, $\kappa$', interpreter='latex')
ylim([0,10^5])
xlim([10^(-8), 10^(-1)])
set(gca,'FontSize',20);
text(10^(-5), 20, '$\hat{\mu} = 0.33$', 'FontSize', 20, Interpreter='latex')
text(10^(-5), 60, '$\hat{\mu} = \mu^*$', 'FontSize', 20, Interpreter='latex')

axis square
    if i ==1
        title('$50\%$ vaccinated in 1 year', Interpreter='latex')
    elseif i == 6
        title('$75\%$ vaccinated in 1 year', Interpreter='latex')
    else
        title('$100\%$ vaccinated in 1 year', Interpreter='latex')
    end   

end


%% graph of optimal and hybrid fatalities for a and b for vac rate of 100% per year
for i = 1:5:11
figure;
semilogx(kappa_vec, fatality_reduction_optimal.A(i, :),'b','Linewidth', 2)
hold on 
semilogx(kappa_vec, fatality_reduction_optimal.B(i, :),'r','Linewidth', 2)

plot(kappa_vec, fatality_reduction_hybrid.A(i, :), '--b','Linewidth', 2)
plot(kappa_vec, fatality_reduction_hybrid.B(i, :), '--r','Linewidth', 2)
xline(10^(-4),'-.k','Linewidth', 2)
xticks([10^(-8),10^(-7),10^(-6),10^(-5),10^(-4),10^(-3),10^(-2),10^(-1)])               
yticks(-20:10:100)
legend('Fatality reduction in A (Optimal)','Fatality reduction in B (Optimal)',...
    'Fatality reduction in A (Hybrid)','Fatality reduction in B (Hybrid)','Interpreter','latex')
legend boxoff
ax = gca;
ax.XAxisLocation = 'origin';    
%ylabel('Reduction in deaths ($\%$)', interpreter='latex')
ylabel({'Excess';'deaths(\%)'}, interpreter='latex')
xlabel('Coupling coefficient, $\kappa$', interpreter='latex')
ylim([-20,100])
xlim([10^(-8), 10^(-1)])
set(gca,'FontSize',20);
text(10^(-5), 20, '$\hat{\mu} = 0.33$', 'FontSize', 20, Interpreter='latex')
text(10^(-5), 60, '$\hat{\mu} = \mu^*$', 'FontSize', 20, Interpreter='latex')

axis square
    if i ==1
        title('$50\%$ vaccinated in 1 year', Interpreter='latex')
    elseif i == 6
        title('$75\%$ vaccinated in 1 year', Interpreter='latex')
    else
        title('$100\%$ vaccinated in 1 year', Interpreter='latex')
    end   

end


%%


%{
%%
    load('vac_don_save_70_perc.mat')

kappa_vec = [10^(-8), 10^(-7),10^(-6),10^(-5), 10^(-4), 10^(-3),...
    10^(-2), 10^(-1)];

vaccination_iter = 1;
for vaccination_rate = 0.5*lambda:0.1*lambda:1.5*lambda 
    parsS.vaccination_rate_baseline = vaccination_rate;
    kappa_iter = 1;
    for kappa = kappa_vec
        parsM.kappa = kappa;
        vac_donated = vac_don_save(vaccination_iter, kappa_iter);
        parsS.VA = parsM.total_vaccines * (1-(vac_donated));
        parsS.VB = parsM.total_vaccines * (vac_donated);
        state_sol_optimal = state_solver(parsM, parsC, parsT, parsS, ...
            initial_state);
        death_optimal_A(vaccination_iter, kappa_iter) = ...
            state_sol_optimal.A(end, end);
        death_optimal_B(vaccination_iter, kappa_iter) = ...
            state_sol_optimal.B(end, end);

        parsS.VA = parsM.total_vaccines * (1);
        parsS.VB = parsM.total_vaccines * (0);
        state_sol_no_donation = state_solver(parsM, parsC, parsT, parsS,...
            initial_state);
        death_no_donation_A(vaccination_iter, kappa_iter) = ...
            state_sol_no_donation.A(end, end);
        death_no_donation_B(vaccination_iter, kappa_iter) = ...
            state_sol_no_donation.B(end, end);

        kappa_iter = kappa_iter + 1;
    end    
    vaccination_iter = vaccination_iter +1;
end    

%%
start_val = 5;
figure;
semilogx(kappa_vec(start_val:end), death_no_donation_A(2, start_val:end))
hold on
semilogx(kappa_vec(start_val:end), death_optimal_A(2, start_val:end))

death_reduction_fraction_matrix_A = (death_no_donation_A - death_optimal_A)...
    ./death_no_donation_A;

death_reduction_fraction_matrix_B = (death_no_donation_B - death_optimal_B)...
    ./death_no_donation_B;

figure(77);

for i = 5:8
    plot(0.5*lambda:0.1*lambda:1.5*lambda, ...
        death_reduction_fraction_matrix_A(:,i)*100, ...
        'Linewidth', 3)
    hold on
end    

legend('$\kappa=10^{-4}$','$\kappa=10^{-3}$','$\kappa=10^{-2}$',...
    '$\kappa=10^{-1}$','Interpreter','latex')
xlabel('$Daily~vaccination~rate~(/day)$','Interpreter','latex')
ylabel('$Reduction~in~fatalities~(\%)$','Interpreter','latex')
set(gca,'FontSize',20);
axis square

figure(78);

for i = 5:8
    if i == 5
        plot(0.5*lambda:0.1*lambda:1.3*lambda, ...
        death_reduction_fraction_matrix_B(1:9,i)*100, ...
        'Linewidth', 3)
        hold on
    else
        plot(0.5*lambda:0.1*lambda:1.5*lambda, ...
        death_reduction_fraction_matrix_B(:,i)*100, ...
        'Linewidth', 3)
    hold on
    end
end    
axis square
%legend('$\kappa=10^{-6}$','$\kappa=10^{-5}$','$\kappa=10^{-4}$',...
%     '$\kappa=10^{-3}$','$\kappa=10^{-2}$',...
%     '$\kappa=10^{-1}$','Interpreter','latex')
legend('$\kappa=10^{-4}$',...
    '$\kappa=10^{-3}$','$\kappa=10^{-2}$',...
    '$\kappa=10^{-1}$','Interpreter','latex')    
xlabel('$Daily~vaccination~rate~(/day)$','Interpreter','latex')
ylabel('$Reduction~in~fatalities~(\%)$','Interpreter','latex')
set(gca,'FontSize',20);


%%
parsS.vaccination_rate_baseline = lambda;
parsM.kappa = 0.01;
parsS.VA = parsM.total_vaccines * (1);
parsS.VB = parsM.total_vaccines * (0);
state_sol_no_donation_sample = state_solver(parsM, parsC, parsT, parsS,...
            initial_state);

vac_donated_sample = vac_don_save(6, 7);
parsS.VA = parsM.total_vaccines * (1-vac_donated_sample);
parsS.VB = parsM.total_vaccines*vac_donated_sample;
state_sol_optimal_sample = state_solver(parsM, parsC, parsT, parsS,...
            initial_state);



%% sample plot for no donation vs optimal donation, pop dynamics

% no donation
figure;
semilogy(state_sol_no_donation_sample.A(:,1),...
    state_sol_no_donation_sample.A(:,2), 'Linewidth', 3)
hold on
semilogy(state_sol_no_donation_sample.A(:,1),...
    state_sol_no_donation_sample.A(:,3), 'Linewidth', 3)
semilogy(state_sol_no_donation_sample.A(:,1),...
    state_sol_no_donation_sample.A(:,4), 'Linewidth', 3)
semilogy(state_sol_no_donation_sample.A(:,1),...
    state_sol_no_donation_sample.A(:,5), 'Linewidth', 3)
semilogy(state_sol_no_donation_sample.A(:,1),...
    state_sol_no_donation_sample.A(:,6), 'Linewidth', 3)
semilogy(state_sol_no_donation_sample.A(:,1),...
    state_sol_no_donation_sample.A(:,7), 'Linewidth', 3)
ylim([10^0, 10^6])
xlim([0, 360])
xticks([0, 60, 120, 180, 240, 300, 360])
legend('$S_A$','$E_A$','$I_A$','$R_A$','$V_A$','$D_A$','Interpreter','latex')
xlabel('$Days$','Interpreter','latex')
ylabel('$Population$','Interpreter','latex')
set(gca,'FontSize',20);
axis square

figure;
semilogy(state_sol_no_donation_sample.A(:,1),...
    state_sol_no_donation_sample.B(:,2), 'Linewidth', 3)
hold on
semilogy(state_sol_no_donation_sample.A(:,1),...
    state_sol_no_donation_sample.B(:,3), 'Linewidth', 3)
semilogy(state_sol_no_donation_sample.A(:,1),...
    state_sol_no_donation_sample.B(:,4), 'Linewidth', 3)
semilogy(state_sol_no_donation_sample.A(:,1),...
    state_sol_no_donation_sample.B(:,5), 'Linewidth', 3)
semilogy(state_sol_no_donation_sample.A(:,1),...
    state_sol_no_donation_sample.B(:,6), 'Linewidth', 3)
semilogy(state_sol_no_donation_sample.A(:,1),...
    state_sol_no_donation_sample.B(:,7), 'Linewidth', 3)
ylim([10^0, 10^6])
xlim([0, 360])
xticks([0, 60, 120, 180, 240, 300, 360])
legend('$S_B$','$E_B$','$I_B$','$R_B$','$V_B$','$D_B$','Interpreter','latex')
xlabel('$Days$','Interpreter','latex')
ylabel('$Population$','Interpreter','latex')
set(gca,'FontSize',20);
axis square

% optimal donation
figure;
semilogy(state_sol_optimal_sample.A(:,1),...
    state_sol_optimal_sample.A(:,2), 'Linewidth', 3)
hold on
semilogy(state_sol_optimal_sample.A(:,1),...
    state_sol_optimal_sample.A(:,3), 'Linewidth', 3)
semilogy(state_sol_optimal_sample.A(:,1),...
    state_sol_optimal_sample.A(:,4), 'Linewidth', 3)
semilogy(state_sol_optimal_sample.A(:,1),...
    state_sol_optimal_sample.A(:,5), 'Linewidth', 3)
semilogy(state_sol_optimal_sample.A(:,1),...
    state_sol_optimal_sample.A(:,6), 'Linewidth', 3)
semilogy(state_sol_optimal_sample.A(:,1),...
    state_sol_optimal_sample.A(:,7), 'Linewidth', 3)
ylim([10^0, 10^6])
xlim([0, 360])
xticks([0, 60, 120, 180, 240, 300, 360])
legend('$S_A$','$E_A$','$I_A$','$R_A$','$V_A$','$D_A$','Interpreter','latex')
xlabel('$Days$','Interpreter','latex')
ylabel('$Population$','Interpreter','latex')
set(gca,'FontSize',20);
axis square

figure;
semilogy(state_sol_optimal_sample.A(:,1),...
    state_sol_optimal_sample.B(:,2), 'Linewidth', 3)
hold on
semilogy(state_sol_optimal_sample.A(:,1),...
    state_sol_optimal_sample.B(:,3), 'Linewidth', 3)
semilogy(state_sol_optimal_sample.A(:,1),...
    state_sol_optimal_sample.B(:,4), 'Linewidth', 3)
semilogy(state_sol_optimal_sample.A(:,1),...
    state_sol_optimal_sample.B(:,5), 'Linewidth', 3)
semilogy(state_sol_optimal_sample.A(:,1),...
    state_sol_optimal_sample.B(:,6), 'Linewidth', 3)
semilogy(state_sol_optimal_sample.A(:,1),...
    state_sol_optimal_sample.B(:,7), 'Linewidth', 3)
ylim([10^0, 10^6])
xlim([0, 360])
xticks([0, 60, 120, 180, 240, 300, 360])
legend('$S_B$','$E_B$','$I_B$','$R_B$','$V_B$','$D_B$','Interpreter','latex')
xlabel('$Days$','Interpreter','latex')
ylabel('$Population$','Interpreter','latex')
set(gca,'FontSize',20);
axis square

%% only infected cases for baseline sim of optimal and no-share strats
% simple plot
figure;
subplot(1,2,1);
plot(state_sol_optimal_sample.A(:, 1), state_sol_optimal_sample.A(:,4), ...
    'k','Linewidth', 3)
hold on 
plot(state_sol_optimal_sample.A(:, 1), state_sol_no_donation_sample.A(:,4), ...
    'k--','Linewidth', 3)
legend('$Optimal$','$No~sharing$','Interpreter','latex')
xlabel('$Days$','Interpreter','latex')
ylabel('$Infected~cases~in~A$','Interpreter','latex')
set(gca,'FontSize',20);
axis square
xlim([0, 360])
xticks([0, 60, 120, 180, 240, 300, 360])

subplot(1,2,2);
plot(state_sol_optimal_sample.A(:, 1), state_sol_optimal_sample.B(:,4), ...
    'k','Linewidth', 3)
hold on 
plot(state_sol_optimal_sample.A(:, 1), state_sol_no_donation_sample.B(:,4), ...
    'k--','Linewidth', 3)
legend('$Optimal$','$No~sharing$','Interpreter','latex')
xlabel('$Days$','Interpreter','latex')
ylabel('$Infected~cases~in~B$','Interpreter','latex')
set(gca,'FontSize',20);
axis square
xlim([0, 360])
xticks([0, 60, 120, 180, 240, 300, 360])

% semilog plot
figure;
subplot(1,2,1);
semilogy(state_sol_optimal_sample.A(:, 1), state_sol_optimal_sample.A(:,4), ...
    'k','Linewidth', 3)
hold on 
semilogy(state_sol_optimal_sample.A(:, 1), state_sol_no_donation_sample.A(:,4), ...
    'k--','Linewidth', 3)
legend('$Optimal$','$No~sharing$','Interpreter','latex')
xlabel('$Days$','Interpreter','latex')
ylabel('$Infected~cases~in~A$','Interpreter','latex')
set(gca,'FontSize',20);
axis square
xlim([0, 360])
xticks([0, 60, 120, 180, 240, 300, 360])
ylim([10^1, 10^5])

subplot(1,2,2);
semilogy(state_sol_optimal_sample.A(:, 1), state_sol_optimal_sample.B(:,4), ...
    'k','Linewidth', 3)
hold on 
semilogy(state_sol_optimal_sample.A(:, 1), state_sol_no_donation_sample.B(:,4), ...
    'k--','Linewidth', 3)
legend('$Optimal$','$No~sharing$','Interpreter','latex')
xlabel('$Days$','Interpreter','latex')
ylabel('$Infected~cases~in~B$','Interpreter','latex')
set(gca,'FontSize',20);
axis square
xlim([0, 360])
xticks([0, 60, 120, 180, 240, 300, 360])
ylim([10^1, 10^5])

% cumulative infections for optimal and no-sharing 
cumulative_inf_A_optimal = zeros(1, length(state_sol_optimal_sample.A(:,1)));

cumulative_inf_A_no_share = cumulative_inf_A_optimal;

cumulative_inf_B_optimal = cumulative_inf_A_optimal;
cumulative_inf_B_no_share = cumulative_inf_A_optimal;

for i = 2:length(state_sol_optimal_sample.A(:,1))
    cumulative_inf_A_optimal(1, i) = cumulative_inf_A_optimal(1, i-1) + ...
        state_sol_optimal_sample.A(i,3)*parsT.dt/parsM.Tinc;
    cumulative_inf_A_no_share(1, i) = cumulative_inf_A_no_share(1, i-1) + ...
        state_sol_no_donation_sample.A(i,3)*parsT.dt/parsM.Tinc;

    cumulative_inf_B_optimal(1, i) = cumulative_inf_B_optimal(1, i-1) + ...
        state_sol_optimal_sample.B(i,3)*parsT.dt/parsM.Tinc;
    cumulative_inf_B_no_share(1, i) = cumulative_inf_B_no_share(1, i-1) + ...
        state_sol_no_donation_sample.B(i,3)*parsT.dt/parsM.Tinc;
end

cumulative_inf_A_optimal = cumulative_inf_A_optimal + state_sol_optimal_sample.A(1, 4);
cumulative_inf_A_no_share = cumulative_inf_A_no_share + state_sol_no_donation_sample.A(1, 4);

cumulative_inf_B_optimal = cumulative_inf_B_optimal + state_sol_optimal_sample.B(1, 4);
cumulative_inf_B_no_share = cumulative_inf_B_no_share + state_sol_no_donation_sample.B(1, 4);

%% cumulative inf plot
% simple plot
figure;
subplot(1,2,1);
plot(state_sol_optimal_sample.A(:, 1), cumulative_inf_A_optimal, ...
    'k','Linewidth', 3)
hold on 
plot(state_sol_optimal_sample.A(:, 1), cumulative_inf_A_no_share, ...
    'k--','Linewidth', 3)
legend('Optimal sharing','No sharing','Interpreter','latex')
xlabel('Days','Interpreter','latex')
ylabel('Cumulative infected cases in A','Interpreter','latex')
set(gca,'FontSize',20);
axis square
xlim([0, 360])
%ylim([10^2, 10^6])
xticks([0, 60, 120, 180, 240, 300, 360])
legend boxoff

subplot(1,2,2);
plot(state_sol_optimal_sample.A(:, 1), cumulative_inf_B_optimal, ...
    'k','Linewidth', 3)
hold on 
plot(state_sol_optimal_sample.A(:, 1), cumulative_inf_B_no_share, ...
    'k--','Linewidth', 3)
xlabel('Days','Interpreter','latex')
ylabel('Cumulative infected cases in B','Interpreter','latex')
set(gca,'FontSize',20);
axis square
legend('Optimal sharing','No sharing','Interpreter','latex')
xlim([0, 360])
%ylim([10^2, 10^6])
xticks([0, 60, 120, 180, 240, 300, 360])
legend boxoff

%% only fatalities
% plot of deaths in A and B for optimal vs no sharing
figure;
subplot(1,2,1);
plot(state_sol_optimal_sample.A(:, 1), state_sol_optimal_sample.A(:, 7), ...
    'k','Linewidth', 3)
hold on 
plot(state_sol_optimal_sample.A(:, 1), state_sol_no_donation_sample.A(:,7), ...
    'k--','Linewidth', 3)
legend('Optimal sharing','No sharing','Interpreter','latex')
xlabel('Days','Interpreter','latex')
ylabel('Fatalities in A','Interpreter','latex')
set(gca,'FontSize',20);
axis square
xlim([0, 360])
%ylim([10^0, 10^4])
xticks([0, 60, 120, 180, 240, 300, 360])
legend boxoff

subplot(1,2,2);
plot(state_sol_optimal_sample.A(:, 1), state_sol_optimal_sample.B(:, 7), ...
    'k','Linewidth', 3)
hold on 
plot(state_sol_optimal_sample.A(:, 1), state_sol_no_donation_sample.B(:,7), ...
    'k--','Linewidth', 3)
xlabel('Days','Interpreter','latex')
ylabel('Fatalities in B','Interpreter','latex')
set(gca,'FontSize',20);
axis square
legend('Optimal sharing','No sharing','Interpreter','latex')
xlim([0, 360])
%ylim([10^0, 10^4])
xticks([0, 60, 120, 180, 240, 300, 360])
legend boxoff



%% another sample kappa = 0.1, low vaccination rate

parsS.vaccination_rate_baseline = 0.5*lambda;
parsM.kappa = 0.1;
parsS.VA = parsM.total_vaccines * (1);
parsS.VB = parsM.total_vaccines * (0);
state_sol_no_donation_sample = state_solver(parsM, parsC, parsT, parsS,...
            initial_state);

vac_donated_sample = vac_don_save(1, 8);
parsS.VA = parsM.total_vaccines * (1-vac_donated_sample);
parsS.VB = parsM.total_vaccines*vac_donated_sample;
state_sol_optimal_sample = state_solver(parsM, parsC, parsT, parsS,...
            initial_state);

%% sample plot for no donation vs optimal donation, pop dynamics

% no donation
figure;
semilogy(state_sol_no_donation_sample.A(:,1),...
    state_sol_no_donation_sample.A(:,2), 'Linewidth', 3)
hold on
semilogy(state_sol_no_donation_sample.A(:,1),...
    state_sol_no_donation_sample.A(:,3), 'Linewidth', 3)
semilogy(state_sol_no_donation_sample.A(:,1),...
    state_sol_no_donation_sample.A(:,4), 'Linewidth', 3)
semilogy(state_sol_no_donation_sample.A(:,1),...
    state_sol_no_donation_sample.A(:,5), 'Linewidth', 3)
semilogy(state_sol_no_donation_sample.A(:,1),...
    state_sol_no_donation_sample.A(:,6), 'Linewidth', 3)
semilogy(state_sol_no_donation_sample.A(:,1),...
    state_sol_no_donation_sample.A(:,7), 'Linewidth', 3)
ylim([10^0, 10^6])
xlim([0, 360])
xticks([0, 60, 120, 180, 240, 300, 360])
legend('$S_A$','$E_A$','$I_A$','$R_A$','$V_A$','$D_A$','Interpreter','latex')
xlabel('$Days$','Interpreter','latex')
ylabel('$Population$','Interpreter','latex')
set(gca,'FontSize',20);

figure;
semilogy(state_sol_no_donation_sample.A(:,1),...
    state_sol_no_donation_sample.B(:,2), 'Linewidth', 3)
hold on
semilogy(state_sol_no_donation_sample.A(:,1),...
    state_sol_no_donation_sample.B(:,3), 'Linewidth', 3)
semilogy(state_sol_no_donation_sample.A(:,1),...
    state_sol_no_donation_sample.B(:,4), 'Linewidth', 3)
semilogy(state_sol_no_donation_sample.A(:,1),...
    state_sol_no_donation_sample.B(:,5), 'Linewidth', 3)
semilogy(state_sol_no_donation_sample.A(:,1),...
    state_sol_no_donation_sample.B(:,6), 'Linewidth', 3)
semilogy(state_sol_no_donation_sample.A(:,1),...
    state_sol_no_donation_sample.B(:,7), 'Linewidth', 3)
ylim([10^0, 10^6])
xlim([0, 360])
xticks([0, 60, 120, 180, 240, 300, 360])
legend('$S_B$','$E_B$','$I_B$','$R_B$','$V_B$','$D_B$','Interpreter','latex')
xlabel('$Days$','Interpreter','latex')
ylabel('$Population$','Interpreter','latex')
set(gca,'FontSize',20);

% optimal donation
figure;
semilogy(state_sol_optimal_sample.A(:,1),...
    state_sol_optimal_sample.A(:,2), 'Linewidth', 3)
hold on
semilogy(state_sol_optimal_sample.A(:,1),...
    state_sol_optimal_sample.A(:,3), 'Linewidth', 3)
semilogy(state_sol_optimal_sample.A(:,1),...
    state_sol_optimal_sample.A(:,4), 'Linewidth', 3)
semilogy(state_sol_optimal_sample.A(:,1),...
    state_sol_optimal_sample.A(:,5), 'Linewidth', 3)
semilogy(state_sol_optimal_sample.A(:,1),...
    state_sol_optimal_sample.A(:,6), 'Linewidth', 3)
semilogy(state_sol_optimal_sample.A(:,1),...
    state_sol_optimal_sample.A(:,7), 'Linewidth', 3)
ylim([10^0, 10^6])
xlim([0, 360])
xticks([0, 60, 120, 180, 240, 300, 360])
legend('$S_A$','$E_A$','$I_A$','$R_A$','$V_A$','$D_A$','Interpreter','latex')
xlabel('$Days$','Interpreter','latex')
ylabel('$Population$','Interpreter','latex')
set(gca,'FontSize',20);

figure;
semilogy(state_sol_optimal_sample.A(:,1),...
    state_sol_optimal_sample.B(:,2), 'Linewidth', 3)
hold on
semilogy(state_sol_optimal_sample.A(:,1),...
    state_sol_optimal_sample.B(:,3), 'Linewidth', 3)
semilogy(state_sol_optimal_sample.A(:,1),...
    state_sol_optimal_sample.B(:,4), 'Linewidth', 3)
semilogy(state_sol_optimal_sample.A(:,1),...
    state_sol_optimal_sample.B(:,5), 'Linewidth', 3)
semilogy(state_sol_optimal_sample.A(:,1),...
    state_sol_optimal_sample.B(:,6), 'Linewidth', 3)
semilogy(state_sol_optimal_sample.A(:,1),...
    state_sol_optimal_sample.B(:,7), 'Linewidth', 3)
ylim([10^0, 10^6])
xlim([0, 360])
xticks([0, 60, 120, 180, 240, 300, 360])
legend('$S_B$','$E_B$','$I_B$','$R_B$','$V_B$','$D_B$','Interpreter','latex')
xlabel('$Days$','Interpreter','latex')
ylabel('$Population$','Interpreter','latex')
set(gca,'FontSize',20);
%}
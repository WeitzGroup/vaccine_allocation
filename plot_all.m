%% Parameters
% model parameters
parsM.Tinc = 4; % length of incubation period (days)
parsM.Tinf = 6; % duration patient is infectious (days)
parsM.etaI = 0.1; % *true* transmission effectivenesss, note transmission rate beta ~ etaI*contactRate
parsM.mu = 1e-2; % case fatality ratio
parsM.c_baseline = 5; % baseline contact rate (/day)
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

%% reduction in fatalities for different mu (Fig. 2 a/b/c)
mu_vec = 0:0.01:1;
parsM.kappa = 10^(-6);
parsS.vaccination_rate_baseline = 1*lambda;

for i =1:length(mu_vec)
    parsS.VA = parsM.total_vaccines * (1-(mu_vec(i)));
    parsS.VB = parsM.total_vaccines * (mu_vec(i));
    state_sol_test= state_solver(parsM, parsT, parsS,initial_state);
    deaths_A(i) = state_sol_test.A(end,end);
    deaths_B(i) = state_sol_test.B(end,end);
end

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
title('$\kappa = 10^{-6}$')



%% graphs for optimal strats
% heatmap showing optimal mu for different kappa and vaccination rates (Fig. 1)
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

figure; % draw contour lines for heatmap of mu = 33 perc and 0
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



%% actual fatalities in no-share, optimal and hybrid
% no-share Fig. S1 (a)
plot_heatmap(xvals, yvals, fatalities_all.baseline.A, 0, 8*10^4);
plot_heatmap(xvals, yvals, fatalities_all.baseline.B, 0, 8*10^4);

% optimal Fig. S1 (b)
plot_heatmap(xvals, yvals, fatalities_all.optimal.A, 0, 8*10^4);
plot_heatmap(xvals, yvals, fatalities_all.optimal.B, 0, 8*10^4);

% hybrid Fig. S1 (c)
plot_heatmap(xvals, yvals, fatality_hybrid.A, 0, 8*10^4);
plot_heatmap(xvals, yvals, fatality_hybrid.B, 0, 8*10^4);


%% compare hybrid to optimal Fig. S2
plot_heatmap(xvals, yvals, 100*fatality_reduction_hybrid_compare_optimal.A, 0, 100);
plot_heatmap(xvals, yvals, -100*fatality_reduction_hybrid_compare_optimal.B, 0, 100);

%% graph of optimal reduction and hybrid reduction for a and b for vac rate of 50, 75, 100% per year
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

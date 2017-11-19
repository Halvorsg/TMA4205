function plotting()
PCG_table = {load('flops_table_PCG.mat'),load('res_table_PCG.mat'),load('time_table_PCG.mat')};

V_cycle_table = {load('flops_table_V_cycle.mat'),load('res_table_V_cycle.mat'),load('time_table_V_cycle.mat')};

CG_table = {load('flops_table_CG.mat'),load('res_table_CG.mat'),load('time_table_CG.mat')};

h = figure(1);
hold on
plot(PCG_table{3}.time,log(PCG_table{2}.norm1_r))
plot(CG_table{3}.time,log(CG_table{2}.norm_r))
plot(V_cycle_table{3}.time,log(V_cycle_table{2}.norm1_r))
legend('PCG','CG','V-Cycle')
xlabel('Time(s)')
ylabel('ln(||r||_2)')
str = sprintf('Convergence of the three methods with \x03bb = %i',1000);
title(str)
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
xlim([0,12])
ylim([-8,13])
saveTightFigure(h,'Time_Convergence_1000')
hold off
h = figure(2);
hold on
plot(PCG_table{1}.FLOPS,log(PCG_table{2}.norm1_r))
plot(CG_table{1}.FLOPS,log(CG_table{2}.norm_r))
plot(V_cycle_table{1}.FLOPS,log(V_cycle_table{2}.norm1_r))
legend('PCG','CG','V-Cycle')
xlabel('Flops')
ylabel('ln(||r||_2)')
str = sprintf('Convergence of the three methods with \x03bb = %i',1000);
title(str)
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
saveTightFigure(h,'Flops_Convergence_1000')
hold off

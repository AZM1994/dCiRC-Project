%%% Compare the experimental and simulated phase of entrainment
%%% Written by Zheming An

clear,close all

organism_names  = ["JW168", "JW224", "JW60",  "JW260", "D117",...
                   "D119",  "D116",  "JW18",  "JW172", "P4483",...
                   "JW220", "JW200", "JW162", "JW176", "JW238"...
                   "JW261", "JW180", "JW169", "JW228", "P4463",...
                   "JW24",  "JW59",  "JW161", "P4469", "JW22"];

table_header = ["\tau_{DD}", "cz", "ke", "s", "a",...
"Active Compression", "Theoretical Compression", "Compression Ratio", "Active Expansion", "Theoretical Expansion",...
"Expansion Ratio","Active Integral Effect", "Theoretical Integral Effect", "Integral Effect Ratio","Initial Phase",...
"Simulated POE", "Simulated ROE", "Experimental POE"];

optimized_x_POE = [
                21.0	0.6477	1.5737	1.40	0.05	-0.27	18.00	48	-0.67	-0.80
                21.0	0.2436	0.7246	0.00	0.05	-0.24	18.00	72	1.07	1.10
                21.0	0.3696	1.5247	0.05	0.10	-0.0321	18.00	48	-0.33	-0.70
                21.2	0.3548	1.9021	0.00	0.10	-0.22	18.00	48	-0.06	0.60
                21.3	1.4254	1.7593	1.05	0.05	-0.17	18.00	48	-0.19	0.30
                21.3	1.0021	1.9621	0.00	0.05	-0.15	18.00	72	1.49	1.00
                21.5	0.2958	0.6504	0.05	0.05	0.02	18.00	72	1.03	0.90
                21.7	0.2658	0.4019	0.20	0.10	0.06	18.00	96	2.08	2.20
                21.7	0.3373	0.3216	1.00	0.05	-0.28	18.00	96	0.96	1.10
                21.7	0.1559	0.1701	0.85	0.60	-0.03	18.00	72	0.95	0.90
                21.8	0.1536	0.1060	1.40	0.40	-0.06	18.00	96	1.22	1.80
                22.0	0.5679	0.5591	0.30	0.10	-0.15	18.00	72	2.10	1.80
                22.2	0.9090	1.0826	0.40	0.25	-0.07	18.00	96	1.17	1.20
                22.2	0.3997	1.1089	0.10	0.25	0.08	18.00	96	1.76	2.10
                22.2	0.3784	1.9092	0.00	0.57	0.01	18.00	96	0.25	0.30
                22.2	0.4937	1.9138	0.05	0.15	0.41	18.00	96	0.40	0.30
                22.3	0.3879	1.0788	0.20	0.40	-0.20	18.00	96	1.28	1.70
                22.5	0.3429	0.7786	0.15	0.20	0.00	18.00	48	0.98	0.40
                22.7	0.6499	1.9883	0.00	0.10	0.33	18.00	72	1.05	1.30
                22.7	0.3258	0.7034	0.26	0.35	0.27	18.00	96	1.10	1.20
                23.2	0.2915	0.9880	0.20	0.10	-0.02	18.00	48	0.81	0.80
                23.5	0.5339	0.6208	0.60	0.30	0.17	18.00	72	1.02	0.90
                23.8	0.7301	1.2864	0.00	0.05	0.24	18.00	48	1.77	1.50
                23.8	0.4537	0.8656	0.60	0.20	0.28	18.00	48	0.21	0.20
                22.0	0.3636	1.1602	1.00	0.14	-0.99	18.00	96	0.38	0.50];

tau_I           = optimized_x_POE(:,1);
c_z             = optimized_x_POE(:,2);
k_e             = optimized_x_POE(:,3);
s               = optimized_x_POE(:,4);
a               = optimized_x_POE(:,5);

simulated_POE   = optimized_x_POE(:,9);
experimental_POE= optimized_x_POE(:,10);

cost = sqrt(sum((simulated_POE - experimental_POE).^2))/length(experimental_POE);

% Plot the experimental and simulated POE
figure(1);
scatter(tau_I, simulated_POE, 600,'filled', 'MarkerFaceAlpha', 1, 'LineWidth',2)
hold on
scatter(tau_I, experimental_POE, 600,'filled', 'd', 'MarkerFaceAlpha', 1, 'LineWidth',2)
P_sim = polyfit(tau_I,simulated_POE,1);
yfit_sim = P_sim(1)*tau_I+P_sim(2);
plot(tau_I,yfit_sim,'Color',[0 0.4470 0.7410],'LineStyle','--','LineWidth',3);
P_exp = polyfit(tau_I,experimental_POE,1);
yfit_exp = P_exp(1)*tau_I+P_exp(2);
plot(tau_I,yfit_exp,'Color',[0.8500 0.3250 0.0980],'LineStyle','--','LineWidth',3);
x_num = [3350/5000 3500/5000];
y_num = [3.5/4.5 2.9/4.5];
annotation('textarrow',x_num,y_num,'String','y = 0.2467 * x - 4.5785', 'FontSize',32,'Color',[0 0.4470 0.7410]) % sim
x_num = [3150/5000 3300/5000];
y_num = [1.9/4.5 2.65/4.5];
annotation('textarrow',x_num,y_num,'String','y = 0.1984 * x - 3.4813', 'FontSize',32,'Color',[0.8500 0.3250 0.0980]) % exp
text(21.1, 2.4, "C.I. = $0.0615$",'Color','k','Interpreter','latex','FontSize',32)
text(tau_I+0.04, simulated_POE, organism_names,'Color','b')
text(tau_I+0.04, experimental_POE, organism_names,'Color','r')
xlabel('$\tau_{DD}$/hours','Interpreter','latex')
ylabel('POE/hours','Interpreter','latex')
[~, hobj,~, ~] = legend('Simulated POE ($\Phi_{sim}$)','Experimental POE ($\Phi_{exp}$)',...
    'Location','best','FontSize',24,'Interpreter','latex');
M = findobj(hobj,'type','patch');
set(M,'MarkerSize',25,'LineWidth',3);
set(gca,'FontSize',32)

hAx = gca;                    % create an axes
hAx.LineWidth = 3;            % set the axis linewidth for box/ticks

% Plot the differences between the experimental and simulated POE
figure(2);
pos2 = [0.1 0.15 0.6 0.8];
subplot('Position',pos2)
color_hah = (simulated_POE - experimental_POE);
mean_diff = mean(simulated_POE - experimental_POE);
colormap winter
scatter(tau_I, simulated_POE - experimental_POE,2000,color_hah, 'filled','MarkerFaceAlpha',1,'LineWidth',2)
hold on
plot([21,24],[mean_diff,mean_diff],'--r','LineWidth',3)
x_num_2 = [2550/5000 2400/5000];
y_num_2 = [1.8/4.5 2.5/4.5];
annotation('textarrow',x_num_2,y_num_2,'String','mean value: -0.031 hours', 'FontSize',28,'Color',[0.8500 0.3250 0.0980]) % exp
text(tau_I-0.06, (simulated_POE - experimental_POE), organism_names,'Color','w','FontWeight','bold','FontSize',12)
xlabel('$\tau_{DD}/hours$','Interpreter','latex')
ylabel('$\Phi_{sim}-\Phi_{exp}/hours$','Interpreter','latex')
set(gca,'FontSize',32)
set(gca,'LineWidth',3)
pos1 = [0.75 0.16 0.2 0.79];
subplot('Position',pos1)
h = boxplot(simulated_POE-experimental_POE,'Labels',{'Box plot of'},'Width',0.6);
xlabel('$\Phi_{sim}-\Phi_{exp}$','Interpreter','latex')
set(gca,'YTickLabel',[]);
set(h,{'linew'},{2})
ylim([-1 1])
h = colorbar;
set(h, 'ylim', [-1 1])
set(gca,'FontSize',32)
set(gca,'LineWidth',3)
caxis([-1 1])


figure(3);
pos3 = [0.1 0.15 0.3 0.8];
subplot('Position',pos3)
cz_with_frq7 = [0.0394;c_z];
simulated_POE_with_frq7 = [1.69;simulated_POE];
scatter(simulated_POE_with_frq7,cz_with_frq7,600,'filled','MarkerFaceAlpha',1,'LineWidth',2)
hold on
plot([-0.97,2.1],[0.2958,0.2958],'--r','LineWidth',2);
plot([-0.97,2.1],[0.5679,0.5679],'--r','LineWidth',2);

text(-0.9, 0.15, 'strong oscillator','Color','k','FontWeight','bold','FontSize',26)
text(-0.9, 1.05, 'weak oscillator','Color','k','FontWeight','bold','FontSize',26)
yticks([0 0.3 0.6 0.9 1.2 1.5])
yticklabels({'0','0.3','0.6','0.9','1.2','1.5'})
ylabel('$c_z$','Interpreter','latex')
xlabel('Simulated POE ($\Phi_{sim}$) /hours','Interpreter','latex')

set(gca,'FontSize',32)
set(gca,'LineWidth',3)
pos4 = [0.45 0.16 0.3 0.79];
subplot('Position',pos4)
h = boxplot(cz_with_frq7,'Labels',{' '},'Width',0.6,'orientation', 'vertical', 'symbol', '', 'whisker',1000);
xlabel('Box plot of $c_z$','Interpreter','latex')
set(gca,'YTickLabel',[]);
set(h,{'linew'},{2})
ylim([0 1.5])

colormap(jet)
h = colorbar;
set(h, 'ylim', [0 1.5], 'YTickLabel',{'0.000','0.296','0.568','1.500'}, ...
               'YTick', [0,0.2958,0.5679,1.5])

set(gca,'FontSize',32)
set(gca,'LineWidth',3)
caxis([0 1.5])

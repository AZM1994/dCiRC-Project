%%% Study the non-uniform response of ecotypes in cycling condition
%%% plot of TAUC and LAUC with the linear regression line
%%% plot LAUC and period mismatch
%%% Written by Zheming An

clear, close all

% list of genotype names
Genotype_table = ["JW168", "JW224", "JW60",  "JW260", "D117",...
                  "D119",  "D116",  "JW18",  "JW172", "P4483",...
                  "JW220", "JW200", "JW162", "JW176", "JW238",...
                  "JW261", "JW180", "JW169", "JW228", "P4463",...
                  "JW24",  "JW59",  "JW161", "P4469", "JW22"];

% load all best fitting parameters and calculated TAUC and LAUC values
CIRC_parameter_table = [
21.0	0.6477	1.5737	1.40	0.05	0.0159	0.2293	0.0693	-0.3177	-4.5874	0.0693	-0.3040	-4.3580	0.0698
21.0	0.2436	0.7246	0.00	0.05	0.1917	0.3834	0.5000	-3.8339	-7.6757	0.4995	-3.6470	-7.2923	0.5001
21.0	0.3696	1.5247	0.05	0.10	0.3601	0.7582	0.4750	-3.6018	-7.5887	0.4746	-3.2461	-6.8306	0.4752
21.2	0.3548	1.9021	0.00	0.10	0.3834	0.7668	0.5000	-3.8344	-7.6751	0.4996	-3.4555	-6.9084	0.5002
21.3	1.4254	1.7593	1.05	0.05	0.0252	0.2422	0.1041	-0.5042	-4.8456	0.1040	-0.4816	-4.6034	0.1046
21.3	1.0021	1.9621	0.00	0.05	0.1917	0.3834	0.5000	-3.8339	-7.6757	0.4995	-3.6470	-7.2923	0.5001
21.5	0.2958	0.6504	0.05	0.05	0.1801	0.3791	0.4750	-3.6013	-7.5892	0.4745	-3.4260	-7.2102	0.4752
21.7	0.2658	0.4019	0.20	0.10	0.2852	0.7129	0.4000	-2.8520	-7.1348	0.3997	-2.5710	-6.4219	0.4004
21.7	0.3373	0.3216	1.00	0.05	0.0272	0.2449	0.1111	-0.5442	-4.8980	0.1111	-0.5197	-4.6531	0.1117
21.7	0.1559	0.1701	0.85	0.60	0.4164	3.0355	0.1372	-0.6950	-5.0597	0.1374	-0.2798	-2.0242	0.1382
21.8	0.1536	0.1060	1.40	0.40	0.1271	1.8347	0.0693	-0.3179	-4.5872	0.0693	-0.1922	-2.7524	0.0698
22.0	0.5679	0.5591	0.30	0.10	0.2356	0.6731	0.3500	-2.3561	-6.7354	0.3498	-2.1244	-6.0623	0.3504
22.2	0.9090	1.0826	0.40	0.25	0.4731	1.5770	0.3000	-1.8934	-6.3109	0.3000	-1.4234	-4.7339	0.3007
22.2	0.3997	1.1089	0.10	0.25	0.8390	1.8646	0.4500	-3.3581	-7.4640	0.4499	-2.5228	-5.5994	0.4505
22.2	0.3784	1.9092	0.00	0.57	2.1844	4.3697	0.4999	-3.8385	-7.6699	0.5005	-1.6563	-3.3001	0.5019
22.2	0.4937	1.9138	0.05	0.15	0.5401	1.1372	0.4750	-3.6022	-7.5882	0.4747	-3.0663	-6.4510	0.4753
22.3	0.3879	1.0788	0.20	0.40	1.1404	2.8515	0.3999	-2.8541	-7.1325	0.4002	-1.7165	-4.2810	0.4009
22.5	0.3429	0.7786	0.15	0.20	0.6209	1.4611	0.4250	-3.1061	-7.3110	0.4248	-2.4890	-5.8499	0.4255
22.7	0.6499	1.9883	0.00	0.10	0.3834	0.7668	0.5000	-3.8344	-7.6751	0.4996	-3.4555	-6.9084	0.5002
22.7	0.3258	0.7034	0.26	0.35	0.8928	2.4133	0.3700	-2.5531	-6.8986	0.3701	-1.6633	-4.4854	0.3708
23.2	0.2915	0.9880	0.20	0.10	0.2852	0.7129	0.4000	-2.8520	-7.1348	0.3997	-2.5710	-6.4219	0.4004
23.5	0.5339	0.6208	0.60	0.30	0.3441	1.6655	0.2066	-1.1478	-5.5531	0.2067	-0.8062	-3.8876	0.2074
23.8	0.7301	1.2864	0.00	0.05	0.1917	0.3834	0.5000	-3.8339	-7.6757	0.4995	-3.6470	-7.2923	0.5001
23.8	0.4537	0.8656	0.60	0.20	0.2294	1.1103	0.2066	-1.1475	-5.5533	0.2066	-0.9210	-4.4430	0.2073
22.0	0.3636	1.1602	1.00	0.14	0.0762	0.6856	0.1111	-0.5443	-4.8979	0.1111	-0.4705	-4.2123	0.1117];

label_list = ["$|\tau_{DD} - T|$", "cz", "$k_e$", "s", "a",...
"Active Compression", "Theoretical Compression", "Compression Ratio","Active Expansion", "Theoretical Expansion",...
"Expansion Ratio","Light-exposed Area Under Curve", "Total Area Under Curve", "Integral Effect Ratio","Initial Phase",...
"Simulated POE", "Simulated ROE", "Experimental POE"];

% Plot TAUC and LAUC
figure(1);
range_011 = [8,12,14,18,21];
range_012 = [1,5,9,10,11,13,15,17,20,22,24,25];
range_013 = [2,3,4,6,7,16,19,23];
colors = jet(90);
i = 12;
j = 13;
x = CIRC_parameter_table(:,i);
y = CIRC_parameter_table(:,j);
corr(x,y)
scatter(CIRC_parameter_table(range_011,i),CIRC_parameter_table(range_011,j),1000,colors(23,:),'filled','MarkerFaceAlpha',0.8)
hold on
scatter(CIRC_parameter_table(range_012,i),CIRC_parameter_table(range_012,j),1000,colors(23,:),'filled','MarkerFaceAlpha',0.8)
scatter(CIRC_parameter_table(range_013,i),CIRC_parameter_table(range_013,j),1000,colors(23,:),'filled','MarkerFaceAlpha',0.8)
hold on
p = polyfit(x,y,1);
yfit = polyval(p,x);
plot(x,yfit,'r','LineWidth',3)

text(CIRC_parameter_table(:,i)-0.08,CIRC_parameter_table(:,j),Genotype_table,'FontSize',14)
xlabel(label_list(i),'Interpreter','latex')
ylabel(label_list(j),'Interpreter','latex')
set(gca,'FontSize',32)
set(gca,'LineWidth',3)
xlim([-4 0])
ylim([-8 0])

%% Plot LAUC and period mismatch
figure(2);
range_011 = [8,12,14,18,21];
range_012 = [1,5,9,10,11,13,15,17,20,22,24,25];
range_013 = [2,3,4,6,7,16,19,23];
colors = jet(90);
i = 12;
j = 1;
x = CIRC_parameter_table(:,i);
y = 24 - CIRC_parameter_table(:,j);
scatter(CIRC_parameter_table(range_011,i),CIRC_parameter_table(range_011,j) - 24,1000,colors(23,:),'filled','MarkerFaceAlpha',0.8)
hold on
scatter(CIRC_parameter_table(range_012,i),CIRC_parameter_table(range_012,j) - 24,1000,colors(23,:),'filled','MarkerFaceAlpha',0.8)
scatter(CIRC_parameter_table(range_013,i),CIRC_parameter_table(range_013,j) - 24,1000,colors(23,:),'filled','MarkerFaceAlpha',0.8)

yticks([-3 -2.5 -2 -1.5 -1 -0.5 0])
yticklabels({'3','2.5','2','1.5','1','0.5','0'})
text(CIRC_parameter_table(:,i)-0.08,CIRC_parameter_table(:,j)-24,Genotype_table,'FontSize',14)
xlabel(label_list(i),'Interpreter','latex')
ylabel(label_list(j),'Interpreter','latex')
set(gca,'FontSize',32)
set(gca,'LineWidth',3)
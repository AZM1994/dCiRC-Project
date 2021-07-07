%%% Plot the convex cost - s - a space to find the best fitting (s,a)
%%% Written by Zheming An
clear, close all
% range of (s,a) pairs
list_s = 0:0.05:2;
list_a = 0:0.05:2; 

%% I
weighted_cost_matrix_I = load('Weighted_cost_matrix_JW18_I.mat');
weighted_cost_matrix_I = weighted_cost_matrix_I.cost_matrix;
[S,A] = meshgrid(list_s, list_a);

figure(1);
s_I = surf(S,A,weighted_cost_matrix_I','FaceAlpha',0.8);
[~,loc_d_I] = min(weighted_cost_matrix_I(:));
[min_s_I,min_a_I] = ind2sub(size(weighted_cost_matrix_I),loc_d_I);
min_cost_I = weighted_cost_matrix_I(min_s_I,min_a_I);
hold on
scatter3(list_s(min_s_I), list_a(min_a_I), min_cost_I, 800, 'r', 'filled')
xlabel('s')
ylabel('a')
zlabel('Weighted Cost JW18 I')
title_string_I = sprintf('Weighted Cost JW18 I \n minimum cost = %0.8f, at (s, a) = (%0.2f, %0.2f)',... 
                          min_cost_I, list_s(min_s_I), list_a(min_a_I));
title(title_string_I)
set(gca,'FontSize',32)
set(gca,'LineWidth',3)



%% II
weighted_cost_matrix_II = load('Weighted_cost_matrix_JW18_II.mat');
weighted_cost_matrix_II = weighted_cost_matrix_II.cost_matrix;

figure(2);
s_II = surf(S,A,weighted_cost_matrix_II','FaceAlpha',0.8);
[~,loc_d_II] = min(weighted_cost_matrix_II(:));
[min_s_II,min_a_II] = ind2sub(size(weighted_cost_matrix_II),loc_d_II);
min_cost_II = weighted_cost_matrix_II(min_s_II,min_a_II);
hold on
scatter3(list_s(min_s_II), list_a(min_a_II), min_cost_II, 800, 'r', 'filled')
xlabel('s')
ylabel('a')
zlabel('Weighted Cost JW18 II')
title_string_II = sprintf('Weighted Cost JW18 II \n minimum cost = %0.8f, at (s, a) = (%0.2f, %0.2f)',... 
                          min_cost_II, list_s(min_s_II), list_a(min_a_II));
title(title_string_II)
set(gca,'FontSize',32)
set(gca,'LineWidth',3)



%% III
weighted_cost_matrix_III = load('Weighted_cost_matrix_JW18_III.mat');
weighted_cost_matrix_III = weighted_cost_matrix_III.cost_matrix;

figure(3);
s_III = surf(S,A,weighted_cost_matrix_III','FaceAlpha',0.8);
[~,loc_d_III] = min(weighted_cost_matrix_III(:));
[min_s_III,min_a_III] = ind2sub(size(weighted_cost_matrix_III),loc_d_III);
min_cost_III = weighted_cost_matrix_III(min_s_III,min_a_III);
hold on
scatter3(list_s(min_s_III), list_a(min_a_III), min_cost_III, 800, 'r', 'filled')
xlabel('s')
ylabel('a')
zlabel('Weighted Cost JW18 III')
title_string_III = sprintf('Weighted Cost JW18 III \n minimum cost = %0.8f, at (s, a) = (%0.2f, %0.2f)',... 
                          min_cost_III, list_s(min_s_III), list_a(min_a_III));
title(title_string_III)
set(gca,'FontSize',32)
set(gca,'LineWidth',3)



%% Sum the costs for all three entries
weighted_cost_matrix_sum = weighted_cost_matrix_I + weighted_cost_matrix_II + weighted_cost_matrix_III;
figure(4);
list_s = 0:0.05:1.95;
list_a = 0:0.05:1.95; 
[S,A] = meshgrid(list_s, list_a);
s_sum = surf(S,A,weighted_cost_matrix_sum(1:40,1:40)','FaceAlpha',0.8);
[~,loc_d_sum] = min(weighted_cost_matrix_sum(:));
[min_s_sum,min_a_sum] = ind2sub(size(weighted_cost_matrix_sum),loc_d_sum);
min_cost_sum = weighted_cost_matrix_sum(min_s_sum,min_a_sum);
hold on
scatter3(list_s(min_s_sum), list_a(min_a_sum), min_cost_sum, 800, 'r', 'filled')
xlabel('s')
ylabel('a')
zlabel('Convex Cost')
ylim([0.05 1.95])
title_string_sum = sprintf('minimum cost = %0.2f, at (s, a) = (%0.2f, %0.2f)',... 
                          min_cost_sum, list_s(min_s_sum), list_a(min_a_sum));
title(title_string_sum)
set(gca,'FontSize',32)
set(gca,'LineWidth',3)
% s_sum.EdgeColor = 'none';

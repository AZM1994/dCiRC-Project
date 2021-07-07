%%% Plot the convex cost - s - a space to find the best fitting (s,a)
%%% Written by Zheming An

clear, close all
list_s = 0:0.05:2;
list_a = 0:0.05:2; 

[S,A] = meshgrid(list_s, list_a);

%% IV
weighted_cost_matrix_IV = load('Weighted_cost_matrix_frq7_IV.mat');
weighted_cost_matrix_IV = weighted_cost_matrix_IV.cost_matrix;

figure(4);
s_IV = surf(S,A,weighted_cost_matrix_IV','FaceAlpha',0.8);
[~,loc_d_IV] = min(weighted_cost_matrix_IV(:));
[min_s_IV,min_a_IV] = ind2sub(size(weighted_cost_matrix_IV),loc_d_IV);
min_cost_IV = weighted_cost_matrix_IV(min_s_IV,min_a_IV);
hold on
scatter3(list_s(min_s_IV), list_a(min_a_IV), min_cost_IV, 800, 'r', 'filled')
xlabel('s')
ylabel('a')
zlabel('Weighted Cost frq7 IV')
title_string_IV = sprintf('Weighted Cost frq7 IV \n minimum cost = %0.8f, at (s, a) = (%0.2f, %0.2f)',... 
                          min_cost_IV, list_s(min_s_IV), list_a(min_a_IV));
title(title_string_IV)
set(gca,'FontSize',32)

%% V
weighted_cost_matrix_V = load('Weighted_cost_matrix_frq7_V.mat');
weighted_cost_matrix_V = weighted_cost_matrix_V.cost_matrix;

figure(5);
s_V = surf(S,A,weighted_cost_matrix_V','FaceAlpha',0.8);
[~,loc_d_V] = min(weighted_cost_matrix_V(:));
[min_s_V,min_a_V] = ind2sub(size(weighted_cost_matrix_V),loc_d_V);
min_cost_V = weighted_cost_matrix_V(min_s_V,min_a_V);
hold on
scatter3(list_s(min_s_V), list_a(min_a_V), min_cost_V, 800, 'r', 'filled')
xlabel('s')
ylabel('a')
zlabel('Weighted Cost frq7 V')
title_string_V = sprintf('Weighted Cost frq7 V \n minimum cost = %0.8f, at (s, a) = (%0.2f, %0.2f)',... 
                          min_cost_V, list_s(min_s_V), list_a(min_a_V));
title(title_string_V)
set(gca,'FontSize',32)

%% VI
weighted_cost_matrix_VI = load('Weighted_cost_matrix_frq7_VI.mat');
weighted_cost_matrix_VI = weighted_cost_matrix_VI.cost_matrix;

figure(6);
s_VI = surf(S,A,weighted_cost_matrix_VI','FaceAlpha',0.8);
[~,loc_d_VI] = min(weighted_cost_matrix_VI(:));
[min_s_VI,min_a_VI] = ind2sub(size(weighted_cost_matrix_VI),loc_d_VI);
min_cost_VI = weighted_cost_matrix_VI(min_s_VI,min_a_VI);
hold on
scatter3(list_s(min_s_VI), list_a(min_a_VI), min_cost_VI, 800, 'r', 'filled')
xlabel('s')
ylabel('a')
zlabel('Weighted Cost frq7 VI')
title_string_VI = sprintf('Weighted Cost frq7 VI \n minimum cost = %0.8f, at (s, a) = (%0.2f, %0.2f)',... 
                          min_cost_VI, list_s(min_s_VI), list_a(min_a_VI));
title(title_string_VI)
set(gca,'FontSize',32)

%% 
% weighted_cost_matrix_sum = weighted_cost_matrix_I + weighted_cost_matrix_II + weighted_cost_matrix_III...
%                          + weighted_cost_matrix_IV + weighted_cost_matrix_V + weighted_cost_matrix_VI;
weighted_cost_matrix_sum = weighted_cost_matrix_IV + weighted_cost_matrix_V + weighted_cost_matrix_VI;
figure(7);
s_sum = surf(S,A,weighted_cost_matrix_sum','FaceAlpha',0.8);
[~,loc_d_sum] = min(weighted_cost_matrix_sum(:));
[min_s_sum,min_a_sum] = ind2sub(size(weighted_cost_matrix_sum),loc_d_sum);
min_cost_sum = weighted_cost_matrix_sum(min_s_sum,min_a_sum);
hold on
scatter3(list_s(min_s_sum), list_a(min_a_sum), min_cost_sum, 800, 'r', 'filled')
xticks(0:0.5:2)
xticklabels({'0','0.5','1','1.5','2'})
yticks(0:0.5:2)
yticklabels({'0','0.5','1','1.5','2'})
xlabel('s')
ylabel('a')
zlabel('Convex Cost')
% title_string_sum = sprintf('Weighted Cost frq7 sum \n minimum cost = %0.8f, at (s, a) = (%0.2f, %0.2f)',... 
%                           min_cost_sum, list_s(min_s_sum), list_a(min_a_sum));
title_string_sum = sprintf('minimum cost = %0.2f, at (s, a) = (%0.2f, %0.2f)',... 
                          min_cost_sum, list_s(min_s_sum), list_a(min_a_sum));
title(title_string_sum)
set(gca,'FontSize',32)
set(gca,'LineWidth',3)
% s_sum.EdgeColor = 'none';

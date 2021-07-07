%%% Plot the convex cost - s - a space to find the best fitting (s,a)
%%% Written by Zheming An

clear, close all
list_s = 0:0.05:2;
list_a = 0:0.05:2; 

%% I
weighted_cost_matrix_I = load('Weighted_cost_matrix_HG2.mat');
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
zlabel('Convex Cost')
title_string_sum = sprintf('minimum cost = %0.2f, at (s, a) = (%0.2f, %0.2f)',... 
                          min_cost_I, list_s(min_s_I), list_a(min_a_I));
title(title_string_sum)
set(gca,'FontSize',32)
set(gca,'LineWidth',3)
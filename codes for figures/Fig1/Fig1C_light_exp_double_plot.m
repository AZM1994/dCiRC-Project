%%% Generating Fig. 1 C for the illustrative example of the dCiRC model,
%%% reflecting both the parametric and non-parametric aspects of the
%%% entrainment under cycling conditions.
%%% Fig. 1 C: Double plot showing the trajectory of the phase differences between the internal
%%% clock and the zeitgeber clock.
%%% Written by Zheming An

clear, close all
%% load phase differences data
analytical_phase_difference = load('analytical_phase_difference.mat');
analytical_z_time = load('analytical_z_time.mat');
analytical_phase_difference = analytical_phase_difference.analytical_phase_difference;
analytical_z_time = analytical_z_time.analytical_z_time;

% parameters
T_interval       = 0.01; % unit in hour, time interval 
zeitgeber_period = 24; % unit in hour

% adjust the phase differences data for the double plot
analytical_phase_difference_adjusted = zeros(length(analytical_phase_difference),1);
analytical_phase_difference_adjusted(1:30) = analytical_phase_difference(1:30);
analytical_phase_difference_adjusted(31:65) = analytical_phase_difference(31:65) + zeitgeber_period;
analytical_phase_difference_adjusted(66:end) = analytical_phase_difference(66:end) + 2 * zeitgeber_period;

%% Double plot
figure(1);
stripe = [360, 480, 360, 480];
X1 = [6 12 18];
X2 = [30 36 42];
Y1 = [stripe;stripe;stripe];
newcolors = [255 255 0;128 128 128]/256;
colororder(newcolors)
area(0:6, zeros(7,1) + 1680, 'FaceColor',[128 128 128]/256,'LineStyle','none','FaceAlpha', 0.6)
hold on
area(18:30, zeros(13,1) + 1680, 'FaceColor',[128 128 128]/256,'LineStyle','none','FaceAlpha', 0.6)
area(42:48, zeros(7,1) + 1680, 'FaceColor',[128 128 128]/256,'LineStyle','none','FaceAlpha', 0.6)
area(X1, Y1, 'LineStyle','none','FaceAlpha', 0.6)
area(X2, Y1, 'LineStyle','none','FaceAlpha', 0.6)
scatter(flip(analytical_phase_difference_adjusted), analytical_z_time*T_interval+12, 60, 'filled', 'MarkerFaceColor', 'b');
plot([24 24],[0 1680],'--k','LineWidth',3,'HandleVisibility','off')
plot([2.96 2.96],[0 1194],'-.k','LineWidth',1,'HandleVisibility','off')
plot([26.96 26.96],[0 354],'-.k','LineWidth',1,'HandleVisibility','off')
plot([0 48],[192 192],'r','LineWidth',3)
plot([0 48],[360 360],'r','LineWidth',3)
plot([0 48],[504 504],'r','LineWidth',3)
plot([0 48],[840 840],'r','LineWidth',3)
plot([0 48],[1032 1032],'r','LineWidth',3)
plot([0 48],[1200 1200],'r','LineWidth',3)
plot([0 48],[1416 1416],'r','LineWidth',3)
text(3.5,50, '$\Psi$', 'Interpreter', 'latex', 'Fontsize',36, 'Color', 'k', 'FontWeight', 'bold')
text(27.5,50, '$\Psi$', 'Interpreter', 'latex', 'Fontsize',36, 'Color', 'k', 'FontWeight', 'bold')

xticks([0 6 12 18 24 30 36 42 48])
xticklabels({'0','6','12','18','24(0)','6','12','18','24'})
yticks(flip(1680-[0 264 480 648 840 1176 1320 1488 1680]))
yticklabels({'70','62','55','49','35','27','20','11','0'})
xlabel('External Time/hours','Interpreter','latex')
ylabel('days of simulation','Interpreter','latex')
set(gca,'fontsize',30);
set(gca,'LineWidth',3)
axis('square')
xlim([0 48])
ylim([0 1680])

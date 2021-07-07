%%% Convex optimization (PSO) algorithm to find the best fitting dCIRC parameters
%%% Written by Zheming An
clear;close all;
tic
global omega_DD T_interval Phase_Data days_P number_of_timesteps cost_range ...
    T_interval_stable Light_Profile alpha_1 alpha_2 alpha_3 alpha_4
%% Load normalized race tube data
HG_data   = load('linear-HG2-0.1.mat');
HG_data   = HG_data.data_saved;
Light_Profile = load('linear-Light_Profile.mat');
Light_Profile = Light_Profile.Light_Profile;
cost_range      = 1:139;
days_P          = HG_data(:,1);
Phase_Data      = HG_data(:,2);
number_of_timesteps = length(Phase_Data);



%% PSO Convex Optimization
for ii = 1 : number_of_timesteps - 1
	T_interval(ii) = HG_data(ii+1,1) - HG_data(ii,1);
end                            
T_interval_stable = max(T_interval);

alpha_1             = 0; % paramter for convex cost function
alpha_2             = 0; % paramter for convex cost function
alpha_3             = 0; % paramter for convex cost function
alpha_4             = 0; % paramter for convex cost function
tau_DD              = 24.72; % free running period for the selected genotype
omega_DD            = 2 * pi ./ tau_DD; % free running velocity for the selected genotype

% set optimization ranges for the dCIRC parameters
% the order is: ["cz", "elastic fator", "s", "a", "initial theta_I"]
% "cz" is the zeitgeber strength
% "s" is the shape factor of the CIRC
% "a" is the asymmetry factor of the CIRC
% "initial theta_I" is the initial phase angle of the internal clock
lower_bound         = [0.0286817541130158,2.79336508461019,0,1.99999721333515,5.45232999897579];
upper_bound         = [0.0286817541130158,2.79336508461019,0,1.99999721333515,5.45232999897579];
% lower_bound         = [0   0   0   0   0];
% upper_bound         = [3   3   2   2   2*pi]; % ["cz", "elastic", "s", "a", "initial theta_I", "initial omega_I"] 3   3   2   5   5
number_of_variables = 5; % number of optimation variables

% set initial condition for the optimation variables
% initial_swarm       = [0.988233909045979,0.000970852225892652,0.999991057580442,0.00102368402763053,5.75005674814211];
initial_swarm       = [0	0	 0	 1   0];  % Default initial conditions

options = optimoptions(@particleswarm,'PlotFcn','pswplotbestf','SwarmSize',100,'HybridFcn',@fmincon,...
          'FunctionTolerance', 10^-6, 'MaxStallIterations', 5, 'InitialSwarmMatrix',initial_swarm); % , 'Display','iter'
[optimized_x ,Fval,exitFlag,Output] = particleswarm(@CIRC_LD_Cost,number_of_variables,lower_bound,upper_bound,options); 

% display the optimized dCIRC parameters
format long
disp(optimized_x)
minimum_cost = CIRC_LD_Cost(optimized_x);

% Plot the best fitting curves
% Plot the clock running velocity during the entrainment
CIRC_LD_Cost_for_plot_CIRC_Curve(optimized_x)

% Plot the best fitting dCIRC curve
[L2_cost, convex_cost] = CIRC_LD_Cost_for_plot_best_fitting(optimized_x, minimum_cost);

toc 



%% Define the convex cost function for optimization
function convex_cost = CIRC_LD_Cost(x)

global omega_DD T_interval Phase_Data number_of_timesteps cost_range ...
    T_interval_stable Light_Profile alpha_1 alpha_2 alpha_3 alpha_4
% preallocate 
ALL_clocks          = zeros(number_of_timesteps,2); 
zeitgeber_clock     = zeros(number_of_timesteps,1); 
Normalized_theta_Z  = zeros(number_of_timesteps,1); 
Normalized_theta_I  = zeros(number_of_timesteps,1); 
CIRC                = zeros(number_of_timesteps,1); 
light               = zeros(number_of_timesteps,1); 
omega_I             = zeros(number_of_timesteps,1); 
domega_I            = zeros(number_of_timesteps,1); 

% Parameters
zeit_t          = 24; % period of zeitgeber clock
C_ZG            = x(1); % zeitgeber strength
elastic_factor  = x(2); % elastic factor
s               = x(3); % shape factor of the CIRC
a               = x(4); % asymmetry factor of the CIRC
initial_theta_I = x(5); % initial phase angle of the internal clock
initial_omega_I = omega_DD; % free running velocity for the selected genotype
CIRC_shift      = 0; % the time shift of the CIRC curve

% The desinated light profiles, given by center of the light explosure, 
% and radius of the light explosure
light_profiles       = [6,9]* 2 * pi / zeit_t;
% convert zeitgeber time to the phase angle, pi / zeit_t corresponds to 30 mins 

light_profiles_index = 1; % choose the light profile
epsilon              = light_profiles(light_profiles_index,1); % radius of the light explosure
light_middle         = light_profiles(light_profiles_index,2); % center of the light explosure

%% Simulate the entrainment process
theta_I_initial = initial_theta_I; 
omega_I_initial = initial_omega_I; 
initial_conditions = [theta_I_initial, omega_I_initial];
ALL_clocks(1,:) = initial_conditions;

% Runge-kutta forth order iterative method
for i = 1:number_of_timesteps - 1
    zeitgeber_clock(i) = (i-1)*2*pi/(zeit_t/T_interval_stable);
    k1= RK4_circadian_CIRC(ALL_clocks(i,:));
    k2_clocks = ALL_clocks(i,:) + 0.5*T_interval(i)*k1;

    k2= RK4_circadian_CIRC(k2_clocks);
    k3_clocks = ALL_clocks(i,:) + 0.5*T_interval(i)*k2;

    k3= RK4_circadian_CIRC(k3_clocks);
    k4_clocks = ALL_clocks(i,:) + T_interval(i)*k3;

    k4= RK4_circadian_CIRC(k4_clocks);
    ALL_clocks(i+1,:) =  ALL_clocks(i,:) + T_interval(i)*(k1 + 2*k2 + 2*k3 + k4) / 6.0; 
end

% Calculate the sinusoidal of phase angle of the internal clock 
theta_I_signal = sin(ALL_clocks(:,1));

% COST
fitting_sincurve_P = theta_I_signal; 

% define the convex cost function
convex_cost = sqrt(sum((fitting_sincurve_P(cost_range) - Phase_Data(cost_range)).^2)) ...
                   + alpha_1 * (C_ZG - 0.3)^2 ...
                   + alpha_2 * (elastic_factor - 1)^2 ...
                   + alpha_3 * s^2 ...
                   + alpha_4 * (a-1)^2;
               
function [ derivative ] = RK4_circadian_CIRC(y)
    theta_I = y(1,1);
    omega_I(i) = y(1,2);
    Normalized_theta_Z(i) = mod(zeitgeber_clock(i),2*pi);
    Normalized_theta_I(i) = mod(theta_I,2*pi);
    % Light function
%     if abs(Normalized_theta_Z(i) - light_middle) < epsilon 
%             light(i) = c_z;
%     else
%         light(i) = 0;
%     end
    light(i) = C_ZG * Light_Profile(i);
    % CIRC function
    if Normalized_theta_I(i) >= CIRC_shift && Normalized_theta_I(i) <= pi + CIRC_shift
        CIRC(i) = sin(Normalized_theta_I(i) - CIRC_shift) + s * sin(2 * (Normalized_theta_I(i) - CIRC_shift));
        if CIRC(i) < 0
            CIRC(i) = 0;
        end
        if a < 1
            CIRC(i) = CIRC(i) * a;
        end
    else
        CIRC(i) = -sin(2 * pi - (Normalized_theta_I(i) - CIRC_shift)) -...
               s * sin(2 * pi - (2 * (Normalized_theta_I(i) - CIRC_shift)));
        if CIRC(i) > 0
            CIRC(i) = 0;
        end
        if a > 1
            CIRC(i) = CIRC(i) / a;
        end
    end

    if s < 1
        max_1 = -0.5669 * s^3 + 1.1431 * s^2 + 0.1703 * s + 0.9963;
        CIRC(i) = CIRC(i) / max_1;
    else
        max_2 = 0.0029 * s^2 + 0.9738 * s + 0.7783;
        CIRC(i) = CIRC(i) / max_2;
    end
    % dCIRC Model ODEs
    dtheta_I = omega_I(i); 
    domega_I(i) = light(i) * CIRC(i) + elastic_factor * (omega_DD - omega_I(i)); 
    derivative = [dtheta_I, domega_I(i)]; 
        
end
end

%% Plot the best fitting curves
function [L2_cost, convex_cost] = CIRC_LD_Cost_for_plot_best_fitting(optimized_x, minimum_cost)

global omega_DD T_interval Phase_Data days_P number_of_timesteps cost_range ...
    T_interval_stable Light_Profile alpha_1 alpha_2 alpha_3 alpha_4
% preallocate 
ALL_clocks                  = zeros(number_of_timesteps,2); 
zeitgeber_clock             = zeros(number_of_timesteps,1); 
Normalized_theta_Z          = zeros(number_of_timesteps,1); 
Normalized_theta_I          = zeros(number_of_timesteps,1); 
CIRC                        = zeros(number_of_timesteps,1); 
light                       = zeros(number_of_timesteps,1); 
omega_I                     = zeros(number_of_timesteps,1); 
domega_I                    = zeros(number_of_timesteps,1); 

% Parameters
zeit_t          = 24; 
amplitude       = 1; 
C_ZG            = optimized_x(1); 
elastic_factor  = optimized_x(2);
s               = optimized_x(3);
a               = optimized_x(4);
initial_theta_I = optimized_x(5); 
initial_omega_I = omega_DD; 
CIRC_shift      = 0;
%% 
light_profiles       = [6,9]* 2 * pi / zeit_t; % dawn dusk allday 
light_profiles_index = 1; 
epsilon              = light_profiles(light_profiles_index,1); % Window of entrainment pi/zeit_t corresponds to 30 mins 
light_middle         = light_profiles(light_profiles_index,2); 

%% Simulations
theta_I_initial = initial_theta_I; 
omega_I_initial = initial_omega_I; 
initial_conditions = [theta_I_initial, omega_I_initial];
ALL_clocks(1,:) = initial_conditions;

for i = 1:number_of_timesteps - 1
    zeitgeber_clock(i) = (i-1)*2*pi/(zeit_t/T_interval_stable);
    k1= RK4_circadian_CIRC(ALL_clocks(i,:));
    k2_clocks = ALL_clocks(i,:) + 0.5*T_interval(i)*k1;

    k2= RK4_circadian_CIRC(k2_clocks);
    k3_clocks = ALL_clocks(i,:) + 0.5*T_interval(i)*k2;

    k3= RK4_circadian_CIRC(k3_clocks);
    k4_clocks = ALL_clocks(i,:) + T_interval(i)*k3;

    k4= RK4_circadian_CIRC(k4_clocks);
    ALL_clocks(i+1,:) =  ALL_clocks(i,:) + T_interval(i)*(k1 + 2*k2 + 2*k3 + k4) / 6.0;             
end
    
    theta_I_signal = sin(ALL_clocks(:,1));

    fitting_sincurve_P = amplitude .* theta_I_signal;
    convex_cost = sqrt(sum((fitting_sincurve_P(cost_range) - Phase_Data(cost_range)).^2)) ...
                   + alpha_1 * (C_ZG - 0.3)^2 ...
                   + alpha_2 * (elastic_factor - 1)^2 ...
                   + alpha_3 * s^2 ...
                   + alpha_4 * (a-1)^2;
	L2_cost = sqrt(sum((fitting_sincurve_P(cost_range) - Phase_Data(cost_range)).^2));
    
    % find the intial angle the internal clock
    [Af, Bf] = min(abs(Phase_Data(cost_range(1:80))));
    actural_initial = fitting_sincurve_P(cost_range(Bf));
    assignin('base','Af',Af)
    assignin('base','Bf',Bf)
    assignin('base','actural_initial',actural_initial)


    
    %% plot the best fitting curves
    figure(5);
        plot(days_P(1:end-1), zeros(1, number_of_timesteps - 1) + 2*pi/zeit_t, '--k', 'LineWidth',0.1);
        xticks(3:12:147)
        xticklabels({'ZT0','ZT12','ZT0','ZT12','ZT0','ZT12','ZT0','ZT12','ZT0','ZT12','ZT0','ZT12','ZT0'})
        yticks(0.2:0.1:0.3)
        yticklabels({' ',' ',' '})
        set(gca,'fontsize',32);
        ax1 = gca; % current axes
        ax1.XAxisLocation = 'top';
        xlim([days_P(1) days_P(end)])
        ylim([min(nonzeros(omega_I)) - 0.03 max(omega_I)])
        ax1_pos = ax1.Position; % position of first axes
        ax2 = axes('Position',[0.130000000000000,0.110000000000000,0.775000000000000,0.816904732046092],...
        'XAxisLocation','bottom',...
        'YAxisLocation','right',...
        'Color','none');
        plot(days_P, Phase_Data, 'k', 'LineWidth',6)
        hold on
        plot(days_P(cost_range), fitting_sincurve_P(cost_range), 'r', 'LineWidth',6)
        area(days_P,zeros(number_of_timesteps,1) + 1, 'FaceColor', 'k', 'FaceAlpha', 1, 'LineStyle','none')
        area(days_P,zeros(number_of_timesteps,1) - 1, 'FaceColor', 'k', 'FaceAlpha', 1, 'LineStyle','none')
        area(days_P, light/C_ZG,...
            'FaceColor', 'y', 'FaceAlpha', 1, 'LineStyle','none')
        area(days_P, -light/C_ZG,...
            'FaceColor', 'y', 'FaceAlpha', 1, 'LineStyle','none')
        alpha(0.3)
        set(gca,'FontSize',36)
        set(gcf,'Position',[550 80 700 450]);
        disp([L2_cost, convex_cost, minimum_cost])
        xlabel('time/hours','Interpreter','latex')
        ylabel('rhythmic signal','Interpreter','latex')
        xlim([days_P(1) days_P(end)])
        ylim([-1 1])
        xticks(0:12:144)
        xticklabels({'0','12','24','36','48','60','72','84','96','108','120','132','144'})
        title('')
        legend('Detrend Data', 'Best Fitting', 'AutoUpdate','off')
        

    % Plot the clock running velocity during the entrainment
    figure(2);
        plot(days_P(1:end-1), zeros(1, number_of_timesteps - 1) + 2*pi/zeit_t, '--k', 'LineWidth',0.1);
        xticks(3:12:147)
        xticklabels({'ZT0','ZT12','ZT0','ZT12','ZT0','ZT12','ZT0','ZT12','ZT0','ZT12','ZT0','ZT12','ZT0'})
        yticks(0.2:0.1:0.3)
        yticklabels({' ',' ',' '})
        set(gca,'fontsize',32);
        ax1 = gca; % current axes
        ax1.XAxisLocation = 'top';
        xlim([days_P(1) days_P(end)])
        ylim([min(nonzeros(omega_I)) - 0.03 max(omega_I)])
        ax2 = axes('Position',[0.130000000000000,0.110000000000000,0.775000000000000,0.830714311766626],...
        'XAxisLocation','bottom',...
        'YAxisLocation','right',...
        'Color','none');
        plot(days_P(1:end-1), omega_I(1:end-1),'r', 'LineWidth',6);
        legend('angular velocity of internal clock ($\omega_I$)', 'Location', 'best', 'AutoUpdate','off','Interpreter','latex')
        hold on
        plot(days_P(1:end-1), zeros(1, number_of_timesteps - 1) + 2*pi/zeit_t, '--k', ...
             days_P(1:end-1), zeros(1, number_of_timesteps - 1) + omega_DD, '--b', 'LineWidth',6);
        area(days_P,zeros(number_of_timesteps,1) + 1, 'FaceColor', 'k', 'FaceAlpha', 1, 'LineStyle','none')
        area(days_P, light/C_ZG,'FaceColor', 'y', 'FaceAlpha', 1, 'LineStyle','none')
        alpha(0.3)
        xlabel('time/hours','Interpreter','latex')
        ylabel('$\omega_I/rad \cdot hour^{-1}$','Interpreter','latex')
        xticks(0:12:144)
        xticklabels({'0','12','24','36','48','60','72','84','96','108','120','132','144'})
        xlim([days_P(1) days_P(end)])
        ylim([min(nonzeros(omega_I)) - 0.03 max(omega_I) + 0.03])
        text(days_P(end)-5, 2*pi/zeit_t+0.002, '$\omega_T=\frac{2\pi}{T}$','Interpreter','latex','FontSize',42)
        text(days_P(end)-5, omega_DD-0.002, '$\omega_{I0}=\frac{2\pi}{\tau_{I0}}$','Interpreter','latex','FontSize',42)
        set(gca,'fontsize',36);
        set(gcf,'Position',[0 480 700 350]);
                
function [ derivative ] = RK4_circadian_CIRC(y)
    theta_I = y(1,1);
    omega_I(i) = y(1,2);
    Normalized_theta_Z(i) = mod(zeitgeber_clock(i),2*pi);
    Normalized_theta_I(i) = mod(theta_I,2*pi);
    
    % Light function
%     if abs(Normalized_theta_Z(i) - light_middle) < epsilon 
%             light(i) = c_z;
%     else
%         light(i) = 0;
%     end
    light(i) = C_ZG * Light_Profile(i);
    
    % CIRC function
    if Normalized_theta_I(i) >= CIRC_shift && Normalized_theta_I(i) <= pi + CIRC_shift
        CIRC(i) = sin(Normalized_theta_I(i) - CIRC_shift) + s * sin(2 * (Normalized_theta_I(i) - CIRC_shift));
        if CIRC(i) < 0
            CIRC(i) = 0;
        end
        if a < 1
            CIRC(i) = CIRC(i) * a;
        end
    else
        CIRC(i) = -sin(2 * pi - (Normalized_theta_I(i) - CIRC_shift)) -...
               s * sin(2 * pi - (2 * (Normalized_theta_I(i) - CIRC_shift)));
        if CIRC(i) > 0
            CIRC(i) = 0;
        end
        if a > 1
            CIRC(i) = CIRC(i) / a;
        end
    end

    if s < 1
        max_1 = -0.5669 * s^3 + 1.1431 * s^2 + 0.1703 * s + 0.9963;
        CIRC(i) = CIRC(i) / max_1;
    else
        max_2 = 0.0029 * s^2 + 0.9738 * s + 0.7783;
        CIRC(i) = CIRC(i) / max_2;
    end
    
    % dCRIC Model ODEs
        dtheta_I = omega_I(i); 
        domega_I(i) = light(i) * CIRC(i) + elastic_factor * (omega_DD - omega_I(i)); 
        derivative = [dtheta_I, domega_I(i)]; 
end
end

%% Plot the best fitting dCIRC curve
function CIRC_LD_Cost_for_plot_CIRC_Curve(optimized_x)
T_interval      = 0.01;
T_initial       = 0;
T_final         = 24;
T_final_plot    = T_final - T_interval;
% preallocate 
number_of_timesteps         = length(T_initial:T_interval:(T_final - T_interval));
ALL_clocks                  = zeros(number_of_timesteps,2);
zeitgeber_clock             = zeros(number_of_timesteps,1);
Normalized_theta_Z          = zeros(number_of_timesteps,1);
Normalized_theta_I          = zeros(number_of_timesteps,1);
CIRC                        = zeros(number_of_timesteps,1);
light                       = zeros(number_of_timesteps,1);
omega_I                     = zeros(number_of_timesteps,1);
domega_I                    = zeros(number_of_timesteps,1);

% Parameters
zeit_t          = 24; 
tau             = 24;
omega           = 2 * pi ./ tau; 
C_ZG            = 0.001; 
elastic_factor  = 0.001;
s               = optimized_x(3); 
a               = optimized_x(4); 
initial_theta_I = 0;
initial_omega_I = 2 * pi ./ zeit_t; 
CIRC_shift      = 0; 
%% 
light_profiles       = [6,9] * 2 * pi / zeit_t;
light_profiles_index = 1; 
epsilon              = light_profiles(light_profiles_index,1); % Window of entrainment pi/zeit_t corresponds to 30 mins 
light_middle         = light_profiles(light_profiles_index,2); 

%% Simulations
theta_I_initial = initial_theta_I; 
omega_I_initial = initial_omega_I; 
initial_conditions = [theta_I_initial, omega_I_initial];
ALL_clocks(1,:) = initial_conditions;

for i = 1:number_of_timesteps - 1
    zeitgeber_clock(i) = (i-1)*2*pi/(zeit_t/T_interval);
    k1= RK4_circadian_CIRC(ALL_clocks(i,:));
    k2_clocks = ALL_clocks(i,:) + 0.5*T_interval*k1;

    k2= RK4_circadian_CIRC(k2_clocks);
    k3_clocks = ALL_clocks(i,:) + 0.5*T_interval*k2;

    k3= RK4_circadian_CIRC(k3_clocks);
    k4_clocks = ALL_clocks(i,:) + T_interval*k3;

    k4= RK4_circadian_CIRC(k4_clocks);
    ALL_clocks(i+1,:) =  ALL_clocks(i,:) + T_interval*(k1 + 2*k2 + 2*k3 + k4) / 6.0;             
end

%% plot the best fitting dCIRC curve
    figure(3); 
        plot(T_initial:T_interval:T_final_plot, CIRC, 'k', 'LineWidth',8)
    hold on
    area(T_initial:T_interval:T_final_plot,zeros(number_of_timesteps,1) + 1, 'FaceColor', 'k', 'FaceAlpha', 1, 'LineStyle','none')
    area(T_initial:T_interval:T_final_plot,zeros(number_of_timesteps,1) - 1, 'FaceColor', 'k', 'FaceAlpha', 1, 'LineStyle','none')
    area(T_initial:T_interval:T_final_plot, light/C_ZG,'FaceColor', 'y', 'FaceAlpha', 1, 'LineStyle','none')
    area(T_initial:T_interval:T_final_plot, -light/C_ZG,'FaceColor', 'y', 'FaceAlpha', 1, 'LineStyle','none')
    alpha(0.3)
    xlim([0 24])
    ylim([-1.005 1.005])
    dim = [0.61 .52 .3 .3];
    str = {['\color{black} \tau=',num2str(tau,'%.1f'),', s=', num2str(s,'%.2f'), ', a=', num2str(a,'%.2f')]};
    annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',35,'EdgeColor','none');
    set(gca,'fontsize',42);
    set(gcf,'Position',[0 80 720 350]);
    xticks(0:6:30)
    xticklabels({'CT0','CT6','CT12','CT18','CT24','CT30'})
    xlabel('internal time/hours','Interpreter','latex')
    ylabel('capacity to change $\tau_I$','Interpreter','latex')
    text(0.5,0.1,'compression','Color','r','FontSize',32,'Rotation',90);
    text(0.5,-0.8,'expansion','Color','b','FontSize',32,'Rotation',90);
    grid on
    grid minor
    
function [ derivative ] = RK4_circadian_CIRC(y)
    theta_I = y(1,1);
    omega_I(i) = y(1,2);
    Normalized_theta_Z(i) = mod(zeitgeber_clock(i),2*pi);
    Normalized_theta_I(i) = mod(theta_I,2*pi);
    
    % Light function
    if abs(Normalized_theta_Z(i) - light_middle) < epsilon 
            light(i) = C_ZG;
    else
        light(i) = 0;
    end
    
    % CIRC function
    if Normalized_theta_I(i) >= CIRC_shift && Normalized_theta_I(i) <= pi + CIRC_shift
        CIRC(i) = sin(Normalized_theta_I(i) - CIRC_shift) + s * sin(2 * (Normalized_theta_I(i) - CIRC_shift));
        if CIRC(i) < 0
            CIRC(i) = 0;
        end
        if a < 1
            CIRC(i) = CIRC(i) * a;
        end
    else
        CIRC(i) = -sin(2 * pi - (Normalized_theta_I(i) - CIRC_shift)) -...
               s * sin(2 * pi - (2 * (Normalized_theta_I(i) - CIRC_shift)));
        if CIRC(i) > 0
            CIRC(i) = 0;
        end
        if a > 1
            CIRC(i) = CIRC(i) / a;
        end
    end

    if s < 1
        max_1 = -0.5669 * s^3 + 1.1431 * s^2 + 0.1703 * s + 0.9963;
        CIRC(i) = CIRC(i) / max_1;
    else
        max_2 = 0.0029 * s^2 + 0.9738 * s + 0.7783;
        CIRC(i) = CIRC(i) / max_2;
    end
    
    % dCIRC Model ODEs
    dtheta_I = omega_I(i); 
    domega_I(i) = light(i) * CIRC(i) + elastic_factor * (omega - omega_I(i)); 
    derivative = [dtheta_I, domega_I(i)]; 
end
end
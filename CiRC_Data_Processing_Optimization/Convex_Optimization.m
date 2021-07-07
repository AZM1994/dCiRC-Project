%%% Convex optimization (PSO) algorithm to find the best fitting dCIRC parameters
%%% Written by Zheming An
clear;close all;
tic

global Normalized_CIRC_data omega_DD T_interval Phase_Data Time_Data number_of_timesteps ...
    cost_range T_interval_stable alpha_1 alpha_2 alpha_3 alpha_4
%% Load normalized race tube data
Normalized_CIRC_data = load('JW18_II_Normalized.mat');
Normalized_CIRC_data = Normalized_CIRC_data.data_saved;
cost_range           = 10:625; % define the range for cost function
Time_Data            = Normalized_CIRC_data(:,1);
Phase_Data           = Normalized_CIRC_data(:,2);
number_of_timesteps  = length(Phase_Data);

%% PSO Convex Optimization
for ii = 1 : number_of_timesteps - 1
	T_interval(ii) = Normalized_CIRC_data(ii+1,1) - Normalized_CIRC_data(ii,1);
end
T_interval_stable = max(T_interval);

% Parameters
alpha_1             = 0.1; % paramter for convex cost function
alpha_2             = 0.1; % paramter for convex cost function
alpha_3             = 0.1; % paramter for convex cost function
alpha_4             = 0.1; % paramter for convex cost function
tau_DD              = 21.7; % free running period for the selected genotype
omega_DD            = 2 * pi ./ tau_DD; % free running velocity for the selected genotype

% set optimization ranges for the dCIRC parameters
% the order is: ["cz", "elastic fator", "s", "a", "initial theta_I"]
% "cz" is the zeitgeber strength
% "s" is the shape factor of the CIRC
% "a" is the asymmetry factor of the CIRC
% "initial theta_I" is the initial phase angle of the internal clock
lower_bound         = [0   0   0   0   0]; % lower bound
upper_bound         = [3   3   2   2   8]; % upper bound
number_of_variables = 5; % number of optimation variables

% set initial condition for the optimation variables
initial_swarm       = [0 0 0 1 0]; % Default initial conditions
options = optimoptions(@particleswarm,'PlotFcn','pswplotbestf','SwarmSize',100,'HybridFcn',@fmincon,...
          'FunctionTolerance', 10^-6, 'MaxStallIterations', 20, 'InitialSwarmMatrix',initial_swarm);
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
    alpha_1 alpha_2 alpha_3 alpha_4
% preallocate 
ALL_clocks          = zeros(number_of_timesteps,3); % matrix contains all clock data
Normalized_theta_Z  = zeros(number_of_timesteps,1);
Normalized_theta_I  = zeros(number_of_timesteps,1);
CIRC                = zeros(number_of_timesteps,1);
light               = zeros(number_of_timesteps,1);
omega_I             = zeros(number_of_timesteps,1);
domega_I            = zeros(number_of_timesteps,1);

% Parameters
zeit_t          = 24; % period of zeitgeber clock
c_z             = x(1); % zeitgeber strength
elastic_factor  = x(2); % elastic factor
s               = x(3); % shape factor of the CIRC
a               = x(4); % asymmetry factor of the CIRC
initial_theta_I = x(5); % initial phase angle of the internal clock
initial_omega_I = omega_DD; % free running velocity for the selected genotype
CIRC_shift      = 0; % the time shift of the CIRC curve

% The desinated light profiles, given by center of the light explosure, 
% and radius of the light explosure
light_profiles       = [1,0;1,12;6,12;1,7;1.5,7.5;3,9;4.5,10.5]* 2 * pi / zeit_t; 
% convert zeitgeber time to the phase angle, pi / zeit_t corresponds to 30 mins 

light_profiles_index = 3; % choose the light profile
epsilon              = light_profiles(light_profiles_index,1); % radius of the light explosure
light_middle         = light_profiles(light_profiles_index,2); % center of the light explosure

%% Simulate the entrainment process
theta_I_initial = initial_theta_I; 
omega_I_initial = initial_omega_I; 
initial_zeit    = 0; % initial phase angle of zeitgeber clock
zeitgeber_initial = initial_zeit;
initial_conditions = [zeitgeber_initial, theta_I_initial, omega_I_initial];
ALL_clocks(1,:) = initial_conditions;

% Runge-kutta forth order iterative method
for i = 1:number_of_timesteps - 1
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
theta_I_signal = sin(ALL_clocks(:,2));
fitting_sin_curve = theta_I_signal;

% define the weight function for cz and elastic factor
if c_z < 0.3369
    cz_weight_func = (0.3369 - c_z)/c_z;
else
    cz_weight_func = (c_z-0.3369)^2;
end

if elastic_factor < 1.7346
    ke_weight_func = 0.02*(1.7346 - elastic_factor)/elastic_factor;
else
    ke_weight_func = (elastic_factor-1.7346)^2;
end

% define the convex cost function
convex_cost = sqrt(sum((fitting_sin_curve(cost_range) - Phase_Data(cost_range)).^2)) ...
               + alpha_1 * cz_weight_func ...
               + alpha_2 * ke_weight_func ...
               + alpha_3 * s/(2.01-s)...
               + alpha_4 * ((a-1.81)./(2.01-a)+(0.19-a)/(a+0.01) + 1.62);
               
function [ derivative ] = RK4_circadian_CIRC(y)
    zeit_clock = y(1,1);
    theta_I = y(1,2);
    omega_I(i) = y(1,3);
    Normalized_theta_Z(i) = mod(zeit_clock,2*pi);
    Normalized_theta_I(i) = mod(theta_I,2*pi);
    % Light function
    if abs(Normalized_theta_Z(i) - light_middle) < epsilon 
            light(i) = c_z;
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
    dzeit = 2*pi/zeit_t;
    dtheta_I = omega_I(i); 
    domega_I(i) = light(i) * CIRC(i) + elastic_factor * (omega_DD - omega_I(i)); 
    derivative = [dzeit, dtheta_I, domega_I(i)]; 
        
end
end

%% Plot the best fitting curves
function [L2_cost, convex_cost] = CIRC_LD_Cost_for_plot_best_fitting(optimized_x, minimum_cost)

global omega_DD T_interval Phase_Data Time_Data number_of_timesteps cost_range ...
    alpha_1 alpha_2 alpha_3 alpha_4
% preallocate
ALL_clocks                  = zeros(number_of_timesteps,3);
Normalized_theta_Z          = zeros(number_of_timesteps,1);
Normalized_theta_I          = zeros(number_of_timesteps,1);
CIRC                        = zeros(number_of_timesteps,1);
light                       = zeros(number_of_timesteps,1);
omega_I                     = zeros(number_of_timesteps,1);
domega_I                    = zeros(number_of_timesteps,1);

% Parameters
zeit_t          = 24;
c_z             = optimized_x(1);
elastic_factor  = optimized_x(2);
s               = optimized_x(3);
a               = optimized_x(4);
initial_theta_I = optimized_x(5);
initial_omega_I = omega_DD;
CIRC_shift      = 0;

% light profiles
light_profiles       = [1,0;1,12;6,12;1,7;1.5,7.5;3,9;4.5,10.5]* 2 * pi / zeit_t; % pi/zeit_t corresponds to 30 mins
light_profiles_index = 3;
epsilon              = light_profiles(light_profiles_index,1);
light_middle         = light_profiles(light_profiles_index,2);

%% Simulations
theta_I_initial = initial_theta_I;
omega_I_initial = initial_omega_I;
initial_zeit = 0;
zeitgeber_initial = initial_zeit;
initial_conditions = [zeitgeber_initial, theta_I_initial, omega_I_initial];
ALL_clocks(1,:) = initial_conditions;

for i = 1:number_of_timesteps - 1
    k1= RK4_circadian_CIRC(ALL_clocks(i,:));
    k2_clocks = ALL_clocks(i,:) + 0.5*T_interval(i)*k1;
    
    k2= RK4_circadian_CIRC(k2_clocks);
    k3_clocks = ALL_clocks(i,:) + 0.5*T_interval(i)*k2;
    
    k3= RK4_circadian_CIRC(k3_clocks);
    k4_clocks = ALL_clocks(i,:) + T_interval(i)*k3;
    
    k4= RK4_circadian_CIRC(k4_clocks);
    ALL_clocks(i+1,:) =  ALL_clocks(i,:) + T_interval(i)*(k1 + 2*k2 + 2*k3 + k4) / 6.0;
end

theta_I_signal = sin(ALL_clocks(:,2));

fitting_sincurve_P = theta_I_signal; 
if c_z < 0.3369 
    cz_weight_func = (0.3369 - c_z)/c_z;
else
    cz_weight_func = (c_z-0.3369)^2;
end

if elastic_factor < 1.7346
    ke_weight_func = 0.02*(1.7346 - elastic_factor)/elastic_factor;
else
    ke_weight_func = (elastic_factor-1.7346)^2;
end

convex_cost = sqrt(sum((fitting_sincurve_P(cost_range) - Phase_Data(cost_range)).^2)) ...
               + alpha_1 * cz_weight_func ...
               + alpha_2 * ke_weight_func ...
               + alpha_3 * s/(2.01-s)...
               + alpha_4 * ((a-1.81)./(2.01-a)+(0.19-a)/(a+0.01) + 1.62);
L2_cost = sqrt(sum((fitting_sincurve_P(cost_range) - Phase_Data(cost_range)).^2));

% find the intial angle the internal clock
[Af, Bf] = min(abs(Phase_Data(cost_range(1:80))));
actural_initial = fitting_sincurve_P(cost_range(Bf));
assignin('base','Af',Af)
assignin('base','Bf',Bf)
assignin('base','actural_initial',actural_initial)
    
    %% plot the best fitting curves
    figure(1);
        plot(Time_Data, Phase_Data, 'k', 'LineWidth',6)
        hold on
        plot(Time_Data(cost_range), fitting_sincurve_P(cost_range), 'r', 'LineWidth',6)
        area(Time_Data,zeros(number_of_timesteps,1) + 1, 'FaceColor', 'k', 'FaceAlpha', 1, 'LineStyle','none')
        area(Time_Data,zeros(number_of_timesteps,1) - 1, 'FaceColor', 'k', 'FaceAlpha', 1, 'LineStyle','none')
        area(Time_Data, light/c_z,...
            'FaceColor', 'y', 'FaceAlpha', 1, 'LineStyle','none')
        area(Time_Data, -light/c_z,...
            'FaceColor', 'y', 'FaceAlpha', 1, 'LineStyle','none')
        alpha(0.3)
        set(gca,'FontSize',42)
        set(gcf,'Position',[550 80 700 450]);
        disp([L2_cost, convex_cost, minimum_cost])
        xlabel('time/hours','Interpreter','latex')
        ylabel('rhythmic signal','Interpreter','latex')
        xlim([Time_Data(1) Time_Data(end)])
        ylim([-1 1])
        xticks(24:12:120)
        xticklabels({'24','36','48','60','72','84','96','108','120'})
        title('')
        legend('Normalized Data', 'Best Fitting Curve', 'AutoUpdate','off')

    % Plot the clock running velocity during the entrainment
    figure(2);
        plot(Time_Data(1:end-1), omega_I(1:end-1),'r', 'LineWidth',6);
        legend('angular velocity of internal clock ($\omega_I$)', 'Location', 'best', 'AutoUpdate','off','Interpreter','latex')
        hold on
        plot(Time_Data(1:end-1), zeros(1, number_of_timesteps - 1) + 2*pi/zeit_t, '--k', ...
             Time_Data(1:end-1), zeros(1, number_of_timesteps - 1) + omega_DD, '--b', 'LineWidth',6);
        area(Time_Data,zeros(number_of_timesteps,1) + 1, 'FaceColor', 'k', 'FaceAlpha', 1, 'LineStyle','none')
        area(Time_Data, light/c_z,'FaceColor', 'y', 'FaceAlpha', 1, 'LineStyle','none')
        alpha(0.3)
        xlabel('time / hours','Interpreter','latex')
        ylabel('$\omega_I/rad \cdot hour^{-1}$','Interpreter','latex')
        xticks(24:12:120)
        xticklabels({'24','36','48','60','72','84','96','108','120'})
        xlim([Time_Data(1) Time_Data(end)])
        ylim([min(nonzeros(omega_I)) - 0.03 max(omega_I) + 0.03])
        text(Time_Data(end)-5, 2*pi/zeit_t-0.01, '$\omega_T=\frac{2\pi}{T}$','Interpreter','latex','FontSize',42)
        text(Time_Data(end)-5, omega_DD+0.01, '$\omega_{I0}=\frac{2\pi}{\tau_{I0}}$','Interpreter','latex','FontSize',42)
        set(gca,'fontsize',42);
        set(gcf,'Position',[0 480 700 350]);
        
function [ derivative ] = RK4_circadian_CIRC(y)
    zeit_clock = y(1,1);
    theta_I = y(1,2);
    omega_I(i) = y(1,3);
    Normalized_theta_Z(i) = mod(zeit_clock,2*pi);
    Normalized_theta_I(i) = mod(theta_I,2*pi);
    
    % Light function
    if abs(Normalized_theta_Z(i) - light_middle) < epsilon 
            light(i) = c_z;
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
    
    % dCRIC Model ODEs
        dzeit = 2*pi/zeit_t;
        dtheta_I = omega_I(i); 
        domega_I(i) = light(i) * CIRC(i) + elastic_factor * (omega_DD - omega_I(i)); 
        derivative = [dzeit, dtheta_I, domega_I(i)];
end
end

%% Plot the best fitting dCIRC curve
function CIRC_LD_Cost_for_plot_CIRC_Curve(optimized_x)
T_interval      = 0.01; 
T_initial       = 0; 
T_final         = 24; 
T_final_plot    = T_final - T_interval; 
%%% preallocate 
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
c_z            = 0.001; 
elastic_factor  = 0.001;
s               = optimized_x(3); 
a               = optimized_x(4); 
initial_theta_I = 0;
initial_omega_I = 2 * pi ./ zeit_t; 
CIRC_shift      = 0; 

% the light profiles
light_profiles       = [1,0;1,12;6,12;1,7;1.5,7.5;3,9;4.5,10.5] * 2 * pi / zeit_t; % dawn dusk allday 
light_profiles_index = 3; 
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
        fill([T_initial:T_interval:T_final_plot flip(T_initial:T_interval:T_final_plot)],...
            [transpose(light/(c_z)) zeros(size(transpose(light/(c_z))))],'y','LineStyle','-.')
        hold on
        plot(T_initial:T_interval:T_final_plot, CIRC, 'k','LineWidth',5) 
        fill([T_initial:T_interval:T_final_plot flip(T_initial:T_interval:T_final_plot)],...
            [transpose(-light/(c_z)) zeros(size(transpose(light/(c_z))))],'y','LineStyle','-.')
        xlabel('time/hours')
        alpha(0.25)
        legend('light','CIRC')

        xticks(24:12:120)
        xticklabels({'24','36','48','60','72','84','96','108','120'})
        ylim([-1.0 1.0])
        title({['\color{black}CIRC(a=', num2str(round(a,2)), ', s=', num2str(round(s,2)),')']})
        set(gca,'fontsize',24);
        set(gcf,'Position',[0 80 720 350]);
        grid on
        grid minor
function [ derivative ] = RK4_circadian_CIRC(y)
    theta_I = y(1,1);
    omega_I(i) = y(1,2);
    Normalized_theta_Z(i) = mod(zeitgeber_clock(i),2*pi);
    Normalized_theta_I(i) = mod(theta_I,2*pi);
    
    % Light function
    if abs(Normalized_theta_Z(i) - light_middle) < epsilon 
            light(i) = c_z;
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
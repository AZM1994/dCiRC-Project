%%% Explore the convex cost - s - a space to find the best fitting (s,a)
%%% combination which minimize the convex cost
%%% Written by Zheming An

clear;close all;
tic
global omega_DD T_interval Phase_Data days_P number_of_timesteps C_ZG ...
       Light_Profile elastic_factor cost_range T_interval_stable s a alpha_1 alpha_2 alpha_3 alpha_4
% calculate cost for all (s,a) combinations
list_s = 0:0.05:2;
list_a = 0:0.05:2;
% calculate cost for a pair of (s,a)
% list_s = 1.15;
% list_a = 0.15;
total_count = length(list_s) * length(list_a);
list_s_cost  = round(list_s * 20 + 1);
list_a_cost  = round(list_a * 20 + 1);
count  = 0;

%% Load data
genotype_name   = 'HG2';
HG_data   = load('linear-HG2-0.1.mat');
HG_data   = HG_data.data_saved;
Light_Profile = load('linear-Light_Profile.mat');
Light_Profile = Light_Profile.Light_Profile;

previous_cost_matrix = load(['Weighted_cost_matrix_',genotype_name,'.mat']);
previous_cost_matrix = previous_cost_matrix.cost_matrix;
cost_matrix  = previous_cost_matrix;
cost_range      = 1:139;
days_P          = HG_data(:,1);
Phase_Data      = HG_data(:,2);
number_of_timesteps = length(Phase_Data);
%% Parameters
for ii = 1 : number_of_timesteps - 1
    T_interval(ii) = HG_data(ii+1,1) - HG_data(ii,1);
end                            
T_interval_stable = max(T_interval);

alpha_1             = 0;
alpha_2             = 0;
alpha_3             = 0;
alpha_4             = 0;
C_ZG                = 0.0287; % best fitting zeitgeber strength
elastic_factor      = 2.7933; % best fitting elastic factor
tau_DD              = 24.72; % free running period of the internal clock
omega_DD            = 2 * pi ./ tau_DD; % free running velocity of the internal clock
lower_bound         = 0; % lower bound of initial phase angle of the internal clock
upper_bound         = 2*pi; % upper bound of initial phase angle of the internal clock
number_of_variables = 1; % number of optimizing variables
initial_swarm       = 5; % initial value of initial phase angle of the internal clock



%% Simulations
for s_index = 1:size(list_s,2)
    s = list_s(s_index);
    for a_index = 1:size(list_a,2)
    a = list_a(a_index);
        count = count + 1;
        [count,total_count]
        [s,a]
        options = optimoptions(@particleswarm,'PlotFcn','pswplotbestf','SwarmSize',10,'HybridFcn',@fmincon,...
                  'FunctionTolerance', 10^-6, 'MaxStallIterations', 3, 'InitialSwarmMatrix',initial_swarm); 
        [optimized_x ,Fval,exitFlag,Output] = particleswarm(@CIRC_LD_Cost,number_of_variables,lower_bound,upper_bound,options); 
        minimum_cost = CIRC_LD_Cost(optimized_x)
        if minimum_cost < previous_cost_matrix(list_s_cost(s_index),list_a_cost(a_index))
            cost_matrix(list_s_cost(s_index),list_a_cost(a_index)) = minimum_cost;
        end
    end 
end

% save the convex cost matrix
save(['Weighted_cost_matrix_',genotype_name,'.mat'],'cost_matrix')

toc


%% Calculate the convex cost
function convex_cost = CIRC_LD_Cost(x)

global omega_DD T_interval Light_Profile Phase_Data number_of_timesteps C_ZG elastic_factor cost_range T_interval_stable s a alpha_1 alpha_2 alpha_3 alpha_4
%%% preallocate 
ALL_clocks          = zeros(number_of_timesteps,2); 
zeitgeber_clock     = zeros(number_of_timesteps,1); 
Normalized_theta_Z  = zeros(number_of_timesteps,1); 
Normalized_theta_I  = zeros(number_of_timesteps,1); 
CIRC                = zeros(number_of_timesteps,1); 
light               = zeros(number_of_timesteps,1); 
omega_I             = zeros(number_of_timesteps,1); 
domega_I            = zeros(number_of_timesteps,1); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zeit_t          = 24; 
amplitude       = 1; 
initial_theta_I = x(1); 
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

    % COST
    fitting_sincurve_P = amplitude .* theta_I_signal; 
    convex_cost = sqrt(sum((fitting_sincurve_P(cost_range) - Phase_Data(cost_range)).^2)) ...
                   + alpha_1 * (C_ZG - 0.5579)^2 ...
                   + alpha_2 * (elastic_factor - 0.9595)^2 ...
                   + alpha_3 * s^2 ...
                   + alpha_4 * (a-1)^2;
               
function [ derivative ] = RK4_circadian_CIRC(y)
    theta_I = y(1,1);
    omega_I(i) = y(1,2);
    Normalized_theta_Z(i) = mod(zeitgeber_clock(i),2*pi);
    Normalized_theta_I(i) = mod(theta_I,2*pi);
    
    % Light profile
    light(i) = C_ZG * Light_Profile(i);
    
    % CiRC function
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
    
    % dCiRC Model ODEs
        dtheta_I = omega_I(i); 
        domega_I(i) = light(i) * CIRC(i) + elastic_factor * (omega_DD - omega_I(i)); 
        derivative = [dtheta_I, domega_I(i)]; 
end
end
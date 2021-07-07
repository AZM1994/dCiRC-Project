%%% Generating Fig. 1 A, and B for the illustrative example of the dCiRC model,
%%% reflecting both the parametric and non-parametric aspects of the
%%% entrainment under cycling conditions.
%%% Fig. 1 A: The trajectory of the internal clock's period change.
%%% Fig. 1 B: The trajectory of the internal clock's velocity change.
%%% Written by Zheming An

function Fig1AB_light_exp_period_velocity
clear, close all

% Parameters
T_interval      = 0.01; % unit in hours
T_initial       = 0; % unit in hours
T_final         = 1680; % unit in hours, time duration of the entrainment experiment
light_period    = 840; % unit in hours, length of each rounds of the zeitgeber cycles
light_off       = 480; % unit in hours, the time point when transit into DD conditions in each round
light_on        = 0; % unit in hours, the time point when transit into LD conditions in each round
T_final_plot    = T_final - T_interval; 

cz              = 0.0025; % zeitgeber strength
elastic_factor  = 0.05; % elastic factor
zeit_t          = 24; % period of zeitgeber clock
tau_DD          = 23; % free running period of the internal clock
omega_DD        = 2 * pi ./ tau_DD; % free runing velocity of the internal clock
CIRC_shift      = 0; % shifting time of the CIRC curve
phi1            = 0; % initial phase angle of the internal clock

% CIRC parameters
S_A_combination = [0,1; 0.5,2; 2,1; 0.5,0.5; 0.6,1.0];
CIRC_type       = 1;
s   = S_A_combination(CIRC_type,1); % shape factor
a   = S_A_combination(CIRC_type,2); % asymmetry factor

% zeitgeber cycles: light profile
% unit in hours, first column is the center time of the light,
% second column is the radius of the light explosure
light_profiles  = [1,0;1,12;6,12;7,1] * 2 * pi / zeit_t; % convert time to phase angles
light_profiles_index = 3; 
epsilon    = light_profiles(light_profiles_index,1); % the radius of the light explosure
light_middle    = light_profiles(light_profiles_index,2); % the center time of the light

% preallocate 
number_of_timesteps         = (T_final - T_initial) / T_interval; 
ALL_clocks                  = zeros(number_of_timesteps,2); % matrix stores all simulated clock parameters
zeitgeber_clock             = zeros(number_of_timesteps,1);
Normalized_theta_Z          = zeros(number_of_timesteps,1);
Normalized_theta_I          = zeros(number_of_timesteps,1);
CIRC                        = zeros(number_of_timesteps,1);
light                       = zeros(number_of_timesteps,1);
omega_I                     = zeros(number_of_timesteps,1);
domega_I                    = zeros(number_of_timesteps,1);
analytical_phase_difference = zeros(floor(T_final/zeit_t),1);



%% Simulations    
theta_I_initial = phi1; % initial phase angle of the interanl clock
omega_I_initial = omega_DD; % initial velocity of the interanl clock
initial_conditions = [theta_I_initial, omega_I_initial];
ALL_clocks(1,:) = initial_conditions;

% Runge-Kutta forth order method
for i = 1:number_of_timesteps - 1
    zeitgeber_clock(i) = (i-1)*2*pi/(zeit_t/T_interval);
    k1= RK4_circadian_CIRC(ALL_clocks(i,:));
    k2_clocks = ALL_clocks(i,:) + 0.5*T_interval*k1';

    k2= RK4_circadian_CIRC(k2_clocks);
    k3_clocks = ALL_clocks(i,:) + 0.5*T_interval*k2';

    k3= RK4_circadian_CIRC(k3_clocks);
    k4_clocks = ALL_clocks(i,:) + T_interval*k3';

    k4= RK4_circadian_CIRC(k4_clocks);
    ALL_clocks(i+1,:) =  ALL_clocks(i,:) + T_interval*(k1' + 2*k2' + 2*k3' + k4') / 6.0;             
end
% sinusoidal of phase angles of the internal and zeitgeber clocks
    analytical_zeitgeber_signal = sin(zeitgeber_clock);
    theta_I_signal = sin(ALL_clocks(:,1));

%% find the time values corresponding to the peaks of the signal
[~, n_time]= findpeaks(theta_I_signal,'MinPeakHeight',0.99);
[~, analytical_z_time]= findpeaks(analytical_zeitgeber_signal,'MinPeakHeight',0.99);

analytical_days = size(analytical_z_time,1);

% Calculate the phase of entrainment
for j = 1:analytical_days
    [~, analytical_phase_index] = min(abs(analytical_z_time(j) - n_time));
    analytical_phase_difference(j) = (analytical_z_time(j) - n_time(analytical_phase_index)) .* T_interval;
end

% save data for the double plot
save('analytical_phase_difference','analytical_phase_difference')
save('analytical_z_time','analytical_z_time')

%% Plot
% Plot tau:
% the period changes of the internal clock during the entrainment
figure(1)
    tau_I = 2 * pi ./ omega_I;
    area(T_initial:T_interval:(T_final_plot), (zeros(length(T_initial:T_interval:(T_final_plot)),1) + 1) * 25, ...
        'FaceColor', 'k', 'FaceAlpha', 0.3, 'LineStyle','none')
    hold on
    area(T_initial:T_interval:(T_final_plot), light/cz * 25, 'FaceColor', 'y', 'FaceAlpha', 1, 'LineStyle','none')

    plot(T_initial:T_interval:(T_final_plot-T_interval), tau_I(1:end-1), 'b', 'LineWidth',1.5);
    plot(T_initial:T_interval:(T_final_plot-T_interval),repelem(zeit_t,number_of_timesteps-1)','--k',...
         T_initial:T_interval:(T_final_plot-T_interval),repelem(2*pi/omega_DD,number_of_timesteps-1)','--k','LineWidth',2,'HandleVisibility','off');
    plot(1680-[192 192],[0 48],'r','LineWidth',3)
    plot(1680-[360 360],[0 48],'r','LineWidth',3)
    plot(1680-[504 504],[0 48],'r','LineWidth',3)
    plot(1680-[840 840],[0 48],'r','LineWidth',3)
    plot(1680-[1032 1032],[0 48],'r','LineWidth',3)
    plot(1680-[1200 1200],[0 48],'r','LineWidth',3)
    plot(1680-[1416 1416],[0 48],'r','LineWidth',3)
    text(1685,tau_DD, '$\tau_{I0}$', 'Interpreter', 'latex','Fontsize',36, 'Color', 'k', 'FontWeight', 'bold')
    text(1685,24, '$T$', 'Interpreter', 'latex', 'Fontsize',36, 'Color', 'k', 'FontWeight', 'bold')

    xlabel('days of simulation','Interpreter','latex')
    ylabel('period $\tau_I$/hours','Interpreter','latex')
    xticks([0 264 480 648 840 1176 1320 1488 1680])
    xticklabels({'0','11','20','27','35','49','55','62','70'})
    xlim([0 T_final])
    ylim([21 25])
    set(gca,'fontsize',36);
    set(gca,'LineWidth',3)

% plot omega:
% the velocity changes of the internal clock during the entrainment
figure(2)
    area(T_initial:T_interval:(T_final_plot), (zeros(length(T_initial:T_interval:(T_final_plot)),1) + 1) * 0.35, ...
        'FaceColor', 'k', 'FaceAlpha', 0.3, 'LineStyle','none')
    hold on
    area(T_initial:T_interval:(T_final_plot), light/cz * 0.35, 'FaceColor', 'y', 'FaceAlpha', 1, 'LineStyle','none')
    plot(T_initial:T_interval:(T_final_plot-T_interval), omega_I(1:end-1), 'b','LineWidth',2);
    plot(T_initial:T_interval:(T_final_plot-T_interval), repelem(omega_DD, number_of_timesteps-1), '--k', 'LineWidth',2);
    plot(T_initial:T_interval:(T_final_plot-T_interval), repelem(2*pi/24, number_of_timesteps-1), '--k', 'LineWidth',2);
    plot(1680-[192 192],[0 48],'r','LineWidth',3)
    plot(1680-[360 360],[0 48],'r','LineWidth',3)
    plot(1680-[504 504],[0 48],'r','LineWidth',3)
    plot(1680-[840 840],[0 48],'r','LineWidth',3)
    plot(1680-[1032 1032],[0 48],'r','LineWidth',3)
    plot(1680-[1200 1200],[0 48],'r','LineWidth',3)
    plot(1680-[1416 1416],[0 48],'r','LineWidth',3)
    text(1685,omega_DD, '$\omega_{I0}$', 'Interpreter', 'latex','Fontsize',36, 'Color', 'k', 'FontWeight', 'bold')
    text(1685,2*pi/24, '$\frac{2\pi}{T}$', 'Interpreter', 'latex', 'Fontsize',36, 'Color', 'k', 'FontWeight', 'bold')

    xlabel('days of simulation','Interpreter','latex')
    ylabel('angular velocity $\omega_I/hour^{-1}$','Interpreter','latex')
    xticks([0 264 480 648 840 1176 1320 1488 1680])
    xticklabels({'0','11','20','27','35','49','55','62','70'})
    xlim([0 T_final])
    ylim([0.25 0.3])
    set(gca,'fontsize',36);
    set(gca,'LineWidth',3)

% Plot phase of entrainment
figure(3);
    area(T_initial:T_interval:(T_final_plot), (zeros(length(T_initial:T_interval:(T_final_plot)),1) + 1) * 12, ...
        'FaceColor', 'k', 'FaceAlpha', 0.3, 'LineStyle','none')
    hold on
    area(T_initial:T_interval:(T_final_plot), (zeros(length(T_initial:T_interval:(T_final_plot)),1) + 1) * -12, ...
        'FaceColor', 'k', 'FaceAlpha', 0.3, 'LineStyle','none')
    area(T_initial:T_interval:(T_final_plot), light/cz * 12, 'FaceColor', 'y', 'FaceAlpha', 1, 'LineStyle','none')
    area(T_initial:T_interval:(T_final_plot), light/cz * -12, 'FaceColor', 'y', 'FaceAlpha', 1, 'LineStyle','none')
    plot(analytical_z_time*T_interval,analytical_phase_difference, 'b', 'LineWidth',4);
    plot(T_initial:T_interval:(T_final_plot-T_interval),repelem(2.96,number_of_timesteps-1)','--k','LineWidth',2,'HandleVisibility','off')
    plot(1680-[192 192],[-12 12],'r','LineWidth',3)
    plot(1680-[360 360],[-12 12],'r','LineWidth',3)
    plot(1680-[504 504],[-12 12],'r','LineWidth',3)
    plot(1680-[840 840],[-12 12],'r','LineWidth',3)
    plot(1680-[1032 1032],[-12 12],'r','LineWidth',3)
    plot(1680-[1200 1200],[-12 12],'r','LineWidth',3)
    plot(1680-[1416 1416],[-12 12],'r','LineWidth',3)
    text(1685,2.96, '$\Psi$', 'Interpreter', 'latex','Fontsize',36, 'Color', 'k', 'FontWeight', 'bold')
    
    xticks([0 264 480 648 840 1176 1320 1488 1680])
    xticklabels({'0','11','20','27','35','49','55','62','70'})
    yticks([-12 -9 -6 -3 0 3 6 9 12])
    yticklabels({'-12','-9','-6','-3','0','3','6','9','12'})
    xlabel('days of simulation','Interpreter','latex')
    ylabel('External Time /hours')
    xlim([0 T_final])
    ylim([-12 12])
    set(gca,'xaxisLocation', 'bottom')
    set(gca,'fontsize',36);
    set(gca,'LineWidth',3)

% Runge-Kutta forth order method and model ODEs
function [ derivative ] = RK4_circadian_CIRC(y)
    theta_I = y(1,1);
    omega_I(i) = y(1,2);
    Normalized_theta_Z(i) = mod(zeitgeber_clock(i),2*pi);
    Normalized_theta_I(i) = mod(theta_I,2*pi);

% Light profile
if abs(Normalized_theta_Z(i) - light_middle) < epsilon 
    if mod(i,light_period/T_interval) > light_on/T_interval && mod(i,light_period/T_interval) < light_off / T_interval
        light(i) = cz;
    else
        light(i) = 0;
    end
else
    light(i) = 0;
end

% dCiRC functions
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
    derivative = [dtheta_I, domega_I(i)]';
end
end
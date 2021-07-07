%%% Predict phase of entrainment with the best fitting dCIRC parameters
%%% Written by Zheming An

clear;close all;
tic

%% Load the best fitting dCIRC parameter set
load_parameter_set = load('best_fitting_dCiRC_parameter_set.mat');
best_fitting_dCiRC_parameter_set = load_parameter_set.optimized_x_POE;
organism_name_list = ["JW168", "JW224", "JW60",  "JW260", "D117",...
                      "D119",  "D116",  "JW18",  "JW172", "P4283",...
                      "JW220", "JW200", "JW162", "JW176", "JW238",...
                      "JW261", "JW180", "JW169", "JW228", "P4463",...
                      "JW24",  "JW59",  "JW161", "P4469", "JW22"];

POE_range_list = {1:2,1:3,1:2,1:2,1:2,1:3,1:2,1:4,2:4,1:3,...
                  2:4,1:3,1:4,2:4,1:4,1:4,2:4,1:2,1:3,1:4,...
                  1:2,1:3,1:2,1:2,1:4};

results = ["genotype","entrained_time","POE"];

for genotype_index = 1:25
    optimized_parameters     = best_fitting_dCiRC_parameter_set(genotype_index,:);
    organism_names  = organism_name_list(genotype_index);
    POE_range       = cell2mat(POE_range_list(genotype_index));
    results_temp    = [organism_names,CIRC_LD_Cost_for_plot_best_fitting(optimized_parameters,POE_range)];
    results         = [results;results_temp];
end
    
    results
    
toc

function POE_ROE_Avg = CIRC_LD_Cost_for_plot_best_fitting(optimized_parameters,POE_range)

    % Parameters
    T_initial       = 0;
    T_interval      = 0.001;
    T_final         = optimized_parameters(8);
    number_of_timesteps = (T_final - T_initial) / T_interval;
%     plot_time       = T_initial : T_interval : (T_final - T_interval);
    tau_DD          = optimized_parameters(1);
    cz              = optimized_parameters(2);
    elastic_factor  = optimized_parameters(3);
    s               = optimized_parameters(4);
    a               = optimized_parameters(5);
    phi_0           = optimized_parameters(6);
    zeit_shift      = optimized_parameters(7);
    zeit_t          = 24;
    [tau_DD, cz, elastic_factor, s, a]
    omega_DD        = 2 * pi ./ tau_DD;
    initial_theta_I = phi_0;
    initial_omega_I = omega_DD;
    CIRC_shift      = 0;

    % define light profiles
    light_profiles  = [6,12;12,12] * 2 * pi / zeit_t;
    light_profiles_index = 1;
    epsilon              = light_profiles(light_profiles_index,1); % radius of the entrainment window 
    light_middle         = light_profiles(light_profiles_index,2); % center of the entrainment window 
    
    %%% preallocate
    ALL_clocks                  = zeros(number_of_timesteps,2);
    zeitgeber_clock             = zeros(number_of_timesteps,1);
    Normalized_theta_Z          = zeros(number_of_timesteps,1);
    Normalized_theta_I          = zeros(number_of_timesteps,1);
    CIRC                        = zeros(number_of_timesteps,1);
    light                       = zeros(number_of_timesteps,1);
    omega_I                     = zeros(number_of_timesteps,1);
    domega_I                    = zeros(number_of_timesteps,1);
    
    %% Simulations
    theta_I_initial = initial_theta_I; 
    omega_I_initial = initial_omega_I; 
    initial_conditions = [theta_I_initial, omega_I_initial];
    ALL_clocks(1,:) = initial_conditions;

    for i = 1:number_of_timesteps - 1
        zeitgeber_clock(i) = (i-1-(18-zeit_shift)/T_interval)*2*pi/(zeit_t/T_interval);
        k1= RK4_circadian_CIRC(ALL_clocks(i,:));
        k2_clocks = ALL_clocks(i,:) + 0.5*T_interval*k1;

        k2= RK4_circadian_CIRC(k2_clocks);
        k3_clocks = ALL_clocks(i,:) + 0.5*T_interval*k2;

        k3= RK4_circadian_CIRC(k3_clocks);
        k4_clocks = ALL_clocks(i,:) + T_interval*k3;

        k4= RK4_circadian_CIRC(k4_clocks);
        ALL_clocks(i+1,:) =  ALL_clocks(i,:) + T_interval*(k1 + 2*k2 + 2*k3 + k4) / 6.0;             
    end
        analytical_zeitgeber_signal = sin(zeitgeber_clock);
        theta_I_signal = sin(ALL_clocks(:,1));

    % find the time values corresponding to the peaks of the signal
        [~, n_time]= findpeaks(theta_I_signal,'MinPeakHeight',0.99);
        [~, analytical_z_time]= findpeaks(analytical_zeitgeber_signal,'MinPeakHeight',0.99);

        analytical_days = size(analytical_z_time,1);
        analytical_phase_difference = zeros(analytical_days,1);

        for j = 1:analytical_days
            [~, analytical_phase_index] = min(abs(analytical_z_time(j) - n_time));
            analytical_phase_difference(j) = (n_time(analytical_phase_index) - analytical_z_time(j)) .* T_interval;
        end
%     % find 95% interval of the final POE
%         final_phase_of_entrainment = analytical_phase_difference(end - 1);
%         entrained_time = [];
%         phase_value_at_entrained_time = [];
%         if abs(analytical_phase_difference(analytical_days-1) - analytical_phase_difference(analytical_days-2)) < 1.5 * 0.01
%             for j = analytical_days-1:-1:1
%                 if abs(analytical_phase_difference(j) - final_phase_of_entrainment) > ...
%                     1.5 * 0.01
%                     phase_value_at_entrained_time = analytical_phase_difference(j+1);
%                     entrained_time = analytical_z_time(j+1) * T_interval;
%                     assignin('base','POE_ROE',[analytical_phase_difference(end),entrained_time])
% %                     entrained_day = j+1;
%                 break
%                 end
%             end
%         end
        entrained_time = round(T_final,2);
        average_POE = round(mean(analytical_phase_difference(POE_range)),2);
%         phase_value_at_entrained_time = analytical_phase_difference(end);
        POE_ROE_Avg = [entrained_time,average_POE];
        assignin('base','POE_ROE_Avg',POE_ROE_Avg)
        assignin('base','analytical_phase_difference',analytical_phase_difference)

    
        
% 	%% plot 
%     % plot phase angles of the internal clock and zeitgeber clock
%         figure(1);
%         plot(plot_time, analytical_zeitgeber_signal, '--k', 'LineWidth',3)
%         hold on
%         plot(plot_time, theta_I_signal', 'b', 'LineWidth',3)
%         area(plot_time, light/cz, 'FaceColor', 'y', 'FaceAlpha', 0.25, 'LineStyle','-.')
%         hold on
%         area(plot_time, -light/cz, 'FaceColor', 'y', 'FaceAlpha', 0.25, 'LineStyle','-.')
%         title_string = sprintf('tauDD = %0.1f, cz = %0.4f, elastic = %0.4f, s = %0.2f, a = %0.2f,\n initial thetaI = %0.4f',... 
%                 tau_DD, cz, elastic_factor, s, a, initial_theta_I);
%         set(gca,'FontSize',24)
%         xlabel('time/hours')
%         title(title_string)
%         legend('Zeitgeber', 'Internal Clock')
% % 
%     % plot the phase differences
%         figure(2);
%             plot(analytical_z_time*T_interval,analytical_phase_difference, 'LineWidth',3);
%             hold on
%             scatter(analytical_z_time*T_interval,analytical_phase_difference,100,'filled');
%             plot(plot_time,repelem(0,number_of_timesteps)','--k','LineWidth',2,'HandleVisibility','off')
%             xlabel('time/hour')
%             ylabel('phase of entrainment')
%             text(analytical_z_time(end)*T_interval-10, analytical_phase_difference(end)+0.3,...
%                  ['POE = ',num2str(analytical_phase_difference(end))],'FontSize',24)
%             text(entrained_time-15, phase_value_at_entrained_time+0.3,...
%                  ['entrained at t = ',num2str(round(entrained_time))],'FontSize',24)
%             set(gca,'xaxisLocation', 'bottom')
%             set(gca,'fontsize',24);
% 
%     
%     % plot the velocity changes of the internal clock during the entrainment
%         figure(3);
%             plot(plot_time(1:end-1), omega_I(1:end-1), 'r', 'LineWidth',2);
%             legend('omegaI')
%             hold on 
%             plot(plot_time,repelem(2*pi/zeit_t,number_of_timesteps)','--k',...
%                  plot_time,repelem(omega_DD,number_of_timesteps)','--k','LineWidth',2,'HandleVisibility','off')
%             area(plot_time, light/cz * max(omega_I(1:end-1)), 'FaceColor', 'y', 'FaceAlpha', 0.5, 'LineStyle','-.')
%             title('Angular speed of internal clock')
%             y_upper_omega = max(omega_I(1:end-1)) + 0.05;
%             y_lower_omega = min(omega_I(1:end-1)) - 0.05;
%             ylim([y_lower_omega y_upper_omega])
        
    function [ derivative ] = RK4_circadian_CIRC(y)
        theta_I = y(1,1);
        omega_I(i) = y(1,2);
        Normalized_theta_Z(i) = mod(zeitgeber_clock(i),2*pi);
        Normalized_theta_I(i) = mod(theta_I,2*pi);
        
        % Light profile
        if abs(Normalized_theta_Z(i) - light_middle) <= epsilon 
                light(i) = cz;
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
        
        % dCiRC Model ODEs
            dtheta_I = omega_I(i);
            domega_I(i) = light(i) * CIRC(i) + elastic_factor * (omega_DD - omega_I(i));
            derivative = [dtheta_I domega_I(i)];
    end
end
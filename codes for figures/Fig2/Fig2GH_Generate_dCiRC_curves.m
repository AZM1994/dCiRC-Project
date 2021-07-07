%%% Generate the best fitting dCiRC curves
%%% Written by Zheming An

function Fig2GH_Generate_dCiRC_curves
clear,close all
T_interval      = 0.01;
T_initial       = 0;
T_final         = 24;

% preallocate 
number_of_timesteps = length(T_initial : T_interval : T_final);
ALL_clocks          = zeros(number_of_timesteps,2); 
zeitgeber_clock     = zeros(number_of_timesteps,1); 
Normalized_theta_Z  = zeros(number_of_timesteps,1);
Normalized_theta_I  = zeros(number_of_timesteps,1); 
CIRC                = zeros(number_of_timesteps,1); 
light               = zeros(number_of_timesteps,1); 
omega_I             = zeros(number_of_timesteps,1); 
domega_I            = zeros(number_of_timesteps,1); 

% Best fitting parameters for generating dCiRC curves
Genotype_names = ["JW168", "JW224", "JW60",  "JW260", "D117",...
                  "D119",  "D116",  "JW18",  "JW172", "P4483",...
                  "JW220", "JW200", "JW162", "JW176", "JW238",...
                  "JW261", "JW180", "JW169", "JW228", "P4463",...
                  "JW24",  "JW59",  "JW161", "P4469", "JW22",...
                  "DBP338", "HG2"];
best_fitting_dCiRC_parameters = [
                21.0	0.6477	1.5737	1.40	0.05
                21.0	0.2436	0.7246	0.00	0.05
                21.0	0.3696	1.5247	0.05	0.10
                21.2	0.3548	1.9021	0.00	0.10
                21.3	1.4254	1.7593	1.05	0.05
                21.3	1.0021	1.9621	0.00	0.05
                21.5	0.2958	0.6504	0.05	0.05
                21.7	0.2658	0.4019	0.20	0.10
                21.7	0.3373	0.3216	1.00	0.05
                21.7	0.1559	0.1701	0.85	0.60
                21.8	0.1536	0.1060	1.40	0.40
                22.0	0.5679	0.5591	0.30	0.10
                22.2	0.9090	1.0826	0.40	0.25
                22.2	0.3997	1.1089	0.10	0.25
                22.2	0.3784	1.9092	0.00	0.57
                22.2	0.4937	1.9138	0.05	0.15
                22.3	0.3879	1.0788	0.20	0.40
                22.5	0.3429	0.7786	0.15	0.20
                22.7	0.6499	1.9883	0.00	0.10
                22.7	0.3258	0.7034	0.26	0.35
                23.2	0.2915	0.9880	0.20	0.10
                23.5	0.5339	0.6208	0.60	0.30
                23.8	0.7301	1.2864	0.00	0.05
                23.8	0.4537	0.8656	0.60	0.20
                22.0	0.3636	1.1602	1.00	0.14
                30.0	0.0394	0.0056	1.40	1.85
                24.72	0.0287	2.7933	0.00	2.00];

cluster_data_saved = zeros(size(best_fitting_dCiRC_parameters,1),9);

for index = 1 % choose any genotypes, e,g, 1 is JW168; 1:27 are all genotypes
    Genotype       = Genotype_names(index);
    tau_DD         = best_fitting_dCiRC_parameters(index,1); % free running period
    s              = best_fitting_dCiRC_parameters(index,4); % shape factor of CiRC
    a              = best_fitting_dCiRC_parameters(index,5); % asymmetry factor of CiRC

    cz             = 0.0001; % zeitgeber strength
    elastic_factor = 0; % elastic factor
    phi1           = 0; % initial phase angle for the internal clock
    phi2           = 2*pi/24; % initial phase angle for the internal clock
    zeit_t         = 24; % zeitgeber period
    tau_DD_sim     = 24; % free running period
    omega_DD_sim   = 2 * pi ./ tau_DD_sim; % free running velocity
    CIRC_shift     = 0;

    % light profiles: first column is the radius of the entrainment window
    % second column is the center of the entrainment window
    light_profiles  = [1,0;1,12;6,12;...
                       1,7;1.5,7.5;3,9;4.5,10.5];
    light_profiles_index = 3;
    epsilon    = light_profiles(light_profiles_index,1) * 2 * pi / zeit_t; % the radius of the entrainment window
    light_middle    = light_profiles(light_profiles_index,2) * 2 * pi / zeit_t; % the center of the entrainment window


    
    %% Simulations
    theta_I_initial = phi1;
    omega_I_initial = phi2;
    initial_conditions = [theta_I_initial, omega_I_initial];
    ALL_clocks(1,:) = initial_conditions;

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

    %% Plot best fitting dCiRC curves
    fh = figure(1);
%     fh.WindowState = 'maximized';
    plot(T_initial:T_interval:T_final, CIRC, 'k', 'LineWidth',8)
    hold on
    area(T_initial:T_interval:T_final,zeros(number_of_timesteps,1) + 1, 'FaceColor', 'k', 'FaceAlpha', 1, 'LineStyle','none')
    area(T_initial:T_interval:T_final,zeros(number_of_timesteps,1) - 1, 'FaceColor', 'k', 'FaceAlpha', 1, 'LineStyle','none')
    area(T_initial:T_interval:T_final, light/cz,'FaceColor', 'y', 'FaceAlpha', 1, 'LineStyle','none')
    area(T_initial:T_interval:T_final, -light/cz,'FaceColor', 'y', 'FaceAlpha', 1, 'LineStyle','none')
    alpha(0.3)

    %% calculate total area under curve (TAUC) and light-exposed area under curve (LAUC)
    Area_trapz_compress = trapz((light_profiles(light_profiles_index,2) - light_profiles(light_profiles_index,1)):T_interval:zeit_t/2,...
        CIRC((light_profiles(light_profiles_index,2) - light_profiles(light_profiles_index,1))/T_interval:zeit_t/(2*T_interval)));
    Area_trapz_compress_total = trapz(T_interval:T_interval:zeit_t/2, CIRC(1:zeit_t/(2*T_interval)));
    Area_trapz_expansion = trapz(zeit_t/2:T_interval:(light_profiles(light_profiles_index,2) + light_profiles(light_profiles_index,1)),...
        CIRC(zeit_t/(2*T_interval):(light_profiles(light_profiles_index,2) + light_profiles(light_profiles_index,1))/T_interval));
    Area_trapz_expansion_total = trapz(zeit_t/2:T_interval:zeit_t, CIRC(zeit_t/(2*T_interval):zeit_t/T_interval));
    Area_trapz_light = trapz((light_profiles(light_profiles_index,2) - light_profiles(light_profiles_index,1)):T_interval:...
        (light_profiles(light_profiles_index,2) + light_profiles(light_profiles_index,1)),...
        CIRC((light_profiles(light_profiles_index,2) - light_profiles(light_profiles_index,1))/T_interval:...
        (light_profiles(light_profiles_index,2) + light_profiles(light_profiles_index,1))/T_interval));
    Area_trapz_light_total = trapz(1 * T_interval : T_interval : zeit_t, CIRC(1 : zeit_t/T_interval));

    Area_compress = sum(CIRC((light_profiles(light_profiles_index,2) - light_profiles(light_profiles_index,1))/T_interval:...
        zeit_t/(2*T_interval)) .* T_interval);
    Area_compress_total = sum(CIRC((1:zeit_t/(2*T_interval))) .* T_interval);
    Area_expansion = sum(CIRC(zeit_t/(2*T_interval):...
        (light_profiles(light_profiles_index,2) + light_profiles(light_profiles_index,1))/T_interval) .* T_interval);
    Area_expansion_total = sum(CIRC((zeit_t/(2*T_interval):zeit_t/T_interval)) .* T_interval);
    Area_light = sum(CIRC((light_profiles(light_profiles_index,2) - light_profiles(light_profiles_index,1))/T_interval:...
        (light_profiles(light_profiles_index,2) + light_profiles(light_profiles_index,1))/T_interval) .* T_interval);
    Area_light_total = sum(CIRC((1:zeit_t/T_interval)) .* T_interval);

    Area_trapz = [Area_trapz_compress;Area_trapz_compress_total;Area_trapz_expansion;Area_trapz_expansion_total;Area_trapz_light;Area_trapz_light_total];
    Area_eular = [Area_compress;Area_compress_total;Area_expansion;Area_expansion_total;Area_light;Area_light_total];
    compress_ratio = Area_trapz_compress/Area_trapz_compress_total;
    expansion_ratio = Area_trapz_expansion/Area_trapz_expansion_total;
    Area_light_ratio = Area_light/Area_light_total;
    trapz_eular = [Area_trapz,Area_eular];
    
    cluster_data = [Area_trapz_compress, Area_trapz_compress_total, compress_ratio, ...
                    Area_trapz_expansion, Area_trapz_expansion_total, expansion_ratio,...
                    Area_light, Area_light_total, Area_light_ratio];
    cluster_data_saved(index,:) = cluster_data;
    cluster_data_02 = [compress_ratio, Area_trapz_compress, Area_trapz_compress_total, Area_light];
    assignin('base', 'trapz_eular', trapz_eular)
    assignin('base', 'cluster_data', cluster_data)
    assignin('base', 'cluster_data_02', cluster_data_02)
    assignin('base', 'cluster_data_saved', cluster_data_saved)

    legend('dCIRC (' + Genotype + ')')
    xlim([0 24])
    ylim([-1.005 1.005])
    dim = [0.61 .52 .3 .3];
    str = {['\color{black} \tau=',num2str(tau_DD,'%.1f'),', s=', num2str(s,'%.2f'), ', a=', num2str(a,'%.2f')]};
    annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',35,'EdgeColor','none');
%     title({['\color{black}CIRC (tau_{DD}=',num2str(tau_DD),', s=', num2str(s), ', a=', num2str(a),')']})
%     title(Genotype)
    set(gca,'fontsize',42);
    xticks(0:6:30)
    xticklabels({'CT0','CT6','CT12','CT18','CT24','CT30'})
    xlabel('internal time/hours','Interpreter','latex')
    ylabel('capacity to change $\tau_I$','Interpreter','latex')
    text(0.5,0.1,'compression','Color','r','FontSize',32,'Rotation',90);
    text(0.5,-0.8,'expansion','Color','b','FontSize',32,'Rotation',90);
    grid on
    grid minor
    set(fh,'Position',[0 0 1280 640]);
    saveas(gcf,'best_fitting_dCiRC_Curves/dCIRC-'+Genotype, 'epsc')
%     saveas(gcf,['dCIRC-Curves-Optimized/03-dCIRC-'+Genotype + '.pdf'])
    close all

end
%     hold off

function [ derivative ] = RK4_circadian_CIRC(y)
    theta_I = y(1,1);
    omega_I(i) = y(1,2);
    Normalized_theta_Z(i) = mod(zeitgeber_clock(i),2*pi);
    Normalized_theta_I(i) = mod(theta_I,2*pi);

    % Light function
    if abs(Normalized_theta_Z(i) - light_middle) < epsilon 
            light(i) = cz;
    else
        light(i) = 0;
    end

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
        domega_I(i) = light(i) * CIRC(i) + elastic_factor * (omega_DD_sim - omega_I(i));
        derivative = [dtheta_I, domega_I(i)]';
    end
end
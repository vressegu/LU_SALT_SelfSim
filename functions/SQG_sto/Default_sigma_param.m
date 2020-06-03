%% Default parameters for sigma
% Do not change

if sigma.sto
    % Band pass filter of the residual
    % (for sigma.type_spectrum = 'SelfSim_from_LS' only) 
    sigma.band_pass_filter = true;
    % Maybe better with sigma.band_pass_filter = false (slightly stronger v')
    
    % Estimation of the minimum possible energetic sigma wave number
    % (based on large-scale spectrum maximum)
    sigma.estim_k_LS = false;
    
    % Homogeneous dissipation associated with the spectrum slope
    sigma.assoc_diff = false;
    
    %% EOF parameters
    if strcmp(sigma.type_spectrum,'EOF') ...
            || strcmp(sigma.type_spectrum,'Euler_EOF')
        % Ratio between the Shanon resolution and filtering frequency used to
        % filter the heterogenous diffusion coefficient
        Smag.dealias_ratio_mask_LS = 1/8;
        
        % Nb day used to learn to EOFs
        sigma.nbDayLearn= 50;
        
        % Time period sampling for the learning data set
        sigma.Delta_T_on_Delta_t = 8;
        
        % Number of EOF (use all EOFs if set to inf)
        sigma.nb_EOF = 200; % ref
    end
    
    %% Possible time smoothing of sigma
    if strcmp(sigma.type_spectrum,'SelfSim_from_LS')        
        sigma.time_smooth.bool = false;
        
        % Correlation time
        sigma.time_smooth.tau = (64/resolution) * 24*3600 / 10 ;
    end
    
    %% Sigma heterogenity parameters
    
    % Heterogeneous energy flux epsilon
    % ( for sigma.type_spectrum = 'SelfSim_from_LS' only )
    sigma.hetero_energy_flux_v2 = sigma.hetero_energy_flux;
    
    % Smagorinsky-like control of dissipation
    sigma.Smag.bool = false;
    
    % Modulation by local V L (estimated from the velocity and from
    % thegradient of the velocity)
    sigma.hetero_modulation = false;
    
    % Modulation by local V^2
    sigma.hetero_modulation_V2 = false;
    
    if sigma.hetero_modulation | sigma.hetero_energy_flux ...
            | sigma.hetero_modulation_V2 | sigma.hetero_modulation_Smag
        % Ratio between the Shanon resolution and filtering frequency used to
        % filter the heterogenous diffusion coefficient
        Smag.dealias_ratio_mask_LS = 1; % default value
        
        % Compute mudulation from filtered (kappa_min) fields
        sigma.hetero_energy_flux_prefilter = true;
        
        % Filter noise modulations (1/4*kappa_min) fields
        sigma.hetero_energy_flux_postfilter = true;
    else
        % Compute mudulation from filtered (kappa_min) fields
        sigma.hetero_energy_flux_prefilter = false;
    end
    
    if sigma.hetero_energy_flux_v2
        if ~ sigma.hetero_energy_flux
            error('sigma.hetero_energy_flux_v2 is always associated with sigma.hetero_energy_flux');
        end
    end
    
    % Force sigma to be diveregence free
    sigma.proj_free_div = true;
    
    % Use a spatial derivation scheme for the herogeneous
    % disspation
    Smag.spatial_scheme = false;
    
    if ( (sigma.Smag.bool + sigma.hetero_modulation + ...
            sigma.hetero_energy_flux + sigma.hetero_modulation_V2 ...
            + sigma.hetero_modulation_Smag) > 1 ) ...
            || ( (sigma.Smag.bool + sigma.assoc_diff ) > 1 )
        error('These heterogeneous parametrizations cannot be combined');
    end
    
    %% For Smagorinsky-like diffusivity/viscosity or Hyper-viscosity,
    if sigma.Smag.bool || sigma.assoc_diff
        % Rate between the smallest wave number of the spatially-unresolved
        % (not simulated) component of sigma dBt and the largest wave
        % number of the simulation
        sigma.kappaMinUnresolved_on_kappaShanon = 1;
        
        % Rate between the largest wave number of the spatially-unresolved
        % (not simulated) component of sigma dBt and the largest wave
        % number of the simulation
        sigma.kappaMaxUnresolved_on_kappaShanon = 8;
        
    end
    
    if sigma.Smag.bool
        % Smagorinsky energy budget (dissipation epsilon)
        % without taking into account the noise intake
        sigma.Smag.epsi_without_noise = false;
        
        % Ratio between the Shanon resolution and filtering frequency used to
        % filter the heterogenous diffusion coefficient
        Smag.dealias_ratio_mask_LS = 1/8;
        
        sigma.Smag.kappamax_on_kappad = 0.5;
        
        % Heterogeneity of the noise
        sigma.Smag.SS_vel_homo = false;
    end
    
    %% Desactivate the noise
    sigma.no_noise = false;
    if sigma.no_noise
        warning('There is no noise here');
    end
    
    %% Variance tensor a_H
    switch sigma.type_spectrum
        case 'SelfSim_from_LS'
            sigma.k_c = 0;
        otherwise
            switch dynamics
                case 'SQG'
                    switch resolution
                        case 64
                            sigma.k_c = 1/(3e2) /4 % 1/(1200 meters)
                        case 128
                            sigma.k_c = 1/(3e2) % 1/(300 meters)
                    end
                case '2D'
                    error(...
                        'The turbulence 2D is not stable under the action of noise');
                otherwise
                    error('Unknown type of dynamics');
            end
            % a_H is calculated latter in the code using
            % a_H = 2 * f_0 / k_c^2
            % where f_0 is the Corilis frequency
    end
    
    %% Reference spectrum slope of sigma dBt
    % (used as an upper bound for the spectrum slope estimation)
    switch dynamics
        case 'SQG'
            sigma.slope_sigma = - 5/3;
        case '2D'
            sigma.slope_sigma = - 3;
        otherwise
            error('Unknown type of dynamics');
    end
    sigma.slope_sigma_ref = sigma.slope_sigma;
    if  sigma.sto & strcmp(sigma.type_spectrum,'BB')
        sigma.slope_sigma = 0;
    end
    
    %% Definition of the function to compute the variance tensor
    if sigma.sto
        eval(['sigma.fct_tr_a = @(m,k1,k2) fct_norm_tr_a_theo_' ...
            sigma.type_spectrum '(m,k1,k2);']);
    end
else
    % Compute mudulation from filtered (kappa_min) fields
    sigma.hetero_energy_flux_prefilter = false;
    
    sigma.hetero_energy_flux_v2 = false;
end

% Normalize before taking the power 1/3
sigma.hetero_energy_flux_averaging_after = true;

% To filter negative part of the energy flux
sigma.kappa_VLS_on_kappa_LS = 1/8;

% Maximum value considered
sigma.kappaLSforEspi_on_kappamin = 1;

% ADSD comparaison of Self Sim and EOF methods
sigma.plot_ADSD_SelfSim_EOF_bool = false;

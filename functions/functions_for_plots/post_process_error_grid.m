function error_vs_t = post_process_error_grid(stochastic_simulation,...
    type_data,resolution,resolution_HR,forcing,sigma,Lap_visco,HV,Smag,...
    N_ech,first_day,advection_duration,day_save)
% plot the same thing that fct_fft_advection_sto
%

last_day_plot_qq = floor(advection_duration/3600/24);


if nargin == 0
    init;
end
%% Main parameters to choose
% Type of dynamics
dynamics = 'SQG';
%dynamics = '2D';

plot_random_IC = false;
random_IC_large = false


if nargin == 0
    
    % Deterministic or random model
    stochastic_simulation = true;
    sigma.sto = stochastic_simulation;
    % Usual SQG model (stochastic_simulation=false)
    % or SQG_MU model (stochastic_simulation=true)
    
    if sigma.sto
        % Type of spectrum for sigma dBt
        % sigma.type_spectrum = 'Band_Pass_w_Slope' % as in GAFD part II
        % sigma.type_spectrum = 'Low_Pass_w_Slope';
        % Spectrum cst for k<km ans slope for k>km
        % sigma.type_spectrum = 'Low_Pass_streamFct_w_Slope';
        % Matern covariance for the streamfunction
        % spectrum = cst. * k2 .* ( 1 + (k/km)^2 )^slope )
        % ~ k2 for k<km ans slope for k>km
        % type_spectrum = 'BB';
        % type_spectrum = 'Bidouille';
        % sigma.type_spectrum = 'SelfSim_from_LS'
        sigma.type_spectrum = 'EOF'
        %  Sigma computed from self similarities from the large scales
        % sigma.type_spectrum = type_spectrum;
        
        % Homogeneous dissipation associated with the spectrum slope
        sigma.assoc_diff = false;
        
        % Smagorinsky-like control of dissipation
        sigma.Smag.bool = false;
        
        %     % Sigma computed from self similarities from the large scales
        %     sigma.SelfSim_from_LS.bool = true;
        
        %     if sigma.SelfSim_from_LS.bool
        %         % Sigma computed from a energy of absolute diffusivity spectrum
        %         % sigma.SelfSim_from_LS.spectrum = 'energy';
        %         sigma.SelfSim_from_LS.spectrum = 'abs_diff';
        %     end
        
        % if strcmp(sigma.type_spectrum,'SelfSim_from_LS')
        % Heterrogeenosu energy flux epsilon
        sigma.hetero_energy_flux = false;
        
        % Modulation by local V L (estimated from the velocity and from
        % thegradient of the velocity)
        sigma.hetero_modulation = false;
        
        % Modulation by local V^2
        sigma.hetero_modulation_V2 = false;
        
        %     %if strcmp(sigma.type_spectrum,'SelfSim_from_LS')
        %     if sigma.hetero_modulation & strcmp(sigma.type_spectrum,'SelfSim_from_LS')
        if sigma.hetero_modulation | sigma.hetero_energy_flux ...
                | sigma.hetero_modulation_V2 || strcmp(sigma.type_spectrum,'EOF')
            % Ratio between the Shanon resolution and filtering frequency used to
            % filter the heterogenous diffusion coefficient
            Smag.dealias_ratio_mask_LS = 1/8;
            % Smag.dealias_ratio_mask_LS = 1/4;
            
        end
        % end
        
        % Force sigma to be diveregence free
        sigma.proj_free_div = true;
        
        if ( (sigma.Smag.bool + sigma.hetero_modulation + ...
                sigma.hetero_energy_flux + sigma.hetero_modulation_V2 ) > 1 ) ...
                || ( (sigma.Smag.bool + sigma.assoc_diff ) > 1 )
            error('These parametrizations cannot be combined');
        end
        
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
        % For Smagorinsky-like diffusivity/viscosity or Hyper-viscosity,
        if sigma.Smag.bool
            % Smagorinsky energy budget (dissipation epsilon)
            % without taking into account the noise intake
            sigma.Smag.epsi_without_noise = false;
            
            % Ratio between the Shanon resolution and filtering frequency used to
            % filter the heterogenous diffusion coefficient
            Smag.dealias_ratio_mask_LS = 1;
            
            % sigma.Smag.kappamax_on_kappad = 0.5; % (better(?))
            % sigma.Smag.kappamax_on_kappad = 1 / 4;
            sigma.Smag.kappamax_on_kappad = 1 / ...
                sigma.kappaMaxUnresolved_on_kappaShanon;
            
            % Heterogeneity of the noise
            sigma.Smag.SS_vel_homo = false;
            
        end
        
        % Desactivate the noise
        sigma.no_noise = false;
        if sigma.no_noise
            warning('There is no noise here');
        end
    end
end

% Number of realizations in the ensemble
if nargin == 0
    N_ech=200;
end
% ( N_ech is automatically set to 1 in deterministic simulations )

if nargin == 0
    
    first_day = 0
    
    % Type of initial condtions
    type_data = 'Vortices'
    
    % Resolution
    resolution = 64
    
    % Resolution of the reference
    resolution_HR = 512;
    % resolution_HR = 1024;
    
    % The number of grid point is resolution^2
    % It has to be an even integer
    
    % Forcing
    
    % Forcing or not
    forcing = false;
end

% Type de forcing
% forcing_type = 'Kolmogorov';
forcing_type = 'Spring'

% Amplitude of the forcing
ampli_forcing = 10;

% Frequency of the forcing
freq_f = [3 2]


%% Deterministic subgrid tensor

if nargin == 0
    % Viscosity
    Lap_visco.bool = false;
    
    % Hyper-viscosity
    HV.bool = true;
    
    if HV.bool
        % HV.order=4;
        HV.order=8;
    end
    
    % Smagorinsky-like diffusivity/viscosity or Hyper-viscosity
    Smag.bool = false;
    
    % For Smagorinsky-like diffusivity/viscosity or Hyper-viscosity,
    if Smag.bool
        if Lap_visco.bool
            
            % Use a spatial derivation scheme for the herogeneous
            % disspation
            Smag.spatial_scheme = false;
            
            % Ratio between the Shanon resolution and filtering frequency used to
            % filter the heterogenous diffusion coefficient
            Smag.dealias_ratio_mask_LS = 1/8;
            warning('Redondant argument that for heterogeneous small-scale velocity')
            
            % Ratio between the Shanon resolution cut-off ( = pi / sqrt( dx*dy) )
            % and the targeted diffusion scale
            Smag.kappamax_on_kappad = 0.5;
            
            % Factor in front of the additional constant dissipation
            % Set to 0 for no additional constant dissipation
            Smag.weight_cst_dissip = 0;
        elseif HV.bool
            % Ratio between the Shanon resolution and filtering frequency used to
            % filter the heterogenous diffusion coefficient
            Smag.dealias_ratio_mask_LS = 1;
            warning('Redondant argument that for heterogeneous small-scale velocity')
            
            % Ratio between the Shanon resolution cut-off ( = pi / sqrt( dx*dy) )
            % and the targeted diffusion scale
            Smag.kappamax_on_kappad = 1.1;% still small oscillations or just pixels?
            
            % Factor in front of the additional constant dissipation
            % Set to 0 for no additional constant dissipation
            
            Smag.weight_cst_dissip = 1/1;
        end
    else
        Smag.kappamax_on_kappad = 0;
    end
end

%% Optional parameters
% Choose to plot
if nargin > 0
    plots_bool = false;
else
    plots_bool = true;
end

% Compute velocity covariance and absolute diffusivity
cov_and_abs_diff = false;

% Choose to plot one-point one-time moments each day
plot_moments = false;

% Choose to plot the dissipation by scale
plot_epsilon_k = true;
if sigma.sto & sigma.hetero_energy_flux
    plot_epsilon_k = true;
end

% Plot dissipations terms
plot_dissip = true;

if nargin == 0
    % Begin simulation from a precomputed field?
    use_save = false;
    % In this case, which day should be used as initialisation
    if use_save
        day_save = 7;
    else
        day_save = 0;
    end
end

dealias_method = 'exp';
% [WIP] Method for mandatory de-aliasing of non-linear terms in
% pseudospectral codes (advection and non-homogeneous stochastic diffusion)
% - 'lowpass': same as in SQGMU 1;
% - '2/3': the classical 2/3 rule (NB: produces Gibb's oscillations);
% - 'exp': high-order exponential filter (Constantin et al., J. Sci.
%   Comput. (2012)).

% Boundaries conditions
dirichlet = false;

if nargin == 0
    % Variance tensor a_H
    if stochastic_simulation
        if strcmp(sigma.type_spectrum , 'SelfSim_from_LS')
            pre_estim_slope=1e-1;
            pre_5 = 5e-2;
            sigma.kappamin_on_kappamax = ...
                (log(1-pre_5)/log(pre_estim_slope))^(2/HV.order);
            sigma.kappamin_on_kappamax_estim_slope = ...
                (log(1-pre_estim_slope)/log(pre_estim_slope))...
                ^(2/HV.order);
            
            sigma.kappaLS_on_kappamax = 1/8;
        else
            switch resolution
                case  128
                    sigma.kappamin_on_kappamax = 1/2;
                case 64
                    sigma.kappamin_on_kappamax = 1/3;
                otherwise
                    error('unknown');
            end
            
            sigma.kappaLS_on_kappamax = 1/8;
        end
    else
        % If the simulation is deterministic, a_H = 0 and only one simulation
        % is performed
        sigma.k_c = inf; % And then a_H = 0
        N_ech=1;
        plot_moments = false;
    end
    
    % Spectrum slope of sigma dBt
    if sigma.sto
        switch dynamics
            case 'SQG'
                sigma.slope_sigma = - 5/3;
            case '2D'
                sigma.slope_sigma = - 3;
            otherwise
                error('Unknown type of dynamics');
        end
    end
    if  sigma.sto & strcmp(sigma.type_spectrum,'BB')
        sigma.slope_sigma = 0;
    end
    
    % Rate between the smallest and the largest wave number of sigma dBt
    if sigma.sto
        if strcmp(sigma.type_spectrum , 'SelfSim_from_LS')
            sigma.kappamin_on_kappamax = 1/2;
            
            sigma.kappaLS_on_kappamax = 1/8;
        else
            sigma.kappamin_on_kappamax = 1/2;
            
            sigma.kappaLS_on_kappamax = 1/8;
        end
        
        % Rate between the largest wave number of sigma dBt and the largest wave
        % number of the simulation
        sigma.kappamax_on_kappaShanon = 1;
    end
end

% Spectrum slope of the initial condition (if type_data = 'Spectrum' )
switch dynamics
    case 'SQG'
        slope_b_ini = - 5/3;
    case '2D'
        slope_b_ini = - 3;
    otherwise
        error('Unknown type of dynamics');
end

% Physical parameters
model = fct_physical_param(dynamics);

% Gather parameters in the structure model
model.sigma = sigma;
if sigma.sto
    eval(['model.sigma.fct_tr_a = @(m,k1,k2) fct_norm_tr_a_theo_' ...
        model.sigma.type_spectrum '(m,k1,k2);']);
end
if strcmp(type_data,'Spectrum')
    model.slope_b_ini = slope_b_ini;
end
model.dynamics=dynamics;
model.type_data=type_data;
model.advection.N_ech=N_ech;
model.advection.plot_epsilon_k = plot_epsilon_k;
model.advection.plot_dissip = plot_dissip;
model.advection.plot_moments = plot_moments;
model.advection.forcing.bool = forcing;
model.advection.forcing.ampli_forcing = ampli_forcing;
model.advection.forcing.freq_f = freq_f;
model.advection.forcing.forcing_type = forcing_type;
model.advection.HV = HV;
model.advection.cov_and_abs_diff = cov_and_abs_diff;
model.advection.Lap_visco = Lap_visco;
model.advection.Smag = Smag;
if nargin == 0
    model.advection.use_save = use_save;
    model.advection.day_save = day_save;
end
model.grid.dealias_method = dealias_method; %de-aliasing method
model.plots = plots_bool;


%% Generating initial buoyancy
[~,model_HR] = fct_buoyancy_init(model,resolution_HR);
[~,model] = fct_buoyancy_init(model,resolution);

%% Set up
plot_modes = true;
plot_error = false;
plot_high_order_moments = false;
nb_modes = 200;

%% Folder with reference
clear subgrid_details
% if model.advection.HV.bool
% add_subgrid_deter = '_HV';
add_subgrid_deter = ['_HV' '_' fct_num2str(model.advection.HV.order/2)];
% if isinf(model.sigma.k_c) % Deterministic case
model_HR.folder.folder_simu = [ 'images/usual_' model.dynamics ...
    add_subgrid_deter '/' model.type_data ];

if model.advection.forcing.bool
    model_HR.folder.folder_simu = [ model_HR.folder.folder_simu ...
        '_forced_turb_' model_HR.advection.forcing.forcing_type];
else
    model_HR.folder.folder_simu = [ model_HR.folder.folder_simu ...
        '_free_turb' ];
end
model_HR.folder.folder_simu = [ model_HR.folder.folder_simu ...
    '/' num2str(model_HR.grid.MX(1)) 'x' num2str(model_HR.grid.MX(2)) ];
model_HR.folder.folder_simu = [ model_HR.folder.folder_simu ...
    '/Low_Pass_fitlered_version_' ...
    num2str(model.grid.MX(1)) 'x' num2str(model.grid.MX(2)) ];

%% Folder with deterministic model with random initial condition
clear subgrid_details
model_randomIC = model;
model_randomIC.sigma.sto = false;
add_subgrid_deter = ['_HV' '_' fct_num2str(model.advection.HV.order/2)];
if plot_random_IC
    model_randomIC.folder.folder_simu = [ 'images/usual_' model.dynamics ...
        '_randomIC' add_subgrid_deter '/' model.type_data ];
else
    model_randomIC.folder.folder_simu = [ 'images/usual_' model.dynamics ...
        add_subgrid_deter '/' model.type_data ];
end
if model.advection.forcing.bool
    model_randomIC.folder.folder_simu = [ model_randomIC.folder.folder_simu ...
        '_forced_turb_' model_randomIC.advection.forcing.forcing_type];
else
    model_randomIC.folder.folder_simu = [ model_randomIC.folder.folder_simu ...
        '_free_turb' ];
end
model_randomIC.folder.folder_simu = [ model_randomIC.folder.folder_simu ...
    '/' num2str(model_randomIC.grid.MX(1)) 'x' num2str(model_randomIC.grid.MX(2)) ];
if plot_random_IC
    if random_IC_large
        model_randomIC.folder.folder_simu = [ model_randomIC.folder.folder_simu ...
            '/large_IC_perturb' ];
    else
        model_randomIC.folder.folder_simu = [ model_randomIC.folder.folder_simu ...
            '/small_IC_perturb' ];
    end
end



%% Folder with deterministic model with random initial condition
clear subgrid_details
model_deter = model;
model_deter.sigma.sto = false;
% if model.advection.HV.bool
% add_subgrid_deter = '_HV';
add_subgrid_deter = ['_HV' '_' fct_num2str(model.advection.HV.order/2)];
model_deter.folder.folder_simu = [ 'images/usual_' model.dynamics ...
    add_subgrid_deter '/' model.type_data ];
if model.advection.forcing.bool
    model_deter.folder.folder_simu = [ model_deter.folder.folder_simu ...
        '_forced_turb_' model_deter.advection.forcing.forcing_type];
else
    model_deter.folder.folder_simu = [ model_deter.folder.folder_simu ...
        '_free_turb' ];
end
model_deter.folder.folder_simu = [ model_deter.folder.folder_simu ...
    '/' num2str(model_deter.grid.MX(1)) 'x' num2str(model_deter.grid.MX(2)) ];

%% Folder to save plots and files
clear subgrid_details
if model.advection.HV.bool
    add_subgrid_deter = ['_HV' '_' fct_num2str(model.advection.HV.order/2)];
elseif model.advection.Lap_visco.bool
    add_subgrid_deter = '_Lap_visco';
else
    add_subgrid_deter = '_no_deter_subgrid';
end
if model.sigma.sto & model.sigma.assoc_diff
    add_subgrid_deter = [add_subgrid_deter '_assoc_diff'];
end
if ( ( model.advection.HV.bool | model.advection.Lap_visco.bool) & ...
        model.advection.Smag.bool ) | ...
        (model.sigma.sto & model.sigma.Smag.bool )
    add_subgrid_deter = [add_subgrid_deter '_Smag'];
    if model.sigma.sto & model.sigma.Smag.bool & ...
            model.sigma.Smag.epsi_without_noise
        add_subgrid_deter = [add_subgrid_deter '_epsi_without_noise'];
    end
elseif model.sigma.sto & model.sigma.hetero_modulation
    add_subgrid_deter = [add_subgrid_deter '_hetero_modulation'];  
    if  isfield(model.sigma,'hetero_energy_flux_prefilter') && ...
            model.sigma.hetero_energy_flux_prefilter
        add_subgrid_deter = [add_subgrid_deter ...
            '_prefilter'];        
    end
    if  isfield(model.sigma,'hetero_energy_flux_postfilter') && ...
            model.sigma.hetero_energy_flux_postfilter
        add_subgrid_deter = [add_subgrid_deter ...
            '_postfilter'];
    end
elseif model.sigma.sto & model.sigma.hetero_modulation_V2
    add_subgrid_deter = [add_subgrid_deter '_hetero_modulation_V2'];  
    if  isfield(model.sigma,'hetero_energy_flux_prefilter') && ...
            model.sigma.hetero_energy_flux_prefilter
        add_subgrid_deter = [add_subgrid_deter ...
            '_prefilter'];        
    end   
    if  isfield(model.sigma,'hetero_energy_flux_postfilter') && ...
            model.sigma.hetero_energy_flux_postfilter
        add_subgrid_deter = [add_subgrid_deter ...
            '_postfilter'];
    end
elseif model.sigma.sto & model.sigma.hetero_energy_flux
    add_subgrid_deter = [add_subgrid_deter '_hetero_energy_flux'];
    if model.sigma.hetero_energy_flux_v2
        add_subgrid_deter = [add_subgrid_deter '_v2'];
    end
    if model.sigma.hetero_energy_flux_averaging_after
        add_subgrid_deter = [add_subgrid_deter '_1_3rd_before_norm'];
    end
    if isfield(model.sigma,'kappa_VLS_on_kappa_LS')
        add_subgrid_deter = [add_subgrid_deter ...
            '_kappa_LS_on_kappa_VLS_' ...
            num2str(1/model.sigma.kappa_VLS_on_kappa_LS)];
    end
    if isfield(model.sigma,'kappaLSforEspi_on_kappamin')
        add_subgrid_deter = [add_subgrid_deter ...
            '_kappamin_on_kappaLSforEspi__' ...
            num2str(1/model.sigma.kappaLSforEspi_on_kappamin)];
    end
    if  isfield(model.sigma,'hetero_energy_flux_prefilter') && ...
            model.sigma.hetero_energy_flux_prefilter
        add_subgrid_deter = [add_subgrid_deter ...
            '_prefilter'];        
    end
    if  isfield(model.sigma,'hetero_energy_flux_postfilter') && ...
            model.sigma.hetero_energy_flux_postfilter
        add_subgrid_deter = [add_subgrid_deter ...
            '_postfilter'];
    end
elseif model.sigma.sto & model.sigma.hetero_modulation_Smag
    add_subgrid_deter = [add_subgrid_deter '_hetero_modulation_Smag'];    
    if  isfield(model.sigma,'hetero_energy_flux_prefilter') && ...
            model.sigma.hetero_energy_flux_prefilter
        add_subgrid_deter = [add_subgrid_deter ...
            '_prefilter'];        
    end
    if  isfield(model.sigma,'hetero_energy_flux_postfilter') && ...
            model.sigma.hetero_energy_flux_postfilter
        add_subgrid_deter = [add_subgrid_deter ...
            '_postfilter'];
    end
end
if model.sigma.sto & model.sigma.no_noise
    add_subgrid_deter = [add_subgrid_deter '_no_noise'];
end

if ~ model.sigma.sto % Deterministic case
    model.folder.folder_simu = [ 'images/usual_' model.dynamics ...
        add_subgrid_deter '/' model.type_data ];
else % Stochastic case
    model.folder.folder_simu = [ 'images/' model.dynamics ...
        '_MU' add_subgrid_deter '/' ...
        'type_spectrum_sigma_' model.sigma.type_spectrum '/' ...
        model.type_data ];
end
if model.advection.forcing.bool
    model.folder.folder_simu = [ model.folder.folder_simu ...
        '_forced_turb_' model.advection.forcing.forcing_type];
else
    model.folder.folder_simu = [ model.folder.folder_simu ...
        '_free_turb' ];
end
model.folder.folder_simu = [ model.folder.folder_simu ...
    '/' num2str(model.grid.MX(1)) 'x' num2str(model.grid.MX(2)) ];
if ( ( model.advection.HV.bool | model.advection.Lap_visco.bool) & ...
        model.advection.Smag.bool)
    subgrid_details = ['kappamax_on_kappad_' ...
        fct_num2str(model.advection.Smag.kappamax_on_kappad) ...
        '_dealias_ratio_mask_LS_' ...
        fct_num2str(model.advection.Smag.dealias_ratio_mask_LS)];
    if model.advection.Smag.spatial_scheme
        subgrid_details = [ subgrid_details '_spatial_scheme'];
    end
    model.folder.folder_simu = [ model.folder.folder_simu ...
        '/' subgrid_details ];
end
if model.sigma.sto
    if model.sigma.Smag.bool
        subgrid_details = ['kappamax_on_kappad_' ...
            fct_num2str(model.sigma.Smag.kappamax_on_kappad) ...
            '_dealias_ratio_mask_LS_' ...
            fct_num2str(model.advection.Smag.dealias_ratio_mask_LS)];
        if model.sigma.Smag.SS_vel_homo
            subgrid_details = [ subgrid_details '_SS_vel_homo'];
        elseif  model.sigma.proj_free_div
            subgrid_details = [ subgrid_details '_proj_free_div'];
        end
        if model.advection.Smag.spatial_scheme
            subgrid_details = [ subgrid_details '_spatial_scheme'];
        end
    elseif ( model.sigma.hetero_modulation ...
            |  model.sigma.hetero_modulation_V2 ...
            |  model.sigma.hetero_modulation_Smag  ...
            | model.sigma.hetero_energy_flux )
        subgrid_details = ['dealias_ratio_mask_LS_' ...
            fct_num2str(model.advection.Smag.dealias_ratio_mask_LS)];
        if  model.sigma.proj_free_div
            subgrid_details = [ subgrid_details '_proj_free_div'];
        end
    end
    if ~ ( exist('subgrid_details','var')==1)
        subgrid_details = [];
    end
    if ~ strcmp(model.sigma.type_spectrum,'EOF')
        subgrid_details = [ subgrid_details ...
            '_kappamin_on_kappamax_' ....
            fct_num2str(model.sigma.kappamin_on_kappamax) ];
        if strcmp(model.sigma.type_spectrum,'Band_Pass_w_Slope')
            subgrid_details = [ subgrid_details ...
                '_on_kc_' ....
                fct_num2str(1/model.sigma.k_c) ];
        elseif strcmp(model.sigma.type_spectrum,'SelfSim_from_LS') ...
                if model.sigma.estim_k_LS
                subgrid_details = [ subgrid_details ...
                    '_estim_k_LS'];
                end
                if model.sigma.time_smooth.bool
                    subgrid_details = [ subgrid_details ...
                        '_time_smooth_'...
                        num2str(24*3600/model.sigma.time_smooth.tau)];
                end
        end
    else
        subgrid_details = [ subgrid_details ...
            '_nbDayLearn_' ...
            fct_num2str(model.sigma.nbDayLearn) ...
            '_Delta_T_on_Delta_t_' ...
            fct_num2str(model.sigma.Delta_T_on_Delta_t) ...;
            '_nb_EOF_' ...
            fct_num2str(model.sigma.nb_EOF)];
    end
    if model.advection.N_ech>1
        subgrid_details = [ subgrid_details ...
            '_N_ech_' ....
            fct_num2str(model.advection.N_ech) ];
    end
    model.folder.folder_simu = [ model.folder.folder_simu ...
        '/' subgrid_details ];
end

% Create the folders
fct_create_folder_plots(model,random_IC_large,plot_random_IC)

% Colormap
load('BuYlRd.mat');
model.folder.colormap = BuYlRd; clear BuYlRd

% Version of matlab
vers = version;
year = str2double(vers(end-5:end-2));
subvers = vers(end-1);
model.folder.colormap_freeze = ...
    (  year < 2014 || ( year == 2014 && strcmp(subvers,'a') ) );

%% Comparaison of LU parametrisation
clear subgrid_details
if model.sigma.sto && strcmp(model.sigma.type_spectrum , 'EOF')
    model_SelfSim = model;
    model_SelfSim.sigma.type_spectrum = 'SelfSim_from_LS';
    model_SelfSim.sigma.estim_k_LS = false;
    model_SelfSim.sigma.time_smooth.bool=false;
    model_SelfSim.sigma.time_smooth.tau = 24*3600 / 10;
    
    if model_SelfSim.advection.HV.bool
        add_subgrid_deter = ['_HV' '_' fct_num2str(model_SelfSim.advection.HV.order/2)];
    elseif model_SelfSim.advection.Lap_visco.bool
        add_subgrid_deter = '_Lap_visco';
    else
        add_subgrid_deter = '_no_deter_subgrid';
    end
    if model_SelfSim.sigma.sto & model_SelfSim.sigma.assoc_diff
        add_subgrid_deter = [add_subgrid_deter '_assoc_diff'];
    end
    if ( ( model_SelfSim.advection.HV.bool | model_SelfSim.advection.Lap_visco.bool) & ...
            model_SelfSim.advection.Smag.bool ) | ...
            (model_SelfSim.sigma.sto & model_SelfSim.sigma.Smag.bool )
        add_subgrid_deter = [add_subgrid_deter '_Smag'];
        if model_SelfSim.sigma.sto & model_SelfSim.sigma.Smag.bool & ...
                model_SelfSim.sigma.Smag.epsi_without_noise
            add_subgrid_deter = [add_subgrid_deter '_epsi_without_noise'];
        end
    elseif model_SelfSim.sigma.sto & model_SelfSim.sigma.hetero_modulation
        add_subgrid_deter = [add_subgrid_deter '_hetero_modulation'];
    elseif model_SelfSim.sigma.sto & model_SelfSim.sigma.hetero_modulation_V2
        add_subgrid_deter = [add_subgrid_deter '_hetero_modulation_V2'];
    elseif model_SelfSim.sigma.sto & model_SelfSim.sigma.hetero_energy_flux
        add_subgrid_deter = [add_subgrid_deter '_hetero_energy_flux'];
    end
    if model_SelfSim.sigma.sto & model_SelfSim.sigma.no_noise
        add_subgrid_deter = [add_subgrid_deter '_no_noise'];
    end
    if ~ model_SelfSim.sigma.sto % Deterministic case
        model_SelfSim.folder.folder_simu = [ 'images/usual_' model_SelfSim.dynamics ...
            add_subgrid_deter '/' model_SelfSim.type_data ];
    else % Stochastic case
        model_SelfSim.folder.folder_simu = [ 'images/' model_SelfSim.dynamics ...
            '_MU' add_subgrid_deter '/' ...
            'type_spectrum_sigma_' model_SelfSim.sigma.type_spectrum '/' ...
            model_SelfSim.type_data ];
    end
    if model_SelfSim.advection.forcing.bool
        model_SelfSim.folder.folder_simu = [ model_SelfSim.folder.folder_simu ...
            '_forced_turb_' model_SelfSim.advection.forcing.forcing_type];
    else
        model_SelfSim.folder.folder_simu = [ model_SelfSim.folder.folder_simu ...
            '_free_turb' ];
    end
    model_SelfSim.folder.folder_simu = [ model_SelfSim.folder.folder_simu ...
        '/' num2str(model_SelfSim.grid.MX(1)) 'x' num2str(model_SelfSim.grid.MX(2)) ];
    if ( ( model_SelfSim.advection.HV.bool | model_SelfSim.advection.Lap_visco.bool) & ...
            model_SelfSim.advection.Smag.bool)
        subgrid_details = ['kappamax_on_kappad_' ...
            fct_num2str(model_SelfSim.advection.Smag.kappamax_on_kappad) ...
            '_dealias_ratio_mask_LS_' ...
            fct_num2str(model_SelfSim.advection.Smag.dealias_ratio_mask_LS)];
        if model_SelfSim.advection.Smag.spatial_scheme
            subgrid_details = [ subgrid_details '_spatial_scheme'];
        end
        model_SelfSim.folder.folder_simu = [ model_SelfSim.folder.folder_simu ...
            '/' subgrid_details ];
    end
    if model_SelfSim.sigma.sto
        if model_SelfSim.sigma.Smag.bool
            subgrid_details = ['kappamax_on_kappad_' ...
                fct_num2str(model_SelfSim.sigma.Smag.kappamax_on_kappad) ...
                '_dealias_ratio_mask_LS_' ...
                fct_num2str(model_SelfSim.advection.Smag.dealias_ratio_mask_LS)];
            if model_SelfSim.sigma.Smag.SS_vel_homo
                subgrid_details = [ subgrid_details '_SS_vel_homo'];
            elseif  model_SelfSim.sigma.proj_free_div
                subgrid_details = [ subgrid_details '_proj_free_div'];
            end
            if model_SelfSim.advection.Smag.spatial_scheme
                subgrid_details = [ subgrid_details '_spatial_scheme'];
            end
        elseif ( model_SelfSim.sigma.hetero_modulation |  model_SelfSim.sigma.hetero_modulation_V2)
            subgrid_details = ['dealias_ratio_mask_LS_' ...
                fct_num2str(model_SelfSim.advection.Smag.dealias_ratio_mask_LS)];
            if  model_SelfSim.sigma.proj_free_div
                subgrid_details = [ subgrid_details '_proj_free_div'];
            end
        end
        if ~ ( exist('subgrid_details','var')==1)
            subgrid_details = [];
        end
        if ~ strcmp(model_SelfSim.sigma.type_spectrum,'EOF')
            subgrid_details = [ subgrid_details ...
                '_kappamin_on_kappamax_' ....
                fct_num2str(model_SelfSim.sigma.kappamin_on_kappamax) ];
            if strcmp(model_SelfSim.sigma.type_spectrum,'Band_Pass_w_Slope')
                subgrid_details = [ subgrid_details ...
                    '_on_kc_' ....
                    fct_num2str(1/model_SelfSim.sigma.k_c) ];
            elseif strcmp(model_SelfSim.sigma.type_spectrum,'SelfSim_from_LS')
                if model_SelfSim.sigma.estim_k_LS
                    subgrid_details = [ subgrid_details ...
                        '_estim_k_LS'];
                end
                if model_SelfSim.sigma.time_smooth.bool
                    subgrid_details = [ subgrid_details ...
                        '_time_smooth_'...
                        num2str(24*3600/model_SelfSim.sigma.time_smooth.tau)];
                end
            end
        else
            subgrid_details = [ subgrid_details ...
                '_nbDayLearn_' ...
                fct_num2str(model_SelfSim.sigma.nbDayLearn) ...
                '_Delta_T_on_Delta_t_' ...
                fct_num2str(model_SelfSim.sigma.Delta_T_on_Delta_t) ...;
                '_nb_EOF_' ...
                fct_num2str(model_SelfSim.sigma.nb_EOF)];
        end
        if model_SelfSim.advection.N_ech>1
            subgrid_details = [ subgrid_details ...
                '_N_ech_' ....
                fct_num2str(model_SelfSim.advection.N_ech) ];
        end
        model_SelfSim.folder.folder_simu = [ model_SelfSim.folder.folder_simu ...
            '/' subgrid_details ];
    end
end

%% Grid

% Spatial grid
model.grid.x = model.grid.dX(1)*(0:model.grid.MX(1)-1);
model.grid.y = model.grid.dX(2)*(0:model.grid.MX(2)-1);

dkxdky = (2*pi)^2 / (prod(model.grid.MX.*model.grid.dX));
dkappa = sqrt(dkxdky);

% Grid in Fourier space
model = init_grid_k (model);

% Ensemble size
N_ech=model.advection.N_ech;

%%

t_last_plot = -inf;
day_last_plot = - inf;
model.advection.step='finite_variation';

% t_ini=1;

folder_ref = model.folder.folder_simu;
name_file = [model.folder.folder_simu '/files/' num2str(day_save) '.mat'];
load(name_file)
model.folder.folder_simu = folder_ref;
if ~isfield(model.sigma,'sto')
    model.sigma.sto = ~isinf( model.sigma.k_c);
end

dt=model.advection.dt_adv;

F_save = [];
F_save2 = [];
bt1_HR_vect = [];
bt1_LR_vect = [];

error_vs_t = [];
error_vs_t_SelfSim = [];
v_day = [];

t_ini = first_day*24*3600/dt;
trigger = false;
dt_loop = dt;
%t_ini=1700000

t_loop = t_ini - 1;
while t_loop*dt_loop <= advection_duration
    t_loop = t_loop +1;
    day_num = (floor(t_loop*dt_loop/24/3600));
    
    if day_num > day_last_plot
        day_num = (floor(t_loop*dt_loop/24/3600));
        day = num2str(day_num);
        day
        day_last_plot = day_num;
        
        time =t_loop*dt_loop;
        
        model.advection.plot_modes = plot_modes;
        model.advection.nb_modes = nb_modes;
        t_last_plot=t_loop;
        id_part=1;
        
        width=1.2e3;
        height=0.5e3;
        if strcmp(model.type_data,'klein')
            width=width/2;
            r_c_ax = 0.5;
        else
            r_c_ax =1/1.5;
        end
        X0=[0 1];
        
        %% Specific day
        %         warning('the day is specified manually');
        %         day = day_choose;
        %         day = num2str(day);
        
        %% Load meth with random IC
        clear fft_b;
        if (model.advection.N_ech>1)
            %         if plot_random_IC & (model.advection.N_ech>1)
            name_file_randomIC = [model_randomIC.folder.folder_simu ...
                '/files/' num2str(day) '.mat'];
            load(name_file_randomIC,'fft_T_adv_part','fft_buoy_part','fft_b');
            if (exist('fft_b','var')==1)
            elseif (exist('fft_buoy_part','var')==1)
                fft_b = fft_buoy_part;
            elseif (exist('fft_T_adv_part','var')==1)
                fft_b = fft_T_adv_part;
            else
                error('Cannot find buoyancy field')
            end
            fft_b_classic = fft_b;clear fft_b;
        end
        
        %% Load determinsitic meth without random IC
        clear fft_b;
        name_file_deter = [model_deter.folder.folder_simu ...
            '/files/' num2str(day) '.mat'];
        load(name_file_deter,'fft_T_adv_part','fft_buoy_part','fft_b');
        if (exist('fft_b','var')==1)
        elseif (exist('fft_buoy_part','var')==1)
            fft_b = fft_buoy_part;
        elseif (exist('fft_T_adv_part','var')==1)
            fft_b = fft_T_adv_part;
        else
            error('Cannot find buoyancy field')
        end
        fft_b_deter = fft_b;clear fft_b;
        model_deter.grid.x=model_deter.grid.x_ref;
        model_deter.grid.y=model_deter.grid.y_ref;
        model_deter.folder.colormap=model.folder.colormap;
        
        %% Load selfSim
        if model.sigma.sto && strcmp(model.sigma.type_spectrum , 'EOF')
            model_SelfSim_ref = model_SelfSim;
            name_file_SelfSim = ...
                [model_SelfSim.folder.folder_simu '/files/' num2str(day) '.mat'];
            clear fft_b_SelfSim fft_b fft_buoy_part;
            if exist( name_file_SelfSim,'file')==2
                load(name_file_SelfSim,'fft_T_adv_part','fft_buoy_part','fft_b');
                if ~(exist('fft_b','var')==1)
                    fft_b_SelfSim = fft_buoy_part; clear fft_buoy_part
                else
                    fft_b_SelfSim = fft_b; clear fft_b
                end
            else
                warning('Cannot find the following file');
                fprintf([ name_file ' \n']);
            end
            
            
        end
        
        %% Load
            
        model_ref = model;
        name_file = [model.folder.folder_simu '/files/' num2str(day) '.mat'];
        clear fft_b fft_buoy_part;

        if exist( name_file,'file')==2
            if eval(day) == first_day
                model_ref_post_process = model;
            end
            load(name_file);
            if eval(day) == first_day
                model = model_ref_post_process;
            end
            model.folder.folder_simu = folder_ref;
            if ~(exist('fft_b','var')==1)
                if (exist('fft_buoy_part','var')==1)
                    fft_b = fft_buoy_part;
                    
                elseif (exist('fft_buoy_part_ref','var')==1)
                    fft_b = fft_buoy_part_ref; clear fft_buoy_part_ref
                else
                    error('Cannot find buoyancy field')
                end
            end
        else
            warning('Cannot find the following file');
            fprintf([ name_file ' \n']);
            return
        end
        
        %% Load reference
        name_file_HR = [model_HR.folder.folder_simu '/files/' num2str(day) '.mat'];
        load(name_file_HR,'fft_buoy_part_ref','spectrum_ref');
        model_HR.grid.x=model_HR.grid.x_ref;
        model_HR.grid.y=model_HR.grid.y_ref;
        model_HR.folder.colormap=model.folder.colormap;
        
        fprintf([ num2str(time/(24*3600)) ' days of advection \n'])
        
        % Colormap
        load('BuYlRd.mat');
        model.folder.colormap = BuYlRd; clear BuYlRd
        
        % Version of matlab
        vers = version;
        year = str2double(vers(end-5:end-2));
        subvers = vers(end-1);
        model.folder.colormap_freeze = ...
            (  year < 2014 || ( year == 2014 && strcmp(subvers,'a') ) );
        %%
        
        if model.advection.N_ech > 1
            model.advection.plot_moments = true;
            fct_plot_post_process(model,fft_b,day);
            pause(0.1);
        end
        
        fct_plot_post_process(model_deter,fft_b_deter,day);
        
        model.advection.plot_moments = false;
        [spectrum,name_plot] = fct_plot_post_process(model,fft_b,day);
        
        model_mean=model;
        model_mean.folder.folder_simu = [ model_mean.folder.folder_simu ...
            '/mean'];
        fct_create_folder_plots(model_mean)
        fct_plot_post_process(model_mean,mean(fft_b,4),day);
        
        if model.sigma.sto && strcmp(model.sigma.type_spectrum , 'EOF')
            color_SelfSim = [1 0 1];
            fct_spectrum_multi(model,fft_b_SelfSim,color_SelfSim);
            fct_spectrum_multi(model,fft_b_deter,'m');
            fct_spectrum_multi(model,fft_buoy_part_ref,'r');
            legend('-5/3',...
                ['LU EOF' num2str(resolution)],...
                ['LU Self.Sim.' num2str(resolution)],...
                ['Deter. ' num2str(resolution)],...
                ['Deter. ' num2str(resolution_HR')]);
        elseif model.sigma.sto
            fct_spectrum_multi(model,fft_b_deter,'m');
            fct_spectrum_multi(model,fft_buoy_part_ref,'r');
            legend('-5/3',['LU ' num2str(resolution)],...
                ['Deter. ' num2str(resolution)],...
                ['Deter. ' num2str(resolution_HR')]);
        else
            fct_spectrum_multi(model,fft_buoy_part_ref,'r');
            legend('-5/3',...
                ['Deter. ' num2str(resolution)],...
                ['Deter. ' num2str(resolution_HR')]);
        end
        eval( ['print -depsc ' model.folder.folder_simu '/Spectrum/' day '.eps']);
        
        %% Distance to reference
        % Spectrum discrepancy
        temp(1,1) = ...
            dkappa * sum(abs(spectrum(:) - spectrum_ref(:))) ;
        % Error
        error2 = 1/prod(model.grid.MX)^2 * ...
            sum(sum(abs(bsxfun(@plus,fft_b, ...
            - fft_buoy_part_ref)).^2 , 1) ,2);
        % Bias
        bias2 = 1/prod(model.grid.MX)^2 * ...
            sum(sum(abs(bsxfun(@plus,mean(fft_b ,4),...
            - fft_buoy_part_ref)).^2 , 1) ,2);
        temp(1,2) = sqrt(bias2);
        % RMSE
        temp(1,3) = sqrt(mean(error2 ,4));
        % min distance
        temp(1,4) = min(sqrt(error2),[],4);
        % Concatenation
        error_vs_t = [ error_vs_t ; temp ];
        v_day = [ v_day eval(day)];
        
        figure111=figure(111);
        close(figure111)
        figure111=figure(111);
        if model.sigma.sto && strcmp(model.sigma.type_spectrum , 'EOF')
            % Spectrum discrepancy
            temp(1,1) = nan ;
            % Error
            error2 = 1/prod(model.grid.MX)^2 * ...
                sum(sum(abs(bsxfun(@plus,fft_b_SelfSim, ...
                - fft_buoy_part_ref)).^2 , 1) ,2);
            % Bias
            bias2 = 1/prod(model.grid.MX)^2 * ...
                sum(sum(abs(bsxfun(@plus,mean(fft_b_SelfSim ,4),...
                - fft_buoy_part_ref)).^2 , 1) ,2);
            temp(1,2) = sqrt(bias2);
            % RMSE
            temp(1,3) = sqrt(mean(error2 ,4));
            % min distance
            temp(1,4) = min(sqrt(error2),[],4);
            % Concatenation
            error_vs_t_SelfSim = [ error_vs_t_SelfSim ; temp ];
            hold on;
            %figure;
            plot(v_day ,[error_vs_t(:,2:4)  error_vs_t_SelfSim(:,2:4)]');
            
            legend('Bias EOF','RMSE EOF','Min. dist. EOF',...
                'Bias SelfSim','RMSE SelfSim','Min. dist. SelfSim');
        else
            plot(v_day ,error_vs_t(:,2:4)');
            legend('Bias','RMSE','Min. dist.');
        end
        drawnow;
        eval( ['print -depsc ' model.folder.folder_simu ...
            '/error_along_time.eps']);
        
        
        color_95 = [0.9 0.9 1];
        color_50 = [0.6 0.6 1];
        color_ref = [0.9 0 0];
        LineWidth_ref = 2;
        LineStyle_ref = '-.';
        taille_police = 12;
        
        iii_day_keep = ( v_day <= last_day_plot_qq );
        
        if model.sigma.sto
            
            if strcmp(model.sigma.type_spectrum , 'EOF')
                plot_error_ensemble_comp_EOF_SelfSim
                %% Qantiles
                qq_plot_SelfSim(:,:,day_num-first_day+1,:) = subsampl_qq_SelfSim;
                qq_plot_EOF(:,:,day_num-first_day+1,:) = subsampl_qq_EOF;
                qq_ref(:,:,day_num-first_day+1) = subsampl_T_adv_part_HR;
                
                if day_num-first_day+1 > 1
                    
                    %% Quantiles
                    %                 qq_plot(:,:,t_loop-t_ini+1) = subsampl_qq;
                    %                 qq_ref(:,:,t_loop-t_ini+1) = subsampl_T_adv_part_HR;
                    for p = 1:4
                        
                        width = 5;
                        height = 3.2;
                        figure30=figure(30);
                        close(figure30)
                        figure30=figure(30);
                        set(figure30,'Units','inches', ...
                            'Position',[0 0 width height], ...
                            'PaperPositionMode','auto');
                        
                        subplot(2,1,1)
                        hold on;
                        h_0_95 = area( v_day(iii_day_keep), [ squeeze( qq_plot_EOF(p,1,iii_day_keep) ), ...
                            squeeze( qq_plot_EOF(p,4,iii_day_keep) - qq_plot_EOF(p,1,iii_day_keep) ) ]);
                        set (h_0_95(1), 'FaceColor', 'none');
                        set (h_0_95(2), 'FaceColor', color_95);
                        set (h_0_95, 'LineStyle', '-', 'LineWidth', 1, 'EdgeColor', 'none');
                        
                        h_0_50 = area( v_day(iii_day_keep), [ squeeze( qq_plot_EOF(p,2,iii_day_keep) ), ...
                            squeeze( qq_plot_EOF(p,3,iii_day_keep) - qq_plot_EOF(p,2,iii_day_keep) ) ]);
                        set (h_0_50(1), 'FaceColor', 'none');
                        set (h_0_50(2), 'FaceColor', color_50);
                        set (h_0_50, 'LineStyle', '-', 'LineWidth', 1, 'EdgeColor', 'none');
                        
                        % Raise current axis to the top layer, to prevent it
                        % from being hidden by the grayed area
                        set (gca, 'Layer', 'top');
                        
                        plot(v_day(iii_day_keep) ,squeeze(qq_ref(p,:,iii_day_keep)), 'Color', color_ref,...
                            'LineStyle', LineStyle_ref, 'LineWidth', LineWidth_ref);
                        hold off
                        set(gca,...
                            'Units','normalized',...
                            'FontUnits','points',...
                            'FontWeight','normal',...
                            'FontSize',taille_police,...
                            'FontName','Times')
                        ylabel('b(m.s$^{-2}$)',...
                            'FontUnits','points',...
                            'interpreter','latex',...
                            'FontSize',taille_police,...
                            'FontName','Times')
                        xlabel('time(day)',...
                            'interpreter','latex',...
                            'FontUnits','points',...
                            'FontWeight','normal',...
                            'FontSize',taille_police,...
                            'FontName','Times')
                        ax_IC = axis; ax_IC(1:2)=[first_day last_day_plot_qq];axis(ax_IC);
                        
                        
                        subplot(2,1,2)
                        hold on;
                        h_0_95 = area( v_day(iii_day_keep), [ squeeze( qq_plot_SelfSim(p,1,iii_day_keep) ), ...
                            squeeze( qq_plot_SelfSim(p,4,iii_day_keep) - qq_plot_SelfSim(p,1,iii_day_keep) ) ]);
                        set (h_0_95(1), 'FaceColor', 'none');
                        set (h_0_95(2), 'FaceColor', color_95);
                        set (h_0_95, 'LineStyle', '-', 'LineWidth', 1, 'EdgeColor', 'none');
                        
                        h_0_50 = area( v_day(iii_day_keep), [ squeeze( qq_plot_SelfSim(p,2,iii_day_keep) ), ...
                            squeeze( qq_plot_SelfSim(p,3,iii_day_keep) - qq_plot_SelfSim(p,2,iii_day_keep) ) ]);
                        set (h_0_50(1), 'FaceColor', 'none');
                        set (h_0_50(2), 'FaceColor', color_50);
                        set (h_0_50, 'LineStyle', '-', 'LineWidth', 1, 'EdgeColor', 'none');
                        
                        % Raise current axis to the top layer, to prevent it
                        % from being hidden by the grayed area
                        set (gca, 'Layer', 'top');
                        
                        
                        plot(v_day(iii_day_keep) ,squeeze(qq_ref(p,:,iii_day_keep)), 'Color', color_ref,...
                            'LineStyle', LineStyle_ref, 'LineWidth', LineWidth_ref);hold off
                        
                        set(gca,...
                            'Units','normalized',...
                            'FontUnits','points',...
                            'FontWeight','normal',...
                            'FontSize',taille_police,...
                            'FontName','Times')
                        ylabel('b(m.s$^{-2}$)',...
                            'FontUnits','points',...
                            'interpreter','latex',...
                            'FontSize',taille_police,...
                            'FontName','Times')
                        xlabel('time(day)',...
                            'interpreter','latex',...
                            'FontUnits','points',...
                            'FontWeight','normal',...
                            'FontSize',taille_police,...
                            'FontName','Times')
%                         title(['Self similar method'],...
%                             'FontUnits','points',...
%                             'FontWeight','normal',...
%                             'interpreter','latex',...
%                             'FontSize',12,...
%                             'FontName','Times')
%                         title(['Ref. VS 95%- and 50%-CI' ...
%                             ' for the self similar method'],...
%                             'FontUnits','points',...
%                             'FontWeight','normal',...
%                             'interpreter','latex',...
%                             'FontSize',12,...
%                             'FontName','Times')
                        ax_IC = axis; ax_IC(1:2)=[first_day last_day_plot_qq];axis(ax_IC);
                        
                        drawnow;
                        eval( ['print -depsc ' model.folder.folder_simu ...
                            '/quantile_along_time_at_pt_' num2str(p) '.eps']);
                    end
                end
            else
                plot_error_ensemble
                %% Quantiles
                qq_plot(:,:,day_num-first_day+1) = subsampl_qq;
                qq_ref(:,:,day_num-first_day+1) = subsampl_T_adv_part_HR;
                
                
                if day_num-first_day+1 > 1
                    for p = 1:4
                        
                        width = 5;
                        height = 3.2;
                        figure30=figure(30);
                        close(figure30)
                        figure30=figure(30);
                        set(figure30,'Units','inches', ...
                            'Position',[0 0 width height], ...
                            'PaperPositionMode','auto');
                        
                        hold on;
                        h_0_95 = area( v_day(iii_day_keep), [ squeeze( qq_plot(p,1,iii_day_keep) ), ...
                            squeeze( qq_plot(p,4,iii_day_keep) - qq_plot(p,1,iii_day_keep) ) ]);
                        set (h_0_95(1), 'FaceColor', 'none');
                        set (h_0_95(2), 'FaceColor', color_95);
                        plot(v_day(iii_day_keep) ,squeeze(qq_ref(p,:,iii_day_keep)), 'Color', color_ref,...
                            'LineStyle', LineStyle_ref, 'LineWidth', LineWidth_ref);
                        
                        h_0_50 = area( v_day, [ squeeze( qq_plot(p,2,iii_day_keep) ), ...
                            squeeze( qq_plot(p,3,iii_day_keep) - qq_plot(p,2,iii_day_keep) ) ]);
                        set (h_0_50(1), 'FaceColor', 'none');
                        set (h_0_50(2), 'FaceColor',color_50 );
                        set (h_0_50, 'LineStyle', '-', 'LineWidth', 1, 'EdgeColor', 'none');
                        
                        % Raise current axis to the top layer, to prevent it
                        % from being hidden by the grayed area
                        set (gca, 'Layer', 'top');
                        
                        plot(v_day(iii_day_keep) ,squeeze(qq_ref(p,:,iii_day_keep)), 'Color', color_ref,...
                            'LineStyle', LineStyle_ref, 'LineWidth', LineWidth_ref);
                        hold off
                        set(gca,...
                            'Units','normalized',...
                            'FontUnits','points',...
                            'FontWeight','normal',...
                            'FontSize',taille_police,...
                            'FontName','Times')
                        ylabel('b(m.s$^{-2}$)',...
                            'FontUnits','points',...
                            'interpreter','latex',...
                            'FontSize',taille_police,...
                            'FontName','Times')
                        xlabel('time(day)',...
                            'interpreter','latex',...
                            'FontUnits','points',...
                            'FontWeight','normal',...
                            'FontSize',taille_police,...
                            'FontName','Times')
%                         title(['Ref. VS 95%- and 50%-CI'],...
%                             'FontUnits','points',...
%                             'FontWeight','normal',...
%                             'interpreter','latex',...
%                             'FontSize',12,...
%                             'FontName','Times')
                        ax_IC = axis; ax_IC(1:2)=[first_day last_day_plot_qq];axis(ax_IC);
                        
                        drawnow;
                        eval( ['print -depsc ' model.folder.folder_simu ...
                            '/quantile_along_time_at_pt_' num2str(p) '.eps']);
                    end
                end
            end
        end
        
        %%
        fft_w = SQG_large_UQ(model, fft_b(:,:,:,1));
        fct_sigma_spectrum_abs_diff_postprocess(...
            model,fft_w,true,day);
        if model.sigma.sto
            switch model.sigma.type_spectrum
                case 'SelfSim_from_LS'
                    if model.sigma.time_smooth.bool
                        subplot(1,2,2);ax = axis;
                        plot_abs_diff_from_sigma_postprocess_add(model,...
                            fft2(sigma_dBt_on_sq_dt), [0.0 0.7 0.0]);
                        %  fft2(sigma_dBt_on_sq_dt), [0.0 0.5 0.0]);
                        subplot(1,2,2);axis(ax);
                    end
                case {'EOF' ,'Euler_EOF'}
                        type_plot = [0.0 0.5 0.0];
                    subplot(1,2,2);ax = axis;

                    plot_abs_diff_from_sigma_postprocess_add(model,...
                        fft2(sigma_dBt_on_sq_dt), type_plot );
                    subplot(1,2,2);axis(ax);
                case  'Band_Pass_w_Slope'
                    fct_sigma_spectrum(model,fft_w,nan,day);
            end
            
            drawnow
            eval( ['print -depsc ' model.folder.folder_simu ...
                '/AbsDiffByScale_sigma_dB_t_PostProcess/'...
                day '.eps']);
            
        end
        
        %% Plot variance tensor
        fct_plot_VarTensor(model,day);
    end
    
end
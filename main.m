%%%%%%%%%%%%%%%%%%%%
%%% Main
%%%%%%%%%%%%%%%%%%%%
init;

% Enable plots or not
plots_bool = true;

%% Main parameters to choose

%% Flow simulation parameters

% Type of dynamics
dynamics = 'SQG';
%dynamics = '2D'; % Not tested. To use at your own risk

% Duration of the simulation (in seconds)
advection_duration = 3600*24*20; % 20 days

% Resolution
resolution = 128;
% The number of grid point is resolution^2
% It has to be an even integer

% Reference resolution
resolution_HR = 8 * resolution % for post processing plots

%% Initial condition

% Type of initial condtion
type_data ='Constantin_case2' % 
% %     case 2 of Constantin et al., 1994,
% %         https://doi.org/10.1088/0951-7715/7/6/001
% %     used in Resseguier et al., 2020a,
% %         https://doi.org/10.5194/npg-27-209-2020
% type_data = 'disym_Vortices' %
% %     2 large dysymmetric  anticyclones and cyclones
% %     used in Resseguier et al., 2020a, 
% %         https://doi.org/10.5194/npg-27-209-2020
% type_data = 'Vortices' 
% %     2 large anticyclones and 2 large cyclones
% %     used in Resseguier et al., 2017b, 
% %         https://doi.org/10.1080/03091929.2017.1312101 
% %     and Resseguier et al., 2020b, 
% %         https://doi.org/10.1007/s11831-020-09437-x 
% type_data ='Vortices2' 
% %     same as 'Vortices' but properly periodized (by Pierre Derian).
% type_data ='Perturbed_vortices' 
% %     Same flow with slight small-scale modifications
% %     used in Resseguier et al., 2020b, 
% %         https://doi.org/10.1007/s11831-020-09437-x 
% type_data ='Spectrum' 
% %     Gaussian random field with a spectrum slope deined by
% %     the variable slop_b_ini (default value  = -5/3)
% %     used in Resseguier et al., 2020b, 
% %         https://doi.org/10.1007/s11831-020-09437-x 
% type_data ='Constantin_case1' 
% %     case 1 of Constantin et al., 1994, 
% %         https://doi.org/10.1088/0951-7715/7/6/001
% type_data ='Zero' 
% %     Field equal to zero everywhere

% Spectrum slope of the initial condition
% (for type_data = 'Spectrum' only )
switch dynamics
    case 'SQG'
        slope_b_ini = - 5/3;
    case '2D'
        slope_b_ini = - 3;
    otherwise
        error('Unknown type of dynamics');
end

% Begin simulation from a precomputed field?
use_save = false
% In this case, which day should be used as initialisation
day_save = 100;

%% Forcing of the flow

% Forcing or not
forcing = false;

% Type de forcing
% forcing_type = 'Kolmogorov';
% % Forcing of the following type
% % F = ampli_forcing * odg_b * 1/T_caract ...
% %             * cos( 2 freq_f(1) pi x/L_x + 2 freq_f(2) pi y/L_y )
% % F = ampli_forcing * odg_b * 1/T_caract * sin( 2 freq_f pi y/L_y )
forcing_type = 'Spring';
% % Forcing of the following type
% % F = -(1/tau)*b + ampli_forcing * odg_b * 1/T_caract ...
% %         * sin( 2 freq_f(1) pi x/L_x )* sin( 2 freq_f pi y/L_y )

% Amplitude of the forcing
ampli_forcing = 10;

% Spatal frequency of the forcing
freq_f = [3 2];

%% Stochastic parameterization of v'= sigma dBt/dt

% Deterministic or random model
stochastic_simulation = false;
sigma.sto = stochastic_simulation;
% Usual SQG model (stochastic_simulation=false)
% or SQG_MU model (stochastic_simulation=true)

if sigma.sto
    % Type of spectrum for v' = sigma dBt/dt
    % sigma.type_spectrum = 'Band_Pass_w_Slope'; 
    % %     spectrum ~ (cst * k^r) for km < k < kM
    % %     proposed and used in Resseguier et al., 2017b, 
    % %         https://doi.org/10.1080/03091929.2017.1312101 
    % %     used in Resseguier et al., 2020b, 
    % %         https://doi.org/10.1007/s11831-020-09437-x 
    sigma.type_spectrum = 'SelfSim_from_LS'
    % %     ADSD method
    %       Sigma computed from self similarities from the large scales
    % %     proposed (section 2.1.1) and used in Resseguier et al., 2020a, 
    % %         https://doi.org/10.1080/03091929.2017.1312101 
    % %     used in Resseguier et al., 2020b, 
    % %         https://doi.org/10.1007/s11831-020-09437-x 
    % sigma.type_spectrum = 'EOF';
    % %     Empirical Orthogonal Functions (EOF)
    % %     ( Need a pre-computed EOF file )
    % %     used (section 2.2) in Resseguier et al., 2020a, 
    % %         https://doi.org/10.5194/npg-27-209-2020
    % %     proposed and used by Cotter et al., 2019, 2020a,b,c
    % %         https://doi.org/10.1137/18M1167929
    % %         arXiv:1802.05711 ; arXiv:1907.11884
    % %         https://doi.org/10.1007/s10955-020-02524-0
    % sigma.type_spectrum = 'Low_Pass_w_Slope';
    % %     Spectrum cst for k<km 
    % %     and  spectrum = (cst * k^r) for k>km
    % sigma.type_spectrum = 'Low_Pass_streamFct_w_Slope';
    % %     Matern covariance for the streamfunction
    % %     spectrum = cst. * k^2 .* ( 1 + (k/km)^2 )^s )
    % %     So, spectrum ~ cst. * k^2 for k<<km 
    % %     and spectrum ~ cst. * k^r for k>>km 
    % sigma.type_spectrum = 'BB';
    % %     White noise : spectrum = cst
    
    % v' heterogenity 
    % Heterogeneity based on the energy flux epsilon espilon_F
    sigma.hetero_energy_flux = false; 
    % Or Heterogeneity based on the Smagorinsky energy flux 
    %   (i.e. the dissipation espilon_D)
    sigma.hetero_modulation_Smag = false;
    % %     both heterogenity methods are proposed (section 2.1.2)
    % %     and used in Resseguier et al., 2020a, 
    % %         https://doi.org/10.1080/03091929.2017.1312101 
end

% Other sigma parameters (do not change)
Default_sigma_param

% Number of realizations in the ensemble
N_ech = 1;
% %     N_ech is automatically set to 1 in deterministic simulations.
if ~sigma.sto
    N_ech = 1;
end

%% Deterministic subgrid tensor

% (Laplacian) viscosity
Lap_visco.bool = false;

% Hyper-viscosity
HV.bool = true;

% Smagorinsky-like diffusivity/viscosity or Hyper-viscosity
Smag.bool = false;

% Other deterministic subgrid tensor parameters (do not change)
Default_deter_subgrid_tensor

%% Plot parameters

% Choose to plot the dissipation by scale
plot_epsilon_k = false;
if sigma.sto && sigma.hetero_energy_flux
    plot_epsilon_k = true;
end

% Choose to plot dissipations terms 
plot_dissip = false;

%% Additional parameters

% Additional parameters for v'= sigma dBt/dt (do not change)
Additional_sigma_param

% Physical parameters
model = fct_physical_param(dynamics);

%% Gather parameters in a structure called "model"
model.sigma = sigma;
if strcmp(type_data,'Spectrum')
    model.slope_b_ini = slope_b_ini;
end
model.dynamics=dynamics;
model.type_data=type_data;
model.advection.N_ech=N_ech;
model.advection.advection_duration=advection_duration;
model.advection.plot_epsilon_k = plot_epsilon_k;
model.advection.plot_dissip = plot_dissip;
model.advection.plot_moments = plot_moments;
model.advection.forcing.bool = forcing;
model.advection.forcing.ampli_forcing = ampli_forcing;
model.advection.forcing.freq_f = freq_f;
model.advection.forcing.forcing_type = forcing_type;
model.advection.HV = HV;
model.advection.cov_and_abs_diff = false;
model.advection.Lap_visco = Lap_visco;
model.advection.Smag = Smag;
model.advection.use_save = use_save;
model.advection.day_save = day_save;
model.grid.dealias_method = dealias_method; %de-aliasing method
model.plots = plots_bool;
if model.sigma.sto
    if ( strcmp(model.sigma.type_spectrum,'EOF') || ...
            strcmp(model.sigma.type_spectrum,'Euler_EOF') ) ...
            && ( model.sigma.nb_EOF > model.advection.N_ech )
        warning(['The number of EOF is larger than the ensemble size.' ...
            ' Some EOFs are hence removed.']);
        model.sigma.nb_EOF = model.advection.N_ech;
    end
end

%% Random generator
rng('default'); %  The default settings are the Mersenne Twister with seed 0.

%% Generating initial buoyancy
[fft_buoy,model] = fct_buoyancy_init(model,resolution);

%% Advection
[fft_buoy_final, model] = fct_fft_advection_sto_mat(model, fft_buoy);

%% Post-process plots
warning('You need to pre-compute a deterministic reference for this part of the code');
if ~ use_save
    first_day = 0;
    day_save = 0;
else
    first_day = 100;
    day_save = 101;
end

if plots_bool
    post_process_error_grid(model.sigma.sto,...
        model.type_data,resolution,resolution_HR,...
        model.advection.forcing.bool,model.sigma,...
        model.advection.Lap_visco,model.advection.HV,...
        model.advection.Smag,model.advection.N_ech,...
        first_day,advection_duration,day_save);
end

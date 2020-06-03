function [fft_b, model] = fct_fft_advection_sto_mat(model,  fft_b)
% Advection of buoyancy using SQG or SQG_MU model
%

tic
%% Folder to save plots and files
folder_to_save

%% Grid

% Spatial grid
model.grid.x = model.grid.dX(1)*(0:model.grid.MX(1)-1);
model.grid.y = model.grid.dX(2)*(0:model.grid.MX(2)-1);

% Grid in Fourier space
model = init_grid_k (model);

% Ensemble size
N_ech=model.advection.N_ech;

%% Initialisation of the spatial fields

% Remove aliasing
fft_b(model.grid.k.ZM(1),:,:,:)=0;
fft_b(:,model.grid.k.ZM(2),:,:)=0;

% Initial large-scale velocity
fft_w = SQG_large_UQ(model, fft_b);
w=real(ifft2(fft_w));

% Create several identical realizations of the intial buoyancy
fft_b = repmat(fft_b(:,:,1),[1 1 1 model.advection.N_ech]);

%% Choice of the variance tensor a
init_variance_tensor

%% Initialize forcing
init_forcing

%% Initialize hyperviscosity
init_hyperviscosity

%% Choice of time step : CFL
model.advection.dt_adv = fct_CFL(model,w);

%% Fourier transform of the kernel \tilde sigma
fct_init_sigma

%% Print informations

if model.plots
    % Printing some information
    fprintf(['The initial condition is ' model.type_data ' \n'])
    fprintf(['1/k_c is equal to ' num2str(1/model.sigma.k_c) ' m \n'])
    fprintf(['Time step : ' num2str(model.advection.dt_adv) ' seconds \n']);
    fprintf(['Time of advection : ' num2str(...
        model.advection.advection_duration/3600/24) ' days \n']);
    fprintf(['Ensemble size : ' num2str(N_ech) ' realizations \n']);
    fprintf(['Resolution : ' num2str(model.grid.MX(1)) ' x ' ...
        num2str(model.grid.MX(2)) ' \n']);
    if model.advection.HV.bool
        str_subgridtensor = 'Hyper-viscosity';
    elseif model.advection.Lap_visco.bool
        str_subgridtensor = 'Laplacian diffusion';
    else
        str_subgridtensor = 'None';
    end
    if ( model.advection.HV.bool | model.advection.Lap_visco.bool) & ...
            model.advection.Smag.bool
        str_subgridtensor = ['Heterogneous ' str_subgridtensor];
    end
    fprintf(['Deterministic subgrid tensor : ' str_subgridtensor ' \n']);
    fprintf(['Model type : ' add_subgrid_deter ' \n']);
    if model.sigma.sto | model.advection.Smag.bool
        fprintf(['Details : ' subgrid_details ' \n']);
    end
    if model.sigma.sto
        fprintf(['type spectrum sigma :' model.sigma.type_spectrum ' \n']);
    end
end

%% Use a saved files of a former simulation ?
if model.advection.use_save
    warning(['The run begin from an older file instead of from the' ...
        'initial condition']);
    day = num2str(model.advection.day_save);
    name_file = [model.folder.folder_simu '/files/' day '.mat']; clear day
    clear fft_b
    model_ref = model;
    
    load(name_file)
    model = model_ref;
    
    if (exist('fft_b','var')==1)
    elseif (exist('fft_buoy_part','var')==1)
        fft_b = fft_buoy_part; clear fft_buoy_part
    elseif (exist('fft_T_adv_part','var')==1)
        fft_b = fft_T_adv_part; clear fft_T_adv_part
    elseif (exist('fft_T_adv_part','var')==1)
        fft_b = fft_T_adv_part; clear fft_T_adv_part
    elseif (exist('fft_tracer_part','var')==1)
        fft_b = fft_tracer_part; clear fft_tracer_part
    elseif (exist('fft_buoy_part_ref','var')==1)
        fft_b = fft_buoy_part_ref; clear fft_buoy_part_ref
    else
        error('Cannot find buoyancy field')
    end
    if model.sigma.sto
        if size(fft_b,4) < model.advection.N_ech
            if size(fft_b,4) == 1
                fft_b = repmat ( fft_b, [ 1 1 1 model.advection.N_ech ]);
                clear w fft_w
            else
                error('The number of realisation of the saved file is too low');
            end
        end
        if size(fft_b,4) > model.advection.N_ech
            warning(['The number of realisation of the saved file is too high.' ...
                ' Some realisations are hence removed.']);
            fft_b(:,:,:,(model.advection.N_ech+1):end) = [];
        end
        if ( strcmp(model.sigma.type_spectrum,'EOF') || ...
                strcmp(model.sigma.type_spectrum,'Euler_EOF') ) ...
                && ( model.sigma.nb_EOF > model.advection.N_ech )
            warning(['The number of EOF is larger than the ensemble size.' ...
                ' Some EOFs are hence removed.']);
            model.sigma.nb_EOF = model.advection.N_ech;
            sigma(:,:,:,:,(model.advection.N_ech+1):end) = [];
        end
    end
    
    % Version of matlab
    vers = version;
    year = str2double(vers(end-5:end-2));
    subvers = vers(end-1);
    model.folder.colormap_freeze = ...
        (  year < 2014 || ( year == 2014 && strcmp(subvers,'a') ) );
    if ~isfield(model,'plots')
        model.plots = true;
    end
    if ~(exist('time','var')==1)
        time =t*model.advection.dt_adv;
    end
    day_last_plot = floor(time/24/3600);
else
    time = 0;
    w_fv = w;
    day_last_plot = -inf;
end

%% Loop over time
tt_last = -inf;
while time < model.advection.advection_duration
    %% Time-correlated velocity
    fft_w = SQG_large_UQ(model, fft_b);
    w = real(ifft2(fft_w));
    
    % Optional ADSD plot comparaison of Self Sim and EOF methods
    Plot_ADSD_SelfSim_EOF
        
    %% Time-uncorrelated velocity (isotropic and homogeneous in space)
    % Statistics of v'
    [model,sigma,coef_modulation,sigma_s] = fct_sigma(model,fft_b,fft_w,sigma_s);
    if ~ model.sigma.sto % Deterministic case
        sigma_dBt_on_sq_dt = 0;
    else % Stochastic case : sampling of v'        
        if ~ ( strcmp(model.sigma.type_spectrum,'EOF') || ...             
                strcmp(model.sigma.type_spectrum,'Euler_EOF'))
            % Fourier transform of white noise
            dBt_C_on_sq_dt = fft2( randn( [ model.grid.MX 1 N_ech]));
            % Multiplication by the Fourier transform of the kernel \tilde \sigma
            fft_sigma_dBt_on_sq_dt = bsxfun(@times,sigma,dBt_C_on_sq_dt);
            clear dBt_C_on_sq_dt
            % Homogeneous velocity field
            sigma_dBt_on_sq_dt = real(ifft2(fft_sigma_dBt_on_sq_dt));
            clear fft_sigma_dBt_on_sq_dt
        else
            sigma_dBt_on_sq_dt = sum( sigma .* ...
                randn( [ 1 1 1 N_ech model.sigma.nb_EOF ]) , 5);
        end
        
        % Heterogeneous small-scale velocity
        sigma_dBt_on_sq_dt = bsxfun(@times, sqrt(coef_modulation) , ...
            sigma_dBt_on_sq_dt);
        
        if model.sigma.proj_free_div & any(coef_modulation(:) ~= 1)
            % nrj_before_proj_div = mean(sigma_dBt_on_sq_dt(:).^2)
            sigma_dBt_on_sq_dt = fct_proj_free_div(model,sigma_dBt_on_sq_dt);
            % nrj_after_proj_div = mean(sigma_dBt_on_sq_dt(:).^2)
        end
    end
    
    %% Specify determinstic heterogeneous subgrid model
    if model.advection.Smag.bool
        if model.advection.Lap_visco.bool
            % Heterogeneous dissipation coefficient
            model.advection.coef_diff = fct_coef_diff(model,fft_b);
            % Coefficient coef_Smag to target a specific diffusive scale
            model.advection.coef_diff = ...
                model.advection.Smag.coef_Smag * ...
                model.advection.coef_diff ;
        elseif model.advection.HV.bool
            % Heterogeneous HV coefficient
            model.advection.coef_diff = ...
                fct_coef_HV(model,fft_b);
        end
        % Maximum dissipation coefficient
        model.advection.HV.maxVal = max(model.advection.coef_diff(:));
    end
    
    %% Dynamical time step
    model.advection.dt_adv = fct_CFL(model,w);
    time = time + model.advection.dt_adv;
    
    %% Adding time-correlated and time decorrelated velocity    
    if  model.sigma.sto & ~model.sigma.no_noise
        w = w + sigma_dBt_on_sq_dt/sqrt(model.advection.dt_adv);
    end
    
    %% Transport of tracer
    if ~ model.sigma.sto
        % Runge-Kutta 4 scheme
        fft_b = RK4_fft_advection(model,fft_b, w);
    else
        % Euler scheme
        fft_b = fft_b ...
            + deriv_fft_advection( ...
            model, fft_b, w) ...
            * model.advection.dt_adv;
        clear model_temp
    end
    
    %% Discard particles which have blown up
    iii = isnan(fft_b) | isinf(abs(fft_b));
    if any(iii(:))
        iii=any(any(any(iii,3),2),1);
        if all(iii(:))
            if model.plots
                error('The simulation has blown up');
            else
                warning('One simulation has blown up');
                fprintf('Folder of the simulation : \n');
                fprintf([ model.folder.folder_simu ' \n']);
                return;
            end
        end
        nb_dead_pcl = sum(iii);
        warning([ num2str(nb_dead_pcl) ' particle(s) on ' num2str(N_ech) ...
            ' have(s) blown up and' ...
            ' are(is) resampled uniformly on' ...
            ' the set of the others particles']);
        N_ech_temp = N_ech - nb_dead_pcl;
        fft_b(:,:,:,iii)=[];
        iii_sample = randi(N_ech_temp,nb_dead_pcl,1);
        for k=1:nb_dead_pcl
            fft_b(:,:,:,end+1) = fft_b(:,:,:,iii_sample(k));
        end
    end
    clear iii
    
    %% Plots and save
    Plot_and_save
    
end
toc
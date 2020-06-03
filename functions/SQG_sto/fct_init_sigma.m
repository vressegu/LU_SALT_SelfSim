%% Initialisation of v' statistics and its kernel \tilde sigma

if model.sigma.sto
    if strcmp(model.sigma.type_spectrum,'SelfSim_from_LS')
        [sigma, ~, tr_a ,....
            model.sigma.slope_sigma,...
            model.sigma.offset_spectrum_a_sigma, ...
            model.sigma.km_LS ]...
            = fct_sigma_spectrum_abs_diff_mat(model,...
            repmat(fft_w,[1 1 1 model.advection.N_ech]),...
            true,num2str(0));
%         sigma = repmat(sigma,[1 1 1 N_ech]);
        model.sigma.km_LS = repmat(model.sigma.km_LS,[1 N_ech]);
        a0 = tr_a/2;
        model.sigma.a0 = a0;
        model.sigma.a0_on_dt = model.sigma.a0 / model.advection.dt_adv;
        missed_var_small_scale_spectrum = 2*a0;
        % Diffusion coefficient
        model.advection.coef_diff = 1/2 * model.sigma.a0;
    elseif ~ ( strcmp(model.sigma.type_spectrum,'EOF') || ...
            strcmp(model.sigma.type_spectrum,'Euler_EOF') )
        % Fourier transform of the kernel \tilde sigma up to a multiplicative
        % constant
        [sigma, ~, missed_var_small_scale_spectrum ] ...
            = fct_sigma_spectrum(model,fft_w);
        % sigma_on_sq_dt will be used to simulate sigma d B_t
        % missed_var_small_scale_spectrum will be used to set the mulitplicative
        % constant
    end
    
    if model.sigma.Smag.bool | model.sigma.assoc_diff
        % Variance tensor of the simulated small-scle velocity
        model.sigma.a0_LS = 1/2 * missed_var_small_scale_spectrum;
        % Variance tensor of the unsimulated small-scle velocity
        model.sigma.a0_SS = 1/2 * fct_norm_tr_a_theo_Band_Pass_w_Slope(...
            model, ...
            model.sigma.kappaMinUnresolved_on_kappaShanon ...
            *(pi/sqrt(prod(model.grid.dX))), ...
            model.sigma.kappaMaxUnresolved_on_kappaShanon ...
            *(pi/sqrt(prod(model.grid.dX))));
        model.sigma.a0 = model.sigma.a0_LS + model.sigma.a0_SS;
        model.sigma.a0_on_dt = model.sigma.a0 / model.advection.dt_adv;
        
        if model.sigma.assoc_diff
            if ~ strcmp(model.sigma.type_spectrum,'SelfSim_from_LS')
                % model.sigma.a0_LS
                if model.sigma.a0_LS < eps
                    error(['This should not happened and' ...
                        'may casue divison by 0']);
                end
                % Variance tensor
                a0_LS = 2 * model.physical_constant.f0 / model.sigma.k_c^2;
                coef_temp = a0_LS / model.sigma.a0_LS; clear a0_LS;
                
                model.sigma.a0_LS = coef_temp * model.sigma.a0_LS;
                model.sigma.a0_SS = coef_temp * model.sigma.a0_SS;
                model.sigma.a0 = coef_temp * model.sigma.a0;
                model.sigma.a0_on_dt = coef_temp * model.sigma.a0_on_dt;
                
                sigma = sqrt(coef_temp) * sigma;
            end
            % Diffusion coefficient
            model.advection.coef_diff = 1/2 * model.sigma.a0;
            
        elseif model.sigma.Smag.bool
            if model.sigma.a0_SS > eps
                if model.sigma.Smag.epsi_without_noise
                    sigma = ...
                        sqrt(2/(model.sigma.a0_LS+model.sigma.a0_SS)) ...
                        * sigma;
                    model.advection.coef_diff = 1;
                else
                    sigma = sqrt(2/model.sigma.a0_SS) ...
                        * sigma;
                    model.advection.coef_diff = 1 + ...
                        model.sigma.a0_LS / model.sigma.a0_SS ;
                end
            elseif strcmp(model.sigma.type_spectrum,'SelfSim_from_LS')
                % The absolute diffusivity diagnosed from the large-scale
                % kinematic spectrum is too weak. It suggests that there
                % are few small scales and no subgrid terms is needed.
                % Moreover, setting subgris terms to zero prevent numerical
                % errors.
                sigma = zeros(size(sigma));
                % sigma_on_sq_dt = zeros(size(sigma_on_sq_dt));
                model.advection.coef_diff = 0;
            else
                error('Unknow case');
            end
            
        else
            error('Unknown case');
        end
    elseif ~ ( strcmp(model.sigma.type_spectrum,'EOF') || ...
            strcmp(model.sigma.type_spectrum,'Euler_EOF') )
        % Muliplicative constant of the kernel \tilde sigma
        model.sigma.a0_on_dt = model.sigma.a0 / model.advection.dt_adv;
        sigma = ...
            sqrt(2*model.sigma.a0/missed_var_small_scale_spectrum) ...
            * sigma;
        % the factor d=2 is for the dimension d of the space R^d
        % the variance of sigma_dBt/dt is tr(a)/dt = 2 a0 /dt
        
    end
    clear missed_var_small_scale_spectrum
    % warning('This formula may changed if k inf > la resolution');
else
    sigma = 0;
end

if model.sigma.sto & model.sigma.time_smooth.bool
    if  strcmp(model.sigma.type_spectrum,'SelfSim_from_LS') 
        sigma_s = 0;
    else
        error('no need for temporal smoothing since sigma is constant in time');
    end
else
    sigma_s = nan;
end
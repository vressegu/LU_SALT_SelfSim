function [model,sigma,coef_modulation,sigma_s,sigma_n_s] = ...
    fct_sigma(model,fft_b,fft_w,sigma_s,fft_w_resolved)
%% Statistics of v'
%
%% Outputs :
% - model : structure which gathers all simulation information
% - sigma : fft of the v' kernel ( \breve{\sigma} )
%         size : Mx x My x 2 x N_ech
%         sigma = fft of the v' kernel ( \breve{\sigma} )
% - coef_modulation  : optional TKE modulation (default : 1)
%   (size : 1 or Mx x My)
% - sigma_s , sigma_n_s : optional (default : not used/specified).
%                           They can be used for smooth sigma in time
%
%% Inputs :
% - model : structure which gathers all simulation information
% - fft_w : fft of the large-scale velocity
%         size : Mx x My x 2 x 1
% - fft_w_resolved : fft of the fully-resolved velocity
%     (optional : for plots comparison only)
%         size : Mx x My x 2 x 1
%

% Plots intermedate results in order to see what's going on
bool_plot = false;
% bool_plot = true;

% The optional Fourier transform of the HR velocity (fft_w_resolved)
% is here plot comparisons only
if nargin < 5
    fft_w_resolved = nan;
end

if model.sigma.sto
    % ADSD parametrization
    if strcmp(model.sigma.type_spectrum,'SelfSim_from_LS')
        [ sigma, ~, model.sigma.tr_a ,....
            model.sigma.slope_sigma,...
            model.sigma.offset_spectrum_a_sigma, ...
            model.sigma.km_LS ]...
            = fct_sigma_spectrum_abs_diff_mat(...
            model,fft_w,bool_plot,0,fft_w_resolved);
        
        % Smoothing in time of sigma
        if model.sigma.time_smooth.bool
            sigma_n_s = sigma;
            d_sigma_s = 1/model.sigma.time_smooth.tau * ...
                ( - sigma_s + sigma_n_s ) ;
            sigma_s = sigma_s + d_sigma_s * model.advection.dt_adv;
            sigma = sigma_s;
        end
        
        % Variance tensor
        model.sigma.a0 = model.sigma.tr_a/2;
        model.sigma.a0_on_dt = model.sigma.a0 / model.advection.dt_adv;
        
        % Diffusion coefficient
        model.advection.coef_diff = permute( 1/2 * model.sigma.a0 , ...
            [1 3 4 2]);
        if model.sigma.assoc_diff | model.sigma.Smag.bool
            % warning('deal with slope when there are several realizations')
            model.sigma.a0_SS = ...
                1/2 * fct_norm_tr_a_theo_Band_Pass_w_Slope(...
                model, ...
                model.sigma.kappaMinUnresolved_on_kappaShanon ...
                *(pi/sqrt(prod(model.grid.dX))), ...
                model.sigma.kappaMaxUnresolved_on_kappaShanon ...
                *(pi/sqrt(prod(model.grid.dX))));
            model.sigma.a0_LS = ...
                model.sigma.a0 ;
            model.sigma.a0 = ...
                model.sigma.a0 ...
                + model.sigma.a0_SS;
            model.sigma.a0_on_dt = model.sigma.a0 ...
                / model.advection.dt_adv;
            % Diffusion coefficient
            model.advection.coef_diff = permute(...
                1/2 * model.sigma.a0 ,[1 3 4 2]);
            
            if model.sigma.Smag.bool
                iii_non_degenerate_a0_SS = (model.sigma.a0_SS > eps);
                %if model.sigma.a0_SS > eps
                if any(iii_non_degenerate_a0_SS)
                    if model.sigma.Smag.epsi_without_noise
                        sigma(:,:,:,iii_non_degenerate_a0_SS) = ...
                            bsxfun(@times,...
                            permute( ...
                            sqrt(2./(...
                            model.sigma.a0_LS(iii_non_degenerate_a0_SS) ...
                            + model.sigma.a0_SS(iii_non_degenerate_a0_SS))) ...
                            , [ 1 3 4 2]) , ...
                            sigma(:,:,:,(iii_non_degenerate_a0_SS)) );
                    else
                        sigma(:,:,:,iii_non_degenerate_a0_SS) = ...
                            bsxfun( @times, ...
                            permute( ...
                            sqrt(2./...
                            model.sigma.a0_SS(iii_non_degenerate_a0_SS)) ...
                            , [ 1 3 4 2]) , ...
                            sigma(:,:,:,iii_non_degenerate_a0_SS) );
                    end
                end
                if any(~iii_non_degenerate_a0_SS)
                    if strcmp(model.sigma.type_spectrum,'SelfSim_from_LS')
                        % The absolute diffusivity diagnosed from the large-scale
                        % kinematic spectrum is too weak. It suggests that there
                        % are few small scales and no subgrid terms is needed.
                        % Moreover, setting subgris terms to zero prevent numerical
                        % errors.
                        sigma(:,:,:,~iii_non_degenerate_a0_SS) = ...
                            zeros(size(sigma(:,:,:,~iii_non_degenerate_a0_SS)));
                        model.advection.coef_diff(~iii_non_degenerate_a0_SS) = 0;
                    else
                        error('Unknow case');
                    end
                end
            end
            
        end
    elseif model.sigma.Smag.bool | model.sigma.assoc_diff
        model.sigma.a0 = model.sigma.a0_LS + model.sigma.a0_SS;
    else
        % Variance tensor
        model.sigma.a0 = 2 * model.physical_constant.f0 ...
            / model.sigma.k_c^2;
    end
    
    
    if model.sigma.hetero_energy_flux
        
        coef_modulation = ...
            fct_epsilon_k_onLine(model,fft_b,fft_w);
    elseif model.sigma.hetero_modulation | ...
            model.sigma.hetero_modulation_V2
        if ( isfield(model.sigma,'hetero_energy_flux_prefilter')  ...
                &&    model.sigma.hetero_energy_flux_prefilter ) ...
                || (isfield(model.sigma,'hetero_energy_flux_postfilter')  ...
                &&    model.sigma.hetero_energy_flux_postfilter)
            error('not coded yet')
        end
        coef_modulation = ...
            fct_coef_estim_AbsDiff_heterogeneous(model,fft_w);
    elseif model.sigma.hetero_modulation_Smag
        % Heterogeneous dissipation coefficient
        if isfield(model.sigma,'hetero_energy_flux_prefilter')  ...
                &&    model.sigma.hetero_energy_flux_prefilter
            % Pre-filtering
            fft_b_for_modulation = bsxfun(@times, ...
                model.grid.k_aa_sigma_hetero_energy_flux.mask, ...
                fft_b);
        else
            fft_b_for_modulation = fft_b;
        end
        coef_modulation = fct_coef_diff(model,fft_b_for_modulation);
        clear fft_b_for_modulation
        
        if isfield(model.sigma,'hetero_energy_flux_postfilter')  ...
                &&    model.sigma.hetero_energy_flux_postfilter
            % Post-filtering
            coef_modulation = real(ifft2(bsxfun(@times, ...
                model.grid.k_aa_sigma_hetero_energy_flux_post.mask, ...
                fft2(coef_modulation))));
        end
        m_coef_modulation = mean(mean(coef_modulation,1),2);
        coef_modulation = bsxfun( @times, ...
            1./m_coef_modulation, coef_modulation);
        clear m_coef_modulation
    elseif model.sigma.Smag.bool
        % Heterogeneous dissipation coefficient
        coef_modulation = fct_coef_diff(model,fft_b);
        % Coefficient coef_Smag to target a specific diffusive scale
        coef_modulation = model.sigma.Smag.coef_Smag * coef_modulation ;
        
        % Coefficient coef_Smag to target a specific diffusive scale
        iii_non_degenerate_a0_SS = (model.sigma.a0_SS > eps);
        %if model_sampl(sampl).sigma.a0_SS > eps
        model.advection.coef_diff = zeros([1 1 1 model.advection.N_ech]);
        if model.sigma.Smag.epsi_without_noise
            model.advection.coef_diff(iii_non_degenerate_a0_SS) = 1;
        else
            % Taking into account the noise in the energy budget
            model.advection.coef_diff(iii_non_degenerate_a0_SS) = ...
                permute( ...
                (1 + model.sigma.a0_LS(iii_non_degenerate_a0_SS) ./ ...
                model.sigma.a0_SS(iii_non_degenerate_a0_SS)) , ...
                [1 4 3 2] );
        end
        if any(~iii_non_degenerate_a0_SS)
            if strcmp(model.sigma.type_spectrum,'SelfSim_from_LS')
                % The absolute diffusivity diagnosed from the large-scale
                % kinematic spectrum is too weak. It suggests that there
                % are few small scales and no subgrid terms is needed.
                % Moreover, setting subgris terms to zero prevent numerical
                % errors.
                model.advection.coef_diff(~iii_non_degenerate_a0_SS) = 0;
            else
                error('Unknow case');
            end
        end
        
        if model.sigma.Smag.SS_vel_homo
            coef_modulation = mean(mean(coef_modulation,2),1);
        end
    else
        coef_modulation = 1;
    end
    %% Variance tensor
    model.advection.coef_diff = ...
        bsxfun( @times, ...
        model.advection.coef_diff,...
        coef_modulation) ;
    if model.sigma.assoc_diff
        
        model.sigma.a0 = ...
            model.sigma.a0_LS ...
            + model.sigma.a0_SS;
        model.sigma.a0_on_dt = ...
            model.sigma.a0 / ...
            model.advection.dt_adv;
        % Diffusion coefficient
        model.advection.coef_diff = bsxfun(@times, ...
            permute( 1/2 * model.sigma.a0 , [ 1 3 4 2]) , ...
            coef_modulation ) ;
    end
    
    
    
    
    % Maximum of the variance tensor
%     a0_temp = nan([model.advection.N_ech 1]);
    if size(coef_modulation,4)==1
        coef_modulation = repmat(coef_modulation,[1 1 1 model.advection.N_ech]);
    end
    
    model.sigma.a0 = bsxfun(@times, ...
        permute(model.sigma.a0 , [ 1 3 4 2]) , ...
        coef_modulation ) ;
    
    if model.sigma.Smag.bool
        iii_non_degenerate_a0_SS = (model.sigma.a0_SS > eps);
        %if model_sampl(sampl).sigma.a0_SS > eps
        if model.sigma.Smag.epsi_without_noise
            model.sigma.a0(iii_non_degenerate_a0_SS) = ...
                model.sigma.a0(iii_non_degenerate_a0_SS) ...
                ./ (model.sigma.a0_SS(iii_non_degenerate_a0_SS) ...
                + model.sigma.a0_LS(iii_non_degenerate_a0_SS));
        else
            % Taking into account the noise in the energy budget
            model.sigma.a0(iii_non_degenerate_a0_SS) = ...
                model.sigma.a0(iii_non_degenerate_a0_SS) ...
                ./ model.sigma.a0_SS(iii_non_degenerate_a0_SS);
        end
        if any(~iii_non_degenerate_a0_SS)
            if strcmp(model.sigma.type_spectrum,'SelfSim_from_LS')
                % The absolute diffusivity diagnosed from the large-scale
                % kinematic spectrum is too weak. It suggests that there
                % are few small scales and no subgrid terms is needed.
                % Moreover, setting subgris terms to zero prevent numerical
                % errors.
                model.sigma.a0(~iii_non_degenerate_a0_SS) = 0;
            else
                error('Unknow case');
            end
        end
    end
    
    a0_temp = max(max( model.sigma.a0 ,[],2), [],1);
    model.sigma.a0 = max(a0_temp(:)) ; clear a0_temp;
    
else
    sigma = 0;
    model.sigma.a0 = 0;
    coef_modulation=nan;
    sigma_s=nan;
    sigma_n_s=nan;
end
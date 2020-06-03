%% Choice of the variance tensor a

if model.sigma.sto & ...
        ( model.sigma.Smag.bool | model.sigma.assoc_diff )
    if strcmp(model.sigma.type_spectrum,'SelfSim_from_LS')
        [~, ~, tr_a ,...
            model.sigma.slope_sigma,...
            model.sigma.offset_spectrum_a_sigma, ...
            model.sigma.km_LS ...
            ]= fct_sigma_spectrum_abs_diff(...
            model,fft_w,true,num2str(0));
        %         a0 = tr_a/2;
        %         model.sigma.a0 = max(a0);
        a0_LS_temp = tr_a/2;
    else
        % Variance tensor of the simulated small-scle velocity
        a0_LS_temp = 1/2 * model.sigma.fct_tr_a(model, ...
            model.sigma.kappamin_on_kappamax ...
            * model.sigma.kappamax_on_kappaShanon ...
            * (pi/sqrt(prod(model.grid.dX))), ...
            model.sigma.kappamax_on_kappaShanon ...
            * (pi/sqrt(prod(model.grid.dX))));
        %     a0_LS_temp = 1/2 * fct_norm_tr_a_theo(model, ...
        %         model.sigma.kappamin_on_kappamax ...
        %         * model.sigma.kappamax_on_kappaShanon ...
        %         * (pi/sqrt(prod(model.grid.dX))), ...
        %         model.sigma.kappamax_on_kappaShanon ...
        %         * (pi/sqrt(prod(model.grid.dX))), ...
        %         ( 3 - model.sigma.slope_sigma )/2 );
        % Variance tensor of the unsimulated small-scle velocity
        %     a0_SS_temp = 1/2 * fct_norm_tr_a_theo(model, ...
    end
    a0_SS_temp = 1/2 * fct_norm_tr_a_theo_Band_Pass_w_Slope(model, ...
        model.sigma.kappaMinUnresolved_on_kappaShanon ...
        * (pi/sqrt(prod(model.grid.dX))), ...
        model.sigma.kappaMaxUnresolved_on_kappaShanon ...
        * (pi/sqrt(prod(model.grid.dX))));
    % Multiplicative foctor for the dissipation coefficient
    if model.sigma.assoc_diff
        % Variance tensor
        if strcmp(model.sigma.type_spectrum,'SelfSim_from_LS')
            model.sigma.a0 = a0_LS_temp + a0_SS_temp ;
        else
            model.sigma.a0 = ( 1 + a0_SS_temp / a0_LS_temp )...
                * 2 * model.physical_constant.f0 / model.sigma.k_c^2;
        end
        % Diffusion coefficient
        model.advection.coef_diff = 1/2 * model.sigma.a0;
    elseif model.sigma.Smag.bool
        if model.sigma.Smag.epsi_without_noise
            coef_diff_temp = 1;
        else
            if a0_SS_temp > eps
                coef_diff_temp = 1 + a0_LS_temp / a0_SS_temp ;
            elseif strcmp(model.sigma.type_spectrum,'SelfSim_from_LS')
                coef_diff_temp = 0;
            else
                error('Unknow case');
            end
        end
        % Heterogeneous dissipation coefficient
        coef_diff_aa_temp = coef_diff_temp * ...
            fct_coef_diff(model,fft_b);
        % Coefficient coef_Smag to target a specific diffusive scale
        model.sigma.Smag.coef_Smag = ...
            ( model.sigma.Smag.kappamax_on_kappad ...
            * sqrt(prod(model.grid.dX))/pi ) ^ 2;
        coef_diff_aa_temp = model.sigma.Smag.coef_Smag * ...
            coef_diff_aa_temp ;
        % Maximum of the total variance tensor
        model.sigma.a0 = 2 * max(coef_diff_aa_temp(:));
        clear a0_LS_temp a0_SS_temp coef_diff_temp coef_diff_aa_temp
    else
        error('Unknown case');
    end
    
elseif model.sigma.sto & ...
        strcmp(model.sigma.type_spectrum,'SelfSim_from_LS')
    [~, ~, tr_a ,....
        model.sigma.slope_sigma,...
        model.sigma.offset_spectrum_a_sigma, ...
        model.sigma.km_LS ] ...
        = fct_sigma_spectrum_abs_diff_mat(...
        model,repmat(fft_w,[1 1 1 model.advection.N_ech]),true,0);
%         ] = fct_sigma_spectrum_abs_diff(model,fft_w,true,num2str(0));
%     % sigma_on_sq_dt = (1/sqrt(model.advection.dt_adv)) * sigma; clear sigma
    a0 = tr_a/2;
    model.sigma.a0 = max(a0);
    %     model.sigma.k_c = ...
    %         sqrt( 2 * model.physical_constant.f0 / model.sigma.a0 );
    %     % Diffusion coefficient
    %     model.advection.coef_diff = 1/2 * model.sigma.a0;
elseif model.sigma.sto & ...
        ( strcmp(model.sigma.type_spectrum,'EOF') || ...
            strcmp(model.sigma.type_spectrum,'Euler_EOF'))
    % sigma = nan;
    
    current_folder = pwd;
    cd(model.folder.folder_simu);
    cd ..
    model.folder.folder_EOF = [ pwd '/folder_EOF'];
    cd(current_folder);
    
    % Load precomputed EOFs and correspond variance tensor
    load([ model.folder.folder_EOF '/EOF.mat'],'EOF');    
    
    EOF = permute(EOF,[1 2 3 5 4]);
    if isinf(model.sigma.nb_EOF)
        model.sigma.nb_EOF = size(EOF,5);
    else
        EOF = EOF(:,:,:,:,1:model.sigma.nb_EOF);
    end    
    %model.sigma.nb_EOF = size(EOF,5);
    sigma = EOF; clear EOF;
    
    a_xx = sum( sigma .* permute(sigma,[1 2 4 3 5]) ,5);
    a0 =  1/2 * ( a_xx(:,:,1,1) + a_xx(:,:,2,2) );
    a0 =  mean(a0(:));
    
    % Filtering the variance tensor at large scales
    sigma_loc = permute(a_xx,[3 4 1 2]);
    sigma_loc = multi_chol(sigma_loc);
    sigma_loc = permute(sigma_loc,[3 4 1 2]);
    sigma_loc_aa = fft2(sigma_loc);
    sigma_loc_aa = real(ifft2( bsxfun(@times, sigma_loc_aa,...
        model.grid.k_aa_LS.mask) ));
    sigma_loc_aa = permute(sigma_loc_aa,[3 4 1 2]);
    a_xx_aa = multiprod(sigma_loc_aa,multitrans(sigma_loc_aa));
    a_xx_aa = permute(a_xx_aa,[3 4 1 2]);
        
    % Plots
    figure(12);
    for d1=1:2
        for d2=1:2
            subplot(2,2,1+d1-1+2*(d2-1))
            imagesc(model.grid.x,model.grid.y,a_xx(:,:,d1,d2)');
            axis equal; axis xy;colorbar;
            title(['$a_{' num2str(d1) num2str(d2) '}$'],...
                'FontUnits','points',...
                'FontWeight','normal',...
                'interpreter','latex',...
                'FontSize',12,...
                'FontName','Times')
        end
    end
    drawnow
    eval( ['print -depsc ' model.folder.folder_simu ...
        '/Variance_tensor.eps']);
    figure(13);
    for d1=1:2
        for d2=1:2
            subplot(2,2,1+d1-1+2*(d2-1))
            imagesc(model.grid.x,model.grid.y,a_xx_aa(:,:,d1,d2)');
            axis equal; axis xy;colorbar;
            title(['Smooth $a_{' num2str(d1) num2str(d2) '}$'],...
                'FontUnits','points',...
                'FontWeight','normal',...
                'interpreter','latex',...
                'FontSize',12,...
                'FontName','Times')
        end
    end
    drawnow
    eval( ['print -depsc ' model.folder.folder_simu ...
        '/Smooth_Variance_tensor.eps']);
    
    % Variance tensor
    [a0 max(a_xx_aa(:)) ]
    a0 = max( [a0 max(a_xx_aa(:)) ]);
    model.sigma.a0 = a0; clear a0;
    % Diffusion coefficient
    a_xx_aa = permute( a_xx_aa, [ 1 2 4 5 3]);
    model.advection.coef_diff = 1/2 * a_xx_aa; clear a_xx;
else
    % Variance tensor
    model.sigma.a0 = 2 * model.physical_constant.f0 / model.sigma.k_c^2;
    % Diffusion coefficient
    model.advection.coef_diff = 1/2 * model.sigma.a0;
end


if model.sigma.sto & model.sigma.hetero_energy_flux
    if strcmp(model.sigma.type_spectrum,'SelfSim_from_LS')
        model.sigma.km_LS = repmat(model.sigma.km_LS,[1 1 1 N_ech]);
    end
    coef_modulation = fct_epsilon_k_onLine(model,fft_b);
elseif model.sigma.sto & model.sigma.hetero_modulation_Smag
    % Heterogeneous dissipation coefficient
    coef_modulation = fct_coef_diff(model,fft_b);
    m_coef_modulation = mean(mean(coef_modulation,1),2);
    coef_modulation = bsxfun( @times, ...
        1./m_coef_modulation, coef_modulation);
    clear m_coef_modulation
    % coef_modulation = fct_coef_diff(model,fft_b);
elseif model.sigma.sto & ...
        (model.sigma.hetero_modulation | model.sigma.hetero_modulation_V2)
    coef_modulation = fct_coef_estim_AbsDiff_heterogeneous(model,fft_w);
else
    coef_modulation = 1;
end
model.sigma.a0 = model.sigma.a0 * max(coef_modulation(:));
% warning(['The CFL may changed here ...
% in the case of heterogreneous variance tensor'])
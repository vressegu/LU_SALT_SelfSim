function model = init_grid_k (model, bool_skip, bool_lowCost)
%% function model = init_grid_k (model)
% Create a grid in the Fourier space
%
% From struct "model", expecting:
%   - grid.MX (grid size);
%   - grid.dX (spatial sampling step);
%   - grid.dealias_method (de-aliasing method for pseudo-spectral codes)
%
% Modified by P. DERIAN 2017-06-05:
%   - added various dealiasing methods (model.grid.dealias_method parameter).
% Modified by P. DERIAN 2017-06-13:
%   - added ZM = PX+1 as a member of model.grid.k.

if nargin < 2
    bool_skip = false;
end
if nargin < 3
    bool_lowCost = false;
end

% check grid size is even
if any( mod(model.grid.MX,2)~=0)
    error('the number of grid points by axis need to be even');
end
% if bool_lowCost
%     if any( mod(model.grid.MX,1)~=0)
%         error('the number of grid points by axis need to be integer');
%     end    
% else
%     if any( mod(model.grid.MX,2)~=0)
%         error('the number of grid points by axis need to be even');
%     end
% end
PX = model.grid.MX/2;
ZM = PX + 1; %index of the single high-freq mode to be zero'ed out.

bool_heterogeneous_model = (~bool_skip) && (...
    ( model.advection.HV.bool | model.advection.Lap_visco.bool) & ...
    model.advection.Smag.bool | ...
    (model.sigma.sto && (model.sigma.Smag.bool | ...
    ( ( model.sigma.hetero_modulation ...
    | model.sigma.hetero_modulation_V2 ...
    | model.sigma.hetero_modulation_Smag ...
    | model.sigma.hetero_energy_flux ) ...
    & strcmp(model.sigma.type_spectrum,'SelfSim_from_LS')) )) | ...
    ( model.sigma.sto && ...
    ( strcmp(model.sigma.type_spectrum,'EOF') || ...
    strcmp(model.sigma.type_spectrum,'Euler_EOF') ) ) );
model.advection.bool_heterogeneous_model = bool_heterogeneous_model;


%% "Unstable" Fourier grid
% for homogeneous diffusion Laplacian(b), hyper-viscosity Laplacian^p(b),
% SQG relationship...
nx = [ 0:(PX(1)-1) 0 (1-PX(1)):-1]; %NB: the central single high-freq is zero'ed-out
ny = [ 0:(PX(2)-1) 0 (1-PX(2)):-1];
kx = (1./model.grid.MX(1)) .* nx;
ky = (1./model.grid.MX(2)) .* ny;
kx = (2.*pi/model.grid.dX(1)) .* kx; %as wavenumbers
ky = (2.*pi/model.grid.dX(2)) .* ky;
d_kx = kx(2)-kx(1);
d_ky = ky(2)-ky(1);
dK = [ d_kx d_ky ];
kx_ref = dK(1)*(0:(model.grid.MX(1)-1));
ky_ref = dK(2)*(0:(model.grid.MX(2)-1));
Period_k = dK .* model.grid.MX ;
% kx_ref = kx;
% ky_ref = ky;
if ~bool_lowCost
    % the 2D grid
    [kx,ky] = ndgrid(kx,ky);
    k2 = kx.^2+ky.^2;
    k2(ZM(1),:) = 0.; %de-alias the single high freq
    k2(:,ZM(2)) = 0.;
    k = sqrt(k2); %the modulus
    % Specific operators
    if ~ bool_skip
        switch model.dynamics
            case 'SQG'
                over_k = 1./k;
            case '2D'
                over_k = -1./k.^2;
            otherwise
                error('Unknown type of dynamics');
        end
        over_k( k==0 ) = 0.;
        on_k2 = 1./k2;
        % on_k2 = 1./k.^2;
        on_k2( k==0 ) = 0.;
    end
    % Projection on the spce of free divergent function
    proj_free_div(:,:,1,1,1) = 1 - kx.*kx./k2;
    proj_free_div(:,:,1,1,2) = 0 - kx.*ky./k2;
    proj_free_div(:,:,2,1,1) = 0 - ky.*kx./k2;
    proj_free_div(:,:,2,1,2) = 1 - ky.*ky./k2;
    proj_free_div(1,1,:,1,:) = 0;
    proj_free_div(ZM(1),:,:,1,:) = 0;
    proj_free_div(:,ZM(2),:,1,:) = 0;
    % proj_free_div(:,:,1,1) = 1 - kx.*kx./k2;
    % proj_free_div(:,:,1,2) = 0 - kx.*ky./k2;
    % proj_free_div(:,:,2,1) = 0 - ky.*kx./k2;
    % proj_free_div(:,:,2,2) = 1 - ky.*ky./k2;
    % proj_free_div(1,1,:,:) = 0;
    % proj_free_div(ZM(1),:,:,:) = 0;
    % proj_free_div(:,ZM(2),:,:) = 0;
    d_kappa = 2*pi/sqrt(prod(model.grid.MX.* model.grid.dX));
    
    %% Masks associated with the rings of iso wave number
    k_idx_k = k;
    k_idx_k(PX(1)+1,:)=inf;
    k_idx_k(:,PX(2)+1)=inf;
    k_idx_k=k_idx_k(:);
    
    %% Wave number
    MX_kappa=model.grid.MX;
    % if model.sigma.ADSD_under_trace
    %     MX_kappa = MX_kappa - 2*model.sigma.local_grid.sslop_patch;
    % end
    M_kappa=min(MX_kappa);
    P_kappa= M_kappa/2;
    % d_kappa = 2*pi/sqrt(prod(MX_kappa.* model.grid.dX));
    kappa= d_kappa * ( 0:(P_kappa-1) ) ;
    
    d_kappa = kappa(2) - kappa(1);
    kappa_shift = [ kappa(2:end) kappa(end)+d_kappa ];
    idx_k = sparse( bsxfun(@le,kappa, k_idx_k ) );
    idx_k = idx_k & sparse( bsxfun(@lt,k_idx_k, kappa_shift ) );
    
    if (~bool_skip)
        %%  Anti-aliased grid
        % for non-linear operations, e.g. buoyancy advection w'*grad(b)
        % and non-homogeneous stochastic diffusion div(a*grad(b))
        if strcmp(model.grid.dealias_method, '2/3')
            % usual 2/3 rule: zero-out the 1/3 highest freqs
            maskx = double( abs(nx) < (2./3.)*PX(1) );
            masky = double( abs(ny) < (2./3.)*PX(2) );
        elseif strcmp(model.grid.dealias_method, 'exp')
            % see "New Numerical Results for the Surface Quasi-Geostrophic
            % Equation", Constantin et al., J. Sci. Comput. (2012).
            alpha = 36.;
            order = 19.;
            maskx = exp(-alpha*( (2./model.grid.MX(1)).*abs(nx) ).^order);
            masky = exp(-alpha*( (2./model.grid.MX(2)).*abs(ny) ).^order);
        elseif strcmp(model.grid.dealias_method, 'lowpass')
            % from legacy SQGMU code by V. Resseguier.
            maskx = fct_unity_approx5(model.grid.MX(1));
            masky = fct_unity_approx5(model.grid.MX(1));
        else
            error('SQGMU:init_grid_k:invalidParameter', ...
                'Unknown de-aliasing method "%s"', model.grid.dealias_method);
        end
        maskxy = maskx'*masky;
        maskxy(ZM(1),:) = 0.; %de-alias the single high freq
        maskxy(:,ZM(2)) = 0.;
        
        % Projection on the spce of free divergent function
        proj_free_div_aa(:,:,1,1,1) = 1 - kx.*kx./k2.*maskxy;
        proj_free_div_aa(:,:,1,1,2) = 0 - kx.*ky./k2.*maskxy;
        proj_free_div_aa(:,:,2,1,1) = 0 - ky.*kx./k2.*maskxy;
        proj_free_div_aa(:,:,2,1,2) = 1 - ky.*ky./k2.*maskxy;
        proj_free_div_aa(1,1,:,1,:) = 0;
        proj_free_div_aa(ZM(1),:,:,1,:) = 0;
        proj_free_div_aa(:,ZM(2),:,1,:) = 0;
        % proj_free_div_aa(:,:,1,1) = 1 - kx.*kx./k2.*maskxy;
        % proj_free_div_aa(:,:,1,2) = 0 - kx.*ky./k2.*maskxy;
        % proj_free_div_aa(:,:,2,1) = 0 - ky.*kx./k2.*maskxy;
        % proj_free_div_aa(:,:,2,2) = 1 - ky.*ky./k2.*maskxy;
        % proj_free_div_aa(1,1,:,:) = 0;
        % proj_free_div_aa(ZM(1),:,:,:) = 0;
        % proj_free_div_aa(:,ZM(2),:,:) = 0;
        
        %%  Anti-aliased grid at larger scale
        
        % for non-linear operations, e.g. buoyancy advection w'*grad(b)
        % and non-homogeneous stochastic diffusion div(a*grad(b))
        ratio_mask_LS = model.advection.Smag.dealias_ratio_mask_LS;
        if strcmp(model.grid.dealias_method, '2/3')
            % usual 2/3 rule: zero-out the 1/3 highest freqs
            maskx_LS = double( abs(nx) < (2./3.)*PX(1)*ratio_mask_LS );
            masky_LS = double( abs(ny) < (2./3.)*PX(2)*ratio_mask_LS );
        elseif strcmp(model.grid.dealias_method, 'exp')
            % see "New Numerical Results for the Surface Quasi-Geostrophic
            % Equation", Constantin et al., J. Sci. Comput. (2012).
            nx_normalized_LS = (2./model.grid.MX(1)).*abs(nx) / ratio_mask_LS;
            ny_normalized_LS = (2./model.grid.MX(2)).*abs(ny) / ratio_mask_LS;
            
            maskx_LS = exp(-alpha*( nx_normalized_LS ).^order);
            masky_LS = exp(-alpha*( ny_normalized_LS ).^order);
            maskx_LS = maskx_LS .* (nx_normalized_LS<=1);
            masky_LS = masky_LS .* (ny_normalized_LS<=1);
            %%
            n_normalized_LS = sqrt(nx_normalized_LS.^2 + (ny_normalized_LS').^2);
            mask_LS = exp(-alpha*( n_normalized_LS ).^(order) );
            mask_LS = mask_LS .* (n_normalized_LS<=1);
            %%
        elseif strcmp(model.grid.dealias_method, 'lowpass')
            % from legacy SQGMU code by V. Resseguier.
            maskx_LS = fct_unity_approx5(model.grid.MX(1)*ratio_mask_LS);
            masky_LS = fct_unity_approx5(model.grid.MX(1)*ratio_mask_LS);
        else
            error('SQGMU:init_grid_k:invalidParameter', ...
                'Unknown de-aliasing method "%s"', model.grid.dealias_method);
        end
        maskxy_LS = maskx_LS'*masky_LS;
        %%
        maskxy_LS = mask_LS;
        %%
        maskxy_LS(ZM(1),:) = 0.; %de-alias the single high freq
        maskxy_LS(:,ZM(2)) = 0.;
        
        
        %%  Anti-aliased grid at larger scale than sigma dBt
        if model.sigma.sto && ...
                ( model.sigma.hetero_modulation | model.sigma.hetero_energy_flux ...
                | model.sigma.hetero_modulation_V2 | model.sigma.hetero_modulation_Smag ) ...
                && model.sigma.hetero_energy_flux_prefilter
            %      ratio_pre_filter = model.sigma.kappamin_on_kappamax ;
            ratio_pre_filter = model.sigma.kappamin_on_kappamax / 1.5 ;
            %      ratio_pre_filter = model.sigma.kappamin_on_kappamax / 2 ;
            %     model.sigma.hetero_energy_flux
            if strcmp(model.grid.dealias_method, '2/3')
                % usual 2/3 rule: zero-out the 1/3 highest freqs
                maskx_sigma_hetero_energy_flux = double( abs(nx) < (2./3.)*PX(1)*ratio_pre_filter );
                masky_sigma_hetero_energy_flux = double( abs(ny) < (2./3.)*PX(2)*ratio_pre_filter );
            elseif strcmp(model.grid.dealias_method, 'exp')
                % see "New Numerical Results for the Surface Quasi-Geostrophic
                % Equation", Constantin et al., J. Sci. Comput. (2012).
                nx_normalized_sigma_hetero_energy_flux = (2./model.grid.MX(1)).*abs(nx) / ratio_pre_filter;
                ny_normalized_sigma_hetero_energy_flux = (2./model.grid.MX(2)).*abs(ny) / ratio_pre_filter;
                
                maskx_sigma_hetero_energy_flux = exp(-alpha*( nx_normalized_sigma_hetero_energy_flux ).^order);
                masky_sigma_hetero_energy_flux = exp(-alpha*( ny_normalized_sigma_hetero_energy_flux ).^order);
                maskx_sigma_hetero_energy_flux = maskx_sigma_hetero_energy_flux .* (nx_normalized_sigma_hetero_energy_flux<=1);
                masky_sigma_hetero_energy_flux = masky_sigma_hetero_energy_flux .* (ny_normalized_sigma_hetero_energy_flux<=1);
                %%
                n_normalized_sigma_hetero_energy_flux = sqrt(nx_normalized_sigma_hetero_energy_flux.^2 + (ny_normalized_sigma_hetero_energy_flux').^2);
                mask_sigma_hetero_energy_flux = exp(-alpha*( n_normalized_sigma_hetero_energy_flux ).^(order) );
                mask_sigma_hetero_energy_flux = mask_sigma_hetero_energy_flux .* (n_normalized_sigma_hetero_energy_flux<=1);
                %%
            elseif strcmp(model.grid.dealias_method, 'lowpass')
                % from legacy SQGMU code by V. Resseguier.
                maskx_sigma_hetero_energy_flux = fct_unity_approx5(model.grid.MX(1)*ratio_pre_filter);
                masky_sigma_hetero_energy_flux = fct_unity_approx5(model.grid.MX(1)*ratio_pre_filter);
            else
                error('SQGMU:init_grid_k:invalidParameter', ...
                    'Unknown de-aliasing method "%s"', model.grid.dealias_method);
            end
            maskxy_sigma_hetero_energy_flux = maskx_sigma_hetero_energy_flux'*masky_sigma_hetero_energy_flux;
            %%
            maskxy_sigma_hetero_energy_flux = mask_sigma_hetero_energy_flux;
            %%
            maskxy_sigma_hetero_energy_flux(ZM(1),:) = 0.; %de-alias the single high freq
            maskxy_sigma_hetero_energy_flux(:,ZM(2)) = 0.;
        end
        
        
        
        %%  Filtering at much larger scale than sigma dBt
        if model.sigma.sto && ...
                ( model.sigma.hetero_modulation | model.sigma.hetero_energy_flux ...
                | model.sigma.hetero_modulation_V2 | model.sigma.hetero_modulation_Smag ) ...
                && model.sigma.hetero_energy_flux_postfilter
            %      ratio_post_filter = model.sigma.kappamin_on_kappamax ;
            ratio_post_filter = model.sigma.kappamin_on_kappamax / 2 ;
            %      ratio_post_filter = model.sigma.kappamin_on_kappamax / 4 ;
            %     model.sigma.hetero_energy_flux
            if strcmp(model.grid.dealias_method, '2/3')
                % usual 2/3 rule: zero-out the 1/3 highest freqs
                maskx_sigma_hetero_energy_flux = double( abs(nx) < (2./3.)*PX(1)*ratio_post_filter  );
                masky_sigma_hetero_energy_flux = double( abs(ny) < (2./3.)*PX(2)*ratio_post_filter );
            elseif strcmp(model.grid.dealias_method, 'exp')
                % see "New Numerical Results for the Surface Quasi-Geostrophic
                % Equation", Constantin et al., J. Sci. Comput. (2012).
                nx_normalized_sigma_hetero_energy_flux = (2./model.grid.MX(1)).*abs(nx) / ratio_post_filter;
                ny_normalized_sigma_hetero_energy_flux = (2./model.grid.MX(2)).*abs(ny) / ratio_post_filter;
                
                maskx_sigma_hetero_energy_flux = exp(-alpha*( nx_normalized_sigma_hetero_energy_flux ).^order);
                masky_sigma_hetero_energy_flux = exp(-alpha*( ny_normalized_sigma_hetero_energy_flux ).^order);
                maskx_sigma_hetero_energy_flux = maskx_sigma_hetero_energy_flux .* (nx_normalized_sigma_hetero_energy_flux<=1);
                masky_sigma_hetero_energy_flux = masky_sigma_hetero_energy_flux .* (ny_normalized_sigma_hetero_energy_flux<=1);
                %%
                n_normalized_sigma_hetero_energy_flux = sqrt(nx_normalized_sigma_hetero_energy_flux.^2 + (ny_normalized_sigma_hetero_energy_flux').^2);
                mask_sigma_hetero_energy_flux = exp(-alpha*( n_normalized_sigma_hetero_energy_flux ).^(order) );
                mask_sigma_hetero_energy_flux = mask_sigma_hetero_energy_flux .* (n_normalized_sigma_hetero_energy_flux<=1);
                %%
            elseif strcmp(model.grid.dealias_method, 'lowpass')
                % from legacy SQGMU code by V. Resseguier.
                maskx_sigma_hetero_energy_flux = fct_unity_approx5(model.grid.MX(1)*ratio_post_filter);
                masky_sigma_hetero_energy_flux = fct_unity_approx5(model.grid.MX(1)*ratio_post_filter);
            else
                error('SQGMU:init_grid_k:invalidParameter', ...
                    'Unknown de-aliasing method "%s"', model.grid.dealias_method);
            end
            maskxy_sigma_hetero_energy_flux_post = maskx_sigma_hetero_energy_flux'*masky_sigma_hetero_energy_flux;
            %%
            maskxy_sigma_hetero_energy_flux_post = mask_sigma_hetero_energy_flux;
            %%
            maskxy_sigma_hetero_energy_flux_post(ZM(1),:) = 0.; %de-alias the single high freq
            maskxy_sigma_hetero_energy_flux_post(:,ZM(2)) = 0.;
        end
        
        % figure(1);plot(kx(:,1),maskx,'.');
        % figure(2);plot(kx(:,1),maskx_LS,'.');
        % figure(1);plot(ky(1,:),masky,'.');
        % figure(2);plot(ky(1,:),masky_LS,'.');
                
        %% Band pass filter for forcing
        maskf = double (k>=1.2*10.0^(-5) & k<=1.4*10.0^(-5));         % SQG double (kmod>=5.0*10.0^(-6) & kmod<=7.0*10.0^(-6));
        maskfc =double (~(k>=1.2*10.0^(-5) & k<=1.4*10.0^(-5)));  % SQG double (~(kmod>=5.0*10.0^(-6) & kmod<=7.0*10.0^(-6)));
        
    end
end
%% Save
% the "unstable" grid
model.grid.k.ZM = ZM; %indices of the single mode to force to zero
model.grid.k.kx_ref = kx_ref;
model.grid.k.ky_ref = ky_ref;
model.grid.k.dK = dK ;
model.grid.k.Period_k = Period_k;
if ~ bool_lowCost
    model.grid.k.kx = kx;
    model.grid.k.ky = ky;
    % model.grid.k.d_kx = d_kx;
    % model.grid.k.d_ky = d_ky;
    model.grid.k.k2 = k2;
    model.grid.k.k = k;
    if ~ bool_skip
        model.grid.k.on_k = over_k;
        model.grid.k.kx_over_ksqr = kx.*(over_k.^2);
        model.grid.k.ky_over_ksqr = ky.*(over_k.^2);
        model.grid.k.on_k2 = on_k2;
    end
    model.grid.k.proj_free_div = proj_free_div;
    model.grid.k.d_kappa = d_kappa;
    model.grid.k.idx_k = idx_k;
    if (~bool_skip)
        % the "anti-aliased" grid
        model.grid.k_aa.ikx = 1i.*maskxy.*kx; %precomputations for de-aliased gradient
        model.grid.k_aa.iky = 1i.*maskxy.*ky;
        model.grid.k_aa.k2 = (maskxy.^2).*k2;
        %model.grid.k_aa.k2 = maskxy.*k2;
        model.grid.k_aa.mask = maskxy;
        model.grid.k_aa.proj_free_div = proj_free_div_aa;
        % the "anti-aliased" grid at larger scales
        model.grid.k_aa_LS.ikx = 1i.*maskxy_LS.*kx; %precomputations for de-aliased gradient
        model.grid.k_aa_LS.iky = 1i.*maskxy_LS.*ky;
        model.grid.k_aa_LS.k2 = (maskxy_LS.^2).*k2;
        %model.grid.k_aa.k2 = maskxy.*k2;
        model.grid.k_aa_LS.mask = maskxy_LS;
        if model.sigma.sto && ...
                ( model.sigma.hetero_modulation | model.sigma.hetero_energy_flux ...
                | model.sigma.hetero_modulation_V2 | model.sigma.hetero_modulation_Smag ) ...
                && model.sigma.hetero_energy_flux_prefilter
            %     model.grid.k_aa_sigma_hetero_energy_flux.ikx = 1i.*maskxy_sigma_hetero_energy_flux.*kx; %precomputations for de-aliased gradient
            %     model.grid.k_aa_sigma_hetero_energy_flux.iky = 1i.*maskxy_sigma_hetero_energy_flux.*ky;
            %     model.grid.k_aa_sigma_hetero_energy_flux.k2 = (maskxy_sigma_hetero_energy_flux.^2).*k2;
            %     %model.grid.k_aa.k2 = maskxy.*k2;
            model.grid.k_aa_sigma_hetero_energy_flux.mask = maskxy_sigma_hetero_energy_flux;
        end
        if model.sigma.sto && ...
                ( model.sigma.hetero_modulation | model.sigma.hetero_energy_flux ...
                | model.sigma.hetero_modulation_V2 | model.sigma.hetero_modulation_Smag ) ...
                && model.sigma.hetero_energy_flux_prefilter
            %     model.grid.k_aa_sigma_hetero_energy_flux.ikx = 1i.*maskxy_sigma_hetero_energy_flux.*kx; %precomputations for de-aliased gradient
            %     model.grid.k_aa_sigma_hetero_energy_flux.iky = 1i.*maskxy_sigma_hetero_energy_flux.*ky;
            %     model.grid.k_aa_sigma_hetero_energy_flux.k2 = (maskxy_sigma_hetero_energy_flux.^2).*k2;
            %     %model.grid.k_aa.k2 = maskxy.*k2;
            model.grid.k_aa_sigma_hetero_energy_flux_post.mask = maskxy_sigma_hetero_energy_flux_post;
        end
        
    end
end
end

function t = fct_unity_approx5(N_t)
% XP must be a 2 x n matrix
% the result is a vector of size n
%

slop_size_ratio=6;

t=ones(1,N_t);
P_t=N_t/2;
sslop=ceil(N_t/slop_size_ratio);
t((P_t-sslop+1):P_t)= (-tanh(-3 + 6/(sslop-1)*(0:(sslop-1)) ) +1)/2;
t((P_t+2):(P_t+1+sslop))= (tanh(-3 + 6/(sslop-1)*(0:(sslop-1)) )+1)/2;

t(P_t+1)=0;

end


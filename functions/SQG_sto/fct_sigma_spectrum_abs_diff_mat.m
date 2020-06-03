function [breve_sigma,f_sigma,trace_a,....
    slope_w_a_estim,TKE,...
    km_LS] = ...
    fct_sigma_spectrum_abs_diff_mat(model,ft_w,bool_plot,day,...
    ft_w_resolved)
%
%% Outputs :
% - breve_sigma : fft of the v' kernel ( \breve{\sigma} )
%         size : Mx x My x 2 x N_ech
%         breve_sigma = fft of the v' kernel ( \breve{\sigma} )
% - f_sigma : Fourier transform of the associted streamfunction
% - trace_a : total absolute diffusivity of v'
%   (generally used to set the CFL)
% - slope_w_a_estim : 
%   KE spectrum (or ADSD) slope estimated from the large-scale velocity ft_w
% - TKE : total absolute diffusivity of v'
% - km_LS : 
%   wave-number where the estimated ADSD (of the form C*k^r)
%   and the large-scale velocity ADSD touch.
%
%% Inputs :
% - model : structure which gathers all simulation information
% - ft_w : fft of the large-scale velocity
%         size : Mx x My x 2 x N_patch
% - bool_plot : boolean to enable intermediate plots
%         size : 1
% - day : index of the current day (for plots)
%         size : 1
% - ft_w_resolved : fft of the fully-resolved velocity
%     (optional : for plots comparison only)
%         size : Mx x My x 2 x 1
%

% Is the fully-resolved velocity ft_w_resolved available ?
bool_resolved = (nargin >= 5) && (~any(isnan(ft_w_resolved(:))));

% bool_plot = true;
% day = '0';

% label of the chosen particle for plot
id_part = 1;

% Spectrum slope
slope_ref = model.sigma.slope_sigma_ref;
% ADSD slope
slope_ref_a = (slope_ref-3)/2;

%% Compute KESD or ADSD
% KESD or ADSD of large-scale velocity
[spectrum_w_a, kappa, d_kappa, k] = ...
    fct_spec_for_sigma(model, model.grid, ft_w, model.advection.N_ech );
if bool_resolved
    % KESD or ADSD of fully-resolved velocity
    [spectrum_w_a_resolved, kappa_resolved] = ...
        fct_spec_for_sigma(model, model.grid_resolved, ft_w_resolved,1 );
end

%% Estimation of slope of the spectral absolute diffusivity of
% % the large-scale velocity

% Compensated absolute diffusity by scale
spectrum_w_a_comp = bsxfun( @times, kappa.^(-slope_ref_a) , ...
    spectrum_w_a );

% Get the large-scale "length scale" km_LS and the spectral window for the
% linear regression
if ~ model.sigma.sto
    model.sigma.kappamin_on_kappamax = 1/2;
    model.sigma.kappaLS_on_kappamax = 1/8;
end
threshold_k = model.sigma.kappamin_on_kappamax_estim_slope;
% threshold_k = model.sigma.kappamin_on_kappamax;
spectrum_w_a_comp_cut = spectrum_w_a_comp( kappa < kappa(end)*threshold_k , :);

if model.sigma.estim_k_LS
    threshold_k_LS = model.sigma.kappaLS_on_kappamax;
    [~,i_first] = max(spectrum_w_a_comp_cut( [ false; ...
        ( kappa(2:end) < kappa(end)*threshold_k_LS )] , :) ,[], 1 );
else
    i_first = 2;
end
i_first = i_first +1;
km_LS = repmat( kappa(i_first), [ 1 model.advection.N_ech]);

mask_iii_k_LS = (1:size(spectrum_w_a_comp_cut,1))' ...
    * ones(1,model.advection.N_ech);
mask_iii_k_LS = bsxfun( @ge,  mask_iii_k_LS , i_first );
spectrum_w_a_comp_for_estim = mask_iii_k_LS .* spectrum_w_a_comp_cut;
kappa_cut = kappa( kappa<kappa(end)*threshold_k );
kkk = bsxfun( @times, mask_iii_k_LS , kappa_cut );
mask_iii_k_LS_one_value = (1:size(spectrum_w_a_comp_cut,1))' ...
    * ones(1,model.advection.N_ech);
mask_iii_k_LS_one_value = bsxfun( @eq,  mask_iii_k_LS_one_value , i_first );
offset_w_a_comp = sum( mask_iii_k_LS_one_value .* spectrum_w_a_comp_cut ,1);

% Linear regression
iii_reliable =~ (isinf(abs(spectrum_w_a_comp_for_estim))|...
    isnan(spectrum_w_a_comp_for_estim)|...
    (spectrum_w_a_comp_for_estim/max(spectrum_w_a_comp_for_estim(:)) ...
    <=eps)|...
    (kkk>kappa(end)*threshold_k));

% Centering in the point ( km_LS, spectrum_w_a_comp(km_LS) )
logspectrum_w_a_comp_for_estim_centered = bsxfun( @plus, ...
    log10(spectrum_w_a_comp_for_estim) , - log10( offset_w_a_comp) );
logkkk_centered = bsxfun( @plus, log10(kkk)  , - log10( km_LS ) );

% Set to zero the unwanted points
siz = size(iii_reliable);
iii_reliable = iii_reliable(:);
logspectrum_w_a_comp_for_estim_centered( ~ iii_reliable ) = 0; % P_kappa_cut*N_ech
logkkk_centered( ~ iii_reliable ) = 0; % P_kappa*N_ech
logspectrum_w_a_comp_for_estim_centered = reshape ( ...
    logspectrum_w_a_comp_for_estim_centered, siz ); % P_kappa_cut x N_ech
logkkk_centered = reshape ( ...
    logkkk_centered, siz ); % P_kappa_cut x N_ech

% Linear regression
beta_num = sum( ...
    logkkk_centered .* logspectrum_w_a_comp_for_estim_centered , 1);
beta_den = sum( ...
    logkkk_centered .* logkkk_centered , 1);
beta = beta_num ./ beta_den; % 1 x N_ech

slope_w_a_estim = beta + slope_ref_a;

% Prevent some possible unstable behaviors of this parametrisation:
% The maximum allowed slope corresponds to a velocity white in space
slope_w_a_estim = min( [zeros([1 model.advection.N_ech]) ; ...
    slope_w_a_estim]);

%% Omnidirectional Spectrum

% KESD or ADSD with estimated slope
reference_spectrum_a_estim = bsxfun(@power, kappa(2:end), ...
    slope_w_a_estim ) ;

% Absolute diffusivity by scale with theoretical slope
if bool_plot
    reference_spectrum_a = kappa(2:end) .^ slope_ref_a ;
    reference_spectrum_a = repmat( reference_spectrum_a , ...
        [ 1 model.advection.N_ech ]);
end

%% Multiplicative constant C
spectrum_w_a_km_LS = spectrum_w_a(i_first,:);
reference_spectrum_a_estim_km_LS = reference_spectrum_a_estim(i_first-1,:);
if bool_plot
    reference_spectrum_a_km_LS = reference_spectrum_a(i_first-1,:);
end

% Apply multiplicative constant
mult_offset_spectrum_a_estim = spectrum_w_a_km_LS ...
    ./ reference_spectrum_a_estim_km_LS;
reference_spectrum_a_estim = bsxfun(@times, ...
    mult_offset_spectrum_a_estim , reference_spectrum_a_estim );

if bool_plot
    mult_offset_spectrum_a = spectrum_w_a_km_LS ...
        ./ reference_spectrum_a_km_LS ;
    reference_spectrum_a = bsxfun(@times, ...
        mult_offset_spectrum_a , reference_spectrum_a );
end

%% Residual ADSD or KESD
residual_spectrum_a = ...
    reference_spectrum_a_estim - spectrum_w_a(2:end,:);

% Clean it
siz = size(residual_spectrum_a);
residual_spectrum_a = residual_spectrum_a(:);
residual_spectrum_a(residual_spectrum_a<0)=0;
residual_spectrum_a = reshape( residual_spectrum_a, siz);

residual_spectrum_a(1:(i_first-1),:)=0;
% residual_spectrum_a(1:2,:)=0;

%% Bi-directionnal KESD or ADSD of v' and streamfunction

% Band-pass filter
k_inf = kappa(end);% Largest wave number
k0 = model.sigma.kappamin_on_kappamax * k_inf;% Smallest wave number
idx1 = (kappa(2:end) <= k0);
idx3 = (kappa(2:end) > k_inf);
idx2 = ~ (idx1 | idx3);

residual_spectrum_a(idx1 | idx3,:)=0;

if model.sigma.band_pass_filter
    unit_approx = fct_unity_approx_(sum(idx2));
    residual_spectrum_a(idx2,:) = bsxfun(@times, unit_approx' , ...
        residual_spectrum_a(idx2,:) ) ;
end

% To remove the 2 pi which appear when we integrate the spectrum over the
% wave-vector angles and add the (2 pi)^2 which appear when we go from k to
% 2*pi*k
% And discretisation to go from continuous to discrete Fourier transform
residual_spectrum_a = (2*pi/prod(model.grid.dX)) * residual_spectrum_a;

% Division by k^2 to get the spectrum of the streamfunction
% And from omnidirectional spectrum to Fourier tranform square modulus
% Division by k in dimension 2
residual_spectrum_a = bsxfun(@times, 1 ./ ( kappa(2:end).^3) , residual_spectrum_a );

%% FFT(kernel) v'
% From square modulus to modulus
f_sigma = sqrt( residual_spectrum_a );

% From 1D function to 2D function
f_sigma = interp1(kappa,[zeros(1,model.advection.N_ech); f_sigma],k);

% Cleaning
f_sigma(k<=k0,:)=0;
f_sigma(k>k_inf,:)=0;
f_sigma=reshape(f_sigma,[model.grid.MX 1 model.advection.N_ech]);

% Antialiasing
PX=model.grid.MX/2;
f_sigma(PX(1)+1,:,:,:)=0;
f_sigma(:,PX(2)+1,:,:)=0;

% Orthogonal gradient (from streamfunction to velocity)
breve_sigma(:,:,1,:)= 1i * bsxfun(@times, - model.grid.k.ky , f_sigma );
breve_sigma(:,:,2,:)= 1i * bsxfun(@times,  + model.grid.k.kx , f_sigma );

%% Bi-directional KESD or ADSD of sigma dBt
ft_sigma=abs(breve_sigma).^2;
ft_sigma=sum(ft_sigma,3);

%% Calcul of energy
% Division by prod(model.grid.MX) because of the Parseval theorem for
% discrete Fourier transform
% Division by prod(model.grid.MX) again in order to the integration
% of the spectrum over the wave number yields the energy of the
% buoyancy averaged (not just integrated) over the space
% One has to multiply by prod(model.grid.MX) because of the variance of
% the white in space noise
TKE = 1/prod(model.grid.MX) * sum(sum(ft_sigma,2),1);
TKE = permute( TKE, [ 1 4 3 2]);
% trace_a = 1/(prod(model.grid.MX)*model.advection.N_ech) * sum(ft_sigma(:));
% % trace_a = 1/(prod(model.grid.MX) * sum(ft_sigma(:));
trace_a = TKE ;

%% Plots
if bool_plot
    vect_pcl = id_part;
    if bool_resolved
        spectrum_w_a_resolved_plot=spectrum_w_a_resolved;
    end
    w=real(ifft2(ft_w));
    for id_part = vect_pcl
        % Choose one realization
        ft_sigma_plot = ft_sigma(:,:,:,id_part);
        reference_spectrum_a_plot=reference_spectrum_a(:,id_part);
        spectrum_w_a_plot=spectrum_w_a(:,id_part);
        reference_spectrum_a_estim_plot=reference_spectrum_a_estim(:,id_part);
        km_LS_plot = km_LS(id_part);
        
        % Compute the omnidirectional spectrum of sigma dBt
        spectrum_a_sigma = model.grid.k.idx_k' * ft_sigma_plot(:);
        
        % Influence of the complex brownian variance
        spectrum_a_sigma = prod(model.grid.MX)*spectrum_a_sigma;
        
        % Division by prod(model.grid.MX) because of the Parseval theorem for
        % discrete Fourier transform
        % Division by prod(model.grid.MX) again in order to the integration
        % of the spectrum over the wave number yields the energy of the
        % buoyancy averaged (not just integrated) over the space
        spectrum_a_sigma = 1/prod(model.grid.MX)^2 * spectrum_a_sigma;
        
        % Energy test
        
        % Division by the wave number step
        spectrum_a_sigma = spectrum_a_sigma / d_kappa;
        
        % Absolute diffusivity of the small-scale velocity (for test)
        trace_a_from_spectrum = d_kappa * sum(spectrum_a_sigma(:));
        
        % Plots
        
        % % Make the plots appear at the same level thant the large-scale velocity spectrum
        % spectrum_a_sigma_plot = spectrum_a_sigma * mult_offset;
        spectrum_a_sigma_plot = spectrum_a_sigma;
        
        taille_police = 12;
        X0 = [1 2];
        
        figure12=figure(12);
        hold off;
        widthtemp = 12;
        heighttemp = 6;
        set(figure12,'Units','inches', ...
            'Position',[X0 widthtemp heighttemp], ...
            'PaperPositionMode','auto');
        loglog(kappa(2:end),reference_spectrum_a_plot,'k')
        hold on;
        loglog(kappa(2:end),reference_spectrum_a_estim_plot,'k--')
        loglog(kappa(2:end),spectrum_w_a_plot(2:end))
        if bool_resolved
            loglog(kappa_resolved(2:end),spectrum_w_a_resolved_plot(2:end),'b--')
        end
        loglog(kappa(2:end),spectrum_a_sigma_plot(2:end),'r.-')
        
        ax=axis;
        ax(4)=max([spectrum_w_a_plot; reference_spectrum_a_estim_plot ; ...
            reference_spectrum_a_plot ; spectrum_a_sigma_plot]);
        ax(4)=ax(4)*2;
        ax(3)=(kappa(2)/kappa(end))*min([max(spectrum_w_a_plot); ...
            max(reference_spectrum_a_estim_plot); max(reference_spectrum_a_plot); ...
            max(spectrum_a_sigma_plot)]);
        ax(3) = min( [ax(3) min(spectrum_a_sigma_plot(spectrum_a_sigma_plot>0)) ...
            min(reference_spectrum_a_plot) min(reference_spectrum_a_estim_plot)]);
        ax(1:2)=kappa([2 end]);
        if ax(4)>0
            axis(ax)
        end
        ax = axis;
        loglog(km_LS_plot*[1 1],...
            [min(reference_spectrum_a_estim_plot) ax(4)],'k--')
        loglog(...
            model.sigma.kappamin_on_kappamax_estim_slope * kappa(end)*[1 1],...
            [min(reference_spectrum_a_estim_plot) ax(4)],'k--')
        loglog(model.sigma.kappamin_on_kappamax * kappa(end)*[1 1],...
            [min(reference_spectrum_a_estim_plot) ax(4)],'k-.')
        hold off
        set(gca,'XGrid','on','XTickMode','manual');
        width = 9;
        height = 3;
        set(figure12,'Units','inches', ...
            'Position',[X0 width height], ...
            'PaperPositionMode','auto');
        set(gca,'YGrid','on')
        set(gca,...
            'Units','normalized',...
            'FontUnits','points',...
            'FontWeight','normal',...
            'FontSize',taille_police,...
            'FontName','Times')
        xlabel('$\kappa \bigl ( rad.m^{-1} \bigr )$',...
            'FontUnits','points',...
            'FontWeight','normal',...
            'FontSize',taille_police,...
            'interpreter','latex',...
            'FontName','Times')
        ylabel('$A(\kappa)$',...
            'FontUnits','points',...
            'interpreter','latex',...
            'FontSize',taille_police,...
            'FontName','Times')
        title('ADSD for $w$ and $\sigma dB_t$',...
            'FontUnits','points',...
            'FontWeight','normal',...
            'interpreter','latex',...
            'FontSize',12,...
            'FontName','Times')
        
    end
    % Save plot
    drawnow
    folder_simu = model.folder.folder_simu;
    eval( ['print -depsc ' folder_simu '/AbsDiffByScale_sigma_dB_t/'...
        day '.eps']);
    %     keyboard;
    
end
end

function t = fct_unity_approx_(N_t)
% Approximation of unity
%

nx = 1:N_t;
nx=nx-mean(nx);
alpha = 36.;
order = 19.;
t = exp(-alpha*( (2./N_t).*abs(nx) ).^order);

end
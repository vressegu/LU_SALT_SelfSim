%% Default parameters for the deterministic subgrid tensor
% Do not change

% Hyper-viscosity order
if HV.bool
    HV.order=8;
end

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
    % Ratio between the Shanon resolution and filtering frequency used to
    % filter the heterogenous diffusion coefficient
    Smag.dealias_ratio_mask_LS = 1/8;
    
    Smag.kappamax_on_kappad = 0;
end

dealias_method = 'exp';
% [WIP] Method for mandatory de-aliasing of non-linear terms in
% pseudospectral codes (advection and non-homogeneous stochastic diffusion)
% - 'lowpass': same as in SQGMU 1;
% - '2/3': the classical 2/3 rule (NB: produces Gibb's oscillations);
% - 'exp': high-order exponential filter (Constantin et al., J. Sci.
%   Comput. (2012)).
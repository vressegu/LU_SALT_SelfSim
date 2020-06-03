%% Initialize hyperviscosity
if model.advection.Lap_visco.bool | model.advection.HV.bool
    % Root mean square Lyapunov
    [dxUx,dyUx] = gradient(w(:,:,1)',model.grid.x_ref,model.grid.y_ref);
    [dxUy,dyUy] = gradient(w(:,:,2)',model.grid.x_ref,model.grid.y_ref);
    s=(dxUx+dyUy)/2;
    d=(dyUx+dxUy)/2;
    lambda = sqrt(s.^2+d.^2);
    model.advection.lambda_RMS = sqrt(mean(lambda(:).^2));
    clear s d dxUx dxUy dyUx dyUy lambda
    
    % Order of the hyperviscosity
    if model.advection.HV.bool
        % model.advection.HV.order=8;
    else
        model.advection.HV.order=2;
    end
    
    % Hyperviscosity coefficient
    model.advection.HV.val= ...
        40 * model.advection.lambda_RMS * ...
        (mean(model.grid.dX)/pi)^model.advection.HV.order;
    
    model.advection.HV.val
    % model.advection.HV.val*(model.grid.MX(1)^model.advection.HV.order)
    model.advection.HV.val/(mean(model.grid.dX)^model.advection.HV.order)
    model.advection.lambda_RMS
    
    % if ~model.advection.HV.bool % DNS
    if model.advection.Lap_visco.bool &&  strcmp(model.dynamics,'2D')
        if isfield(model.advection, 'forcing') && model.advection.forcing.bool
            model.carct.L_caract = 1/ ...
                sqrt(mean( ( [ model.advection.forcing.Lx ...
                model.advection.forcing.Ly ] ...
                .\ model.advection.forcing.freq_f ).^2));
            warning('Re should be computed with the amplitude of forcing!');
        else
            fft_grad_b = fct_grad(model,fft_b);
            model.carct.L_caract = sqrt( ...
                sum(abs(fft_b(:)).^2) / sum(abs(fft_grad_b(:)).^2) );
        end
        model.carct.U_caract = model.odg_b * model.carct.L_caract;
        Re = model.carct.U_caract * model.carct.L_caract ...
            / model.advection.HV.val
        model.carct.Re = Re;
        model.carct.l_Kolmogorov = Re^(-1/2) * model.carct.L_caract;
        if model.carct.l_Kolmogorov < min(model.grid.dX)
            error('The simulation is under resolved');
        end
    end
    
    if model.advection.Smag.bool
        model.advection.Smag.coef_Smag = ...
            ( model.advection.Smag.kappamax_on_kappad ...
            * sqrt(prod(model.grid.dX))/pi ) ...
            ^ (5*model.advection.HV.order/4 -1/2);
        
        
        % Heterogeneous HV or diffusivity/viscosity coefficient
        if model.advection.Lap_visco.bool
            % Heterogeneous dissipation coefficient
            coef_diff_aa = fct_coef_diff(model, fft_b);
            % Coefficient coef_Smag to target a specific diffusive scale
            coef_diff_aa = model.advection.Smag.coef_Smag * coef_diff_aa ;
            % Maximum dissipation coefficient
            model.advection.HV.maxVal = max(coef_diff_aa(:));
        elseif model.advection.HV.bool
            % Heterogeneous HV coefficient
            coef_HV_aa = fct_coef_HV(model, fft_b);
            
            % Maximum HV coefficient
            model.advection.HV.maxVal = max(coef_HV_aa(:));
        else
            error('Unknown deterministic subgrid tensor');
        end
    else
        model.advection.HV.maxVal = model.advection.HV.val;
    end
else
    model.advection.HV.val= 0;
    model.advection.HV.order = 2;
    model.advection.HV.maxVal = 0;
end
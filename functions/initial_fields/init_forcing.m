%% Initialize forcing

if isfield(model.advection, 'forcing') && model.advection.forcing.bool
    model.advection.forcing.Lx = model.grid.MX(1) * model.grid.dX(1);
    model.advection.forcing.Ly = model.grid.MX(2) * model.grid.dX(2);
    [X_forcing,Y_forcing]=ndgrid(model.grid.x,model.grid.y);
    
    switch model.dynamics
        case 'SQG'
            U_caract = ...
                model.odg_b / model.physical_constant.buoyancy_freq_N;
            U_caract = U_caract /5;
        case '2D'
            U_caract = ...
                model.odg_b * model.advection.forcing.Ly;
            U_caract = U_caract /20;
            U_caract = U_caract /4;
        otherwise
            error('Unknown type of dynamics');
    end
    U_caract = U_caract /3;
    
    model.advection.forcing.on_T =  U_caract / model.advection.forcing.Ly;
    
    model.advection.forcing.ampli_forcing = ...
        model.advection.forcing.ampli_forcing * model.odg_b ;
    ampli_scale = 1;
    
    switch model.advection.forcing.forcing_type
        case 'Kolmogorov'
            F = ampli_scale * ...
                model.advection.forcing.ampli_forcing * ...
                cos( 2 * pi / model.advection.forcing.Lx ...
                * model.advection.forcing.freq_f(1) * X_forcing ...
                + 2 * pi / model.advection.forcing.Ly ...
                * model.advection.forcing.freq_f(2) * Y_forcing );
            F = F / model.advection.forcing.on_T;
        case {'Spring','Hetero_Spring'}
            F = ampli_scale * ...
                model.advection.forcing.ampli_forcing ...
                * sin( 2 * pi / model.advection.forcing.Lx ...
                * model.advection.forcing.freq_f(1) * X_forcing ) ...
                .* sin( 2 * pi / model.advection.forcing.Ly ...
                * model.advection.forcing.freq_f(2) * Y_forcing );
    end
    
    
    if strcmp(model.advection.forcing.forcing_type, 'Hetero_Spring')
        model.advection.forcing.F = F;
    else
        model.advection.forcing.F = fft2(F);
    end
    
    if strcmp(model.type_data, 'Zero')
        fft_w = SQG_large_UQ(model,  ...
            model.odg_b / model.advection.forcing.ampli_forcing ...
            * model.advection.forcing.F);
        if strcmp( model.advection.forcing.forcing_type,'Kolmogorov')
            fft_w = fft_w * model.advection.forcing.on_T;
        end
        w = real(ifft2(fft_w));
    end
    
    clear Lx Ly X_forcing Y_forcing on_T U_caract ampli_scale
end
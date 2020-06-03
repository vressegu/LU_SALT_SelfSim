%% Additional parameters for sigma
% Do not change

if ~ sigma.sto
    % If the simulation is deterministic, a_H = 0 and only one simulation
    % is performed
    sigma.k_c = inf; % And then a_H = 0
    N_ech=1;
    plot_moments = false;
else
    % Plot buoyancy satistical moments
    plot_moments = true;
    
    % Rate between the smallest and the largest wave number of sigma dBt
    switch sigma.type_spectrum
        case {'SelfSim_from_LS','EOF','Euler_EOF'}
            pre_estim_slope=1e-1;
            pre_5 = 5e-2;
            sigma.kappamin_on_kappamax = ...
                (log(1-pre_5)/log(pre_estim_slope))^(2/HV.order);
            sigma.kappamin_on_kappamax_estim_slope = ...
                (log(1-pre_estim_slope)/log(pre_estim_slope))...
                ^(2/HV.order);
            sigma.kappaLS_on_kappamax = 1/8;
        otherwise
            switch resolution
                case  128
                    sigma.kappamin_on_kappamax = 1/2;
                case 64
                    sigma.kappamin_on_kappamax = 1/3;
                otherwise
                    error('unknown');
            end
    end
    
    % Rate between the largest wave number of sigma dBt and the largest wave
    % number of the simulation
    sigma.kappamax_on_kappaShanon = 1;
end
%% ADSD comparaison of Self Sim and EOF methods

if model.sigma.plot_ADSD_SelfSim_EOF_bool
    day_num = (floor(time/24/3600));
    %     if (first_time == time) && ...
    if (day_num == 105) && plt_first_time && ...
            strcmp(model.sigma.type_spectrum,'SelfSim_from_LS')
        plt_first_time = false;
        color_sexy(1,:) = [0.8500, 0.3250, 0.0980];
        color_sexy(2,:) = [0.9290, 0.6940, 0.1250];
        color_sexy(3,:) = [0.4940, 0.1840, 0.5560];
        color_sexy(4,:) = [0.4660, 0.6740, 0.1880] 	;
        color_sexy(5,:) = [0.3010, 0.7450, 0.9330];
        color_sexy(6,:) = [0.6350, 0.0780, 0.1840];
        LineStyle_ =[ '-' ];
        Marker_ = {'none','o','+','*','x','s'};
        
        %     plot_abs_diff_from_sigma_postprocess(model,fft2(sigma_dBt_on_sq_dt))
        fct_sigma_spectrum_abs_diff_postprocess(model,fft2(w),true,'0',...
            color_sexy(6,:));
        
        % Load precomputed EOFs and correspond variance tensor
        model.folder.folder_EOF = [ pwd '/images/SQG_MU_HV_4/' ...
            'type_spectrum_sigma_EOF/' ...
            'disym_Vortices_forced_turb_Spring/64x64/folder_EOF'];
        load([ model.folder.folder_EOF '/EOF.mat'],'EOF');
        EOF = permute(EOF,[1 2 3 5 4]);
        sigma = EOF; clear EOF;
        
        nb_EOF_v = [8000 2000 200 20 2];
        for k=1:length(nb_EOF_v)
            nb_EOF = nb_EOF_v(k);
            sigma = sigma(:,:,:,:,1:nb_EOF);
            sigma_dBt_on_sq_dt = sum( sigma .* ...
                randn( [ 1 1 1 10 nb_EOF ]) , 5);
            plot_abs_diff_from_sigma_postprocess_add(model,...
                fft2(sigma_dBt_on_sq_dt), color_sexy(end-k,:), ...
                LineStyle_, Marker_{k+1} );
        end
        pl = findobj(gca,'Type','line');
        figure10=figure(10);
        width = 9;
        height = 3;
        set(figure10,'Units','inches', ...
            'Position',[0 0  width height], ...
            'PaperPositionMode','auto');
        axP = get(gca,'Position');
        lgd=legend(pl([10 9 5 4 3 2 1 ]),...
            {'w',...
            '$\sigma \dot{B} \  $   Self.Sim.',...
            '$\sigma \dot{B} \ $   8000 EOFs',...
            '$\sigma \dot{B} \ $   2000 EOFs',...
            '$\sigma \dot{B} \ $   200 EOFs',...
            '$\sigma \dot{B} \ $   20 EOFs',...
            '$\sigma \dot{B} \ $   2 EOFs'},...
            'Interpreter','latex',...
            'Location','northwestoutside');
        %         'Location','northeastoutside');
        % %         'Location','northeastoutside');
        set(gca, 'Position', axP);
        set(gcf,'children',flipud(get(gcf,'children')));
        v=sigma_dBt_on_sq_dt;
        2*mean(v(:).^2)
        drawnow
        
        folder_simu = model.folder.folder_simu;
        eval( ['print -depsc ' folder_simu ...
            '/Comp_ADSD_SelfSim_EOF.eps']);
        keyboard;
    end
end

% Other parameters
taille_police = 12;
id_part=1;
type_data = model.type_data;
% if plot_random_IC
%     if  random_IC_large
%         folder_simu = [ model.folder.folder_simu ...
%             '/large_IC_perturb' ];
%     else
%         folder_simu = [ model.folder.folder_simu ...
%             '/small_IC_perturb' ];
%     end
% else
%     folder_simu = [ model.folder.folder_simu ...
%         '/no_IC_perturb' ];
% end
folder_simu = [ model.folder.folder_simu ...
    '/comp_EOF_SelfSim' ];

% folder_simu = model.folder.folder_simu;
plot_moments = model.advection.plot_moments;
%plot_epsilon_k = model.advection.plot_epsilon_k;
map = model.folder.colormap;

fft_b_classic = fft_b_SelfSim;

T_adv_part_HR = real(ifft2(fft_buoy_part_ref));
% if  plot_random_IC
T_adv_part_classic = real(ifft2(fft_b_classic));
% else
%     T_adv_part_classic = nan(size(fft_b));
% end
T_adv_part_EOF = real(ifft2(fft_b));
T_adv_part_SelfSim = real(ifft2(fft_b_SelfSim));
mean_T_EOF = mean(T_adv_part_EOF,4);
mean_T_SelfSim = mean(T_adv_part_SelfSim,4);
mean_T_classic = mean(T_adv_part_classic,4);

% if eval(day) == 15
%     keyboard;
% end


%% Spatial Error
% odgT = max(abs(T_adv_part_HR(:)));
odgT = sqrt(mean(T_adv_part_HR(:).^2));

% if cropped
%     T_adv_part_HR = fft_buoy_part_ref;
%     T_adv_part_classic = fft_b_classic;
%     T_adv_part = fft_b;
%     mean_T = mean(T_adv_part,4);
% end

error_classic = abs(T_adv_part_HR - T_adv_part_classic);
error_sto_EOF = abs(bsxfun(@plus, T_adv_part_HR, - T_adv_part_EOF ));
error_mean_EOF =  abs(T_adv_part_HR - mean_T_EOF);
error_sto_SelfSim = abs(bsxfun(@plus, T_adv_part_HR, - T_adv_part_SelfSim ));
error_mean_SelfSim =  abs(T_adv_part_HR - mean_T_SelfSim);
error_mean_classic =  abs(T_adv_part_HR - mean_T_classic);
%         error_classic = T_adv_part_HR - T_adv_part_classic;
%         error_sto_EOF = bsxfun(@plus, T_adv_part_HR, - T_adv_part );
%         error_mean_EOF =  T_adv_part_HR - mean_T;
std_T_EOF = std(T_adv_part_EOF,0,4);
estim_error_EOF = 1.96 * std_T_EOF;
std_T_SelfSim = std(T_adv_part_SelfSim,0,4);
estim_error_SelfSim = 1.96 * std_T_SelfSim;
std_T_classic = std(T_adv_part_classic,0,4);
estim_error_classic = 1.96 * std_T_classic;


scale = 1.5*max(abs(error_classic(:))/odgT);
% odgT = 1;

ax = ( [min(abs([error_classic(:); estim_error(:) ])/odgT) ...
    max(abs([error_classic(:); estim_error(:) ])/odgT)]);
ax(2)=ax(2)/2;

%% One realisation
width = 3.7;
%             width = 3.5;
height = 3;
close(figure(17))
figure17=figure(17);
set(figure17,'Units','inches', ...
    'Position',[X0(1) X0(2) 2*width height], ...
    'PaperPositionMode','auto');
%             figure(17);
% colormap jet
subplot(1,2,1)
imagesc(model.grid.x,model.grid.y,(abs(error_classic(:,:,:,id_part)')/odgT));
axis xy;axis equal;caxis(ax);colorbar
set(gca,...
    'Units','normalized',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',taille_police,...
    'FontName','Times')
ylabel('y(m)',...
    'FontUnits','points',...
    'interpreter','latex',...
    'FontSize',taille_police,...
    'FontName','Times')
xlabel('x(m)',...
    'interpreter','latex',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',taille_police,...
    'FontName','Times')
title('\hspace{1.5cm} One real. Self.Sim.',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'interpreter','latex',...
    'FontSize',12,...
    'FontName','Times')
%             title('error of classical SQG')
drawnow
subplot(1,2,2)
imagesc(model.grid.x,model.grid.y,(abs(error_sto_EOF(:,:,:,id_part)')/odgT));
axis xy;axis equal;
%             ax=caxis;
caxis(ax);
colorbar
set(gca,...
    'Units','normalized',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',taille_police,...
    'FontName','Times')
ylabel('y(m)',...
    'FontUnits','points',...
    'interpreter','latex',...
    'FontSize',taille_police,...
    'FontName','Times')
xlabel('x(m)',...
    'interpreter','latex',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',taille_police,...
    'FontName','Times')
title('\hspace{1.5cm} One real. EOF',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'interpreter','latex',...
    'FontSize',12,...
    'FontName','Times')
%             title('error one realization')
drawnow;
eval( ['print -depsc ' folder_simu '/spatial_error/' day '.eps']);

%% Bias
close(figure(21))
figure21=figure(21);
%             colormap hot
% colormap jet
set(figure21,'Units','inches', ...
    'Position',[X0(1) X0(2) 2*width height], ...
    'PaperPositionMode','auto');


subplot(1,2,1)
imagesc(model.grid.x,model.grid.y,(abs(error_mean_classic')/odgT));
axis xy;axis equal;caxis(ax);colorbar
set(gca,...
    'Units','normalized',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',taille_police,...
    'FontName','Times')
ylabel('y(m)',...
    'FontUnits','points',...
    'interpreter','latex',...
    'FontSize',taille_police,...
    'FontName','Times')
xlabel('x(m)',...
    'interpreter','latex',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',taille_police,...
    'FontName','Times')
title('\hspace{1.5cm} Bias Self.Sim.',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'interpreter','latex',...
    'FontSize',12,...
    'FontName','Times')
%             title('error of classical SQG')
drawnow

subplot(1,2,2)
imagesc(model.grid.x,model.grid.y,(abs(error_mean_EOF')/odgT));
axis xy;axis equal;caxis(ax);colorbar
%             title('error mean')
set(gca,...
    'Units','normalized',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',taille_police,...
    'FontName','Times')
ylabel('y(m)',...
    'FontUnits','points',...
    'interpreter','latex',...
    'FontSize',taille_police,...
    'FontName','Times')
xlabel('x(m)',...
    'interpreter','latex',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',taille_police,...
    'FontName','Times')
title('\hspace{1.5cm} Bias EOF',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'interpreter','latex',...
    'FontSize',12,...
    'FontName','Times')


eval( ['print -depsc ' folder_simu '/spatial_bias/' day '.eps']);



%% Estimation of the bias
close(figure(19))
figure19=figure(19);
%             colormap hot
% colormap jet
set(figure19,'Units','inches', ...
    'Position',[X0(1) X0(2) 3*width height], ...
    'PaperPositionMode','auto');



subplot(1,3,1)
imagesc(model.grid.x,model.grid.y,(abs(error_mean_EOF')/odgT));
axis xy;axis equal;caxis(ax);colorbar
%             title('error mean')
set(gca,...
    'Units','normalized',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',taille_police,...
    'FontName','Times')
ylabel('y(m)',...
    'FontUnits','points',...
    'interpreter','latex',...
    'FontSize',taille_police,...
    'FontName','Times')
xlabel('x(m)',...
    'interpreter','latex',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',taille_police,...
    'FontName','Times')
title('\hspace{1.5cm} Bias EOF',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'interpreter','latex',...
    'FontSize',12,...
    'FontName','Times')
subplot(1,3,2)
imagesc(model.grid.x,model.grid.y,(abs(estim_error')/odgT));
axis xy;axis equal;caxis(ax);colorbar
%             title('error estimated')
set(gca,...
    'Units','normalized',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',taille_police,...
    'FontName','Times')
ylabel('y(m)',...
    'FontUnits','points',...
    'interpreter','latex',...
    'FontSize',taille_police,...
    'FontName','Times')
xlabel('x(m)',...
    'interpreter','latex',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',taille_police,...
    'FontName','Times')
title('\hspace{1.5cm} Estim. error EOF',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'interpreter','latex',...
    'FontSize',12,...
    'FontName','Times')
drawnow;


subplot(1,3,3)
imagesc(model.grid.x,model.grid.y,(abs(estim_error_classic')/odgT));
axis xy;axis equal;caxis(ax);colorbar
set(gca,...
    'Units','normalized',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',taille_police,...
    'FontName','Times')
ylabel('y(m)',...
    'FontUnits','points',...
    'interpreter','latex',...
    'FontSize',taille_police,...
    'FontName','Times')
xlabel('x(m)',...
    'interpreter','latex',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',taille_police,...
    'FontName','Times')
title('\hspace{1.5cm} Estim. error Self.Sim.',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'interpreter','latex',...
    'FontSize',12,...
    'FontName','Times')
%             title('error of classical SQG')
drawnow

eval( ['print -depsc ' folder_simu '/Estim_spatial_bias/' day '.eps']);


%% Spectral error
LineWidth = 1.3;
MarkerSize = 8;
Color1=[0.8 0.1 0.1];
%             Color1=[0.8 0.0 0.1];
Color2=[0.1 0.0 0.8];

close(figure(18))
figure18=figure(18);
[~,name_plot] = fct_spectrum53_HR_error( ...
    model_randomIC,fft2(error_mean_classic)/odgT,'m+-');
% set(name_plot,'LineWidth',LineWidth,...
%     'MarkerSize',MarkerSize,...
%     'Color',Color2);
hold on;
% [~,name_plot] = fct_spectrum53_HR_error( ...
%     model,fft2(error_sto_EOF(:,:,:,1))/odgT,'--b');
% set(name_plot,'LineWidth',LineWidth,...
%     'MarkerSize',MarkerSize,...
%     'Color',Color1);
[~,name_plot] = fct_spectrum53_HR_error( ...
    model,fft2(error_mean_EOF)/odgT,'b+-');
set(name_plot,'LineWidth',LineWidth,...
    'MarkerSize',MarkerSize,...
    'Color',Color1);
[~,name_plot] = fct_spectrum53_HR_error( ...
    model_randomIC,fft2(estim_error_classic)/odgT,'mo-');
% set(name_plot,'LineWidth',LineWidth,...
%     'MarkerSize',MarkerSize,...
%     'Color',Color2);
[~,name_plot] = fct_spectrum53_HR_error( ...
    model,fft2(estim_error)/odgT,'bo-');
set(name_plot,'LineWidth',LineWidth,...
    'MarkerSize',MarkerSize,...
    'Color',Color1);

%             Lx = model.grid.MX(1)*model.grid.dX(1);
%             dx_LR = Lx/128;
%             kappa_LR = pi / dx_LR;
%             ax = axis;
%             plot(kappa_LR *[1 1],ax(3:4),'g')
hold off;
set(gca,'XGrid','on','XTickMode','manual');
width = 12;
%width = 9;
height = 3;
set(figure18,'Units','inches', ...
    'Position',[X0(1) X0(2) width height], ...
    'PaperPositionMode','auto');
set(gca,'YGrid','on')


set(gca,...
    'Units','normalized',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',taille_police,...
    'FontName','Times')
% ylabel('$|\hat{e}(\kappa)|^2$',...
ylabel('$E(\kappa) \  (rad^{-1}.m)$',...
    'FontUnits','points',...
    'interpreter','latex',...
    'FontSize',taille_police,...
    'FontName','Times')
xlabel('$\kappa (rad.m^{-1})$',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',taille_police,...
    'interpreter','latex',...
    'FontName','Times')
title('Spectrum of the errors and its estimation',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'interpreter','latex',...
    'FontSize',12,...
    'FontName','Times')


legend('Bias Self.Sim.','Bias EOF',...
    'Estim. error Self.Sim.','Estim. error EOF',...
    'Location','westoutside');

drawnow;
eval( ['print -depsc ' folder_simu '/Estim_spectral_error/' day '.eps']);

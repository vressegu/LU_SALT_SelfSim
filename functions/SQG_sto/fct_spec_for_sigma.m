function [spectrum_w_a, kappa, d_kappa, k] = ...
    fct_spec_for_sigma(model, grid_, ft_w, N_ech)
% Compute KESD or ADSD  of ft_w (for the definition of sigma)
% 
%% Outputs :
% - spectrum_w_a : KESD or ADSD  of ft_w
%         size : Mx x My x 2 x N_patch          
% - kappa : wave-number on 1D axis
%         size : 1 x min(Mx,My)/2       
% - d_kappa : wave-number step
%         size : 1       
% - k : wave-number on 2D grid
%         size : Mx x My    
%
%% Inputs :
% - model : structure which gathers all simulation information
% - grid_ : structure which gathers grid information
% - ft_w : fft of the (large-scale or fully-resolved) velocity
%         size : Mx x My x 2 x N_patch
% - N_ech : number of realizations (or number of patches)
%         size : 1
%

% Squared norm
ft_w2=abs(ft_w).^2;
ft_w2=sum(ft_w2,3);
% for the (1/2) in the KE definition
ft_w2 = 0.5 * ft_w2;

% Get parameters
MX=grid_.MX;
dX=grid_.dX;
[s1 s2 s3 s4] = size(ft_w2);
if any([s1 s2 s3 s4]~=[ MX 1 N_ech] )
    error('wrong size');
end
clear s1 s2 s3 s4
if any( mod(MX,2)~=0)
    error('the number of grid points by axis need to be even');
end
PX=MX/2;
ft_w2(PX(1)+1,:,:,:)=0;
ft_w2(:,PX(2)+1,:,:)=0;

%% Wave vector
k = grid_.k.k;
k(PX(1)+1,:)=inf;
k(:,PX(2)+1)=inf;
k=k(:);

%% Wave number
MX_kappa=grid_.MX;
M_kappa=min(MX_kappa);
P_kappa= M_kappa/2;
d_kappa = 2*pi/sqrt(prod(MX_kappa.* grid_.dX));
kappa= d_kappa * ( 0:(P_kappa-1) ) ;

%% Masks associated with the rings of iso wave number
kappa = kappa';

%% Spectrum
% Integration over the rings of iso wave number
ft_w2 = reshape(ft_w2 , [prod(grid_.MX) N_ech ]);
spectrum_w = grid_.k.idx_k' * ft_w2; %  [ P_kappa N_ech ]

% Division by prod(grid_.MX) because of the Parseval theorem for
% discrete Fourier transform
% Division by prod(grid_.MX) again in order to the integration
% of the spectrum over the wave number yields the energy of the
% velocity averaged (not just integrated) over the space
spectrum_w = 1/prod(grid_.MX)^2 * spectrum_w;

% Division by the wave number step
d_kappa = kappa(2)-kappa(1);
spectrum_w = spectrum_w / d_kappa;

% From velocity spectrum to absolute diffusity by scale
spectrum_w_a = zeros(size(spectrum_w));
spectrum_w_a(2:end,:) = bsxfun( @times, kappa(2:end).^(-3/2) , ...
    spectrum_w(2:end,:).^(1/2) );

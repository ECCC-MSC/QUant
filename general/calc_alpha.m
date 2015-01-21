function [c, k, x, alpha_w, alpha_s, alpha_s_scat, alpha_s_visc] = calc_alpha(T, S, z, freq, as, isquartz, M, inDB)
%calculate the attenuation due to scattering and due to propagation in the
%medium in units of dB/m
% M is concentration in kg/m3 or g/L
% as is particle radius in m
% freq is frequency in Hz

rho = 998; %g/cm^3 water density

c = (1449.2 + 4.6*T - 0.055*T.^2 + .00029*T.^3 + (1.34 - 0.01*T)*(S-35)+0.016*z); % SPEED OF SOUND IN m/s
% from page 3 clay and medwin
k = 2*pi*freq./c; 
x = k*as;

if isquartz == 1
    rhop = 2650; % density of quartz, used as the density of the sediment
    zeta = 1;
    gamma_kappa = -.93; % from Hay, 1983 p. 7526
    gamma_rho = .77; % ditto
else
    alertSteph = 'must enter scatterer properties'
    return
end


Ka = (gamma_kappa^2 + (gamma_rho^2)/3)/6; % in thorne and hanes 2002 they set this to 0.18, here I calculate 0.177, nearly the same

% component due to scattering
sigma = rhop/rho;
nu = 1.3E-6; % kinematic viscosity of water, 1.3 times 10^-6 m2/s
beta = sqrt(2*pi*freq./(2*nu));
delta = (1 + 9./(2*beta*as))/2;
s = (9./(4*beta*as)).*(1 + 1./(beta*as));
wiggle = (k.*((sigma - 1).^2)./(2*rhop)).*(s./(s.^2 + (sigma + delta).^2));

mean(wiggle);
alpha_s_scat =M.*Ka.*(x.^4) ./ ((rhop*as)./(1 + (4/3)*Ka*x.^4 + zeta*x.^2)) ; % equation 14 of Sheng and Hay 1988
alpha_s_visc = M.*wiggle;
alpha_s = alpha_s_scat + alpha_s_visc;



% from Fisher and Simmons 1977
alpha_w = (55.9  - 2.37.*T + 4.77E-2.*T.^2 - 3.48E-4.*T.^3)*(1E-15).*(freq^2); % m^-1


if inDB == 1
    alpha_w = 8.686*alpha_w; % to be in dB/m
    alpha_s = 8.686*alpha_s;
end
% How do reaction affinities influence control coefficients and chemical noise?
%
% See Article Figure 7

% -----------------------------------------------------------------------------
% It is assumed that the fluxes in the reference state come in units of mM/s
% for the present calculations, they are transformed into mol/s, so a cell volume 
% has to be given
% -----------------------------------------------------------------------------

% -----------------------------------------------------------------------------
% load variables network, network_CoHid, v, mu, result

model_dir  = [es_BASEDIR '/../resources/models-article/Hepatonet_CCM_Fluxes'];
model_name = 'Hepatonet_Aerobic_ATP_Regeneration'; 
ref_state  = 'Hepatonet_Aerobic_ATP_Regeneration_reference_state'; 

cd(model_dir); load(model_name); load(ref_state);

% For clarity, omit NAD compounds from graphics 
network_CoHid.graphics_par.omitmetabolites       = {'NAD+ [c]', 'NAD+ [m]', 'NADH [c]', 'NADH [m]'};

% -----------------------------------------------------------------------------
% Check for stable steady state

eigmax = max(real(eig(reference_state.control.M)));
if eigmax > 0, error('Unstable steady state'); end
network.kinetics = reference_state.kinetics;
volume           = 10000 * 10^-18;      % 10 000  um^3 -> transform into [m^3]
v_mM_per_s       = reference_state.v;            % 1 = 1 mM/min; transform [mM/min] into mM/s
v                = v_mM_per_s * volume; % transform [mM/s] into mol/s
c                = reference_state.c;
mu               = reference_state.mu;

basename   = 'energetics_control_and_noise';

% -----------------------------------------------------------------------------
% Bounds for omega: 
%  slow timescale 150 min (long experiment) -> 10000 s
%  fast timescale reaction events  (enzyme 1000 copies/cell, kcat = 10) 1/10000 s
%  compute broader omega range (in order to make computed integrals over frequencies reliabale)

omega_fast = 2 * pi * 1;         % 1/sec
omega_slow = 2 * pi * 1/(1000); % 1/(1000) 1/sec (16.7 minutes period)
omega_list = 2 * pi * 10.^[-4:0.1:4];
tau_list   = 10.^[-4:0.1:4];

[Sigma_c, Sigma_j, Sigma_c_list, Sigma_j_list, Sigma_c_specific_1Hertz, Sigma_j_specific_1Hertz] = energetics_control_and_noise(network, network_CoHid, c, v, mu, omega_list, tau_list, omega_fast, omega_slow, [], basename, es_options, es_constraints, volume);

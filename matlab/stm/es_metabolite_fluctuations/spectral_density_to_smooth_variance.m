function V = spectral_density_to_smooth_variance(omega_list,S_list,tau_list)

% V = spectral_density_to_smooth_variance(omega_list,S_list,tau_list)
%
% Variances of fluctuations after smoothing with Gaussian kernels of widths tau
% (computed from spectral densities)
%
% Format of S_list: 
%  Rows: elements
%  Columns: frequencies (given in omega_list, increasing)

if size(S_list,1) > 1,
  for it = 1:size(S_list,1),
    V(it,:) =   spectral_density_to_smooth_variance(omega_list,S_list(it,:),tau_list);
  end
  return
end
  
clear V

% the values in omega_list must be increasing!
if sum(omega_list(2:end)-omega_list(1:end-1)<0), error('The values in omega_list must be increasing'); end

% plot(omega_list,S_list); set(gca,'XScale','log','YScale','log');

slope_1  = [S_list(2)-S_list(1)]/[omega_list(2)-omega_list(1)];
S_list_0 = S_list(1) - omega_list(1) * slope_1;

omega_list_aug = [0; column(omega_list)];
S_list_aug     = [S_list_0; column(S_list)];

for it = 1:length(tau_list),

tau = tau_list(it);

smoothing_function = exp(-omega_list_aug.^2 * tau^2 / 2);

S_smoothed = S_list_aug .* smoothing_function.^2;

% values in interval midpoints
S_smoothed_average = 1/2 * [S_smoothed(2:end)+S_smoothed(1:end-1)];
d_omega_list_aug = omega_list_aug(2:end)-omega_list_aug(1:end-1);

V(1,it) = sum(S_smoothed_average.*d_omega_list_aug);

end


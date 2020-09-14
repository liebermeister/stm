function result = es_sample_elasticities(N, W, ind_ext, es_constraints, es_options, c0, c, u, J, Keq, mu, K, Kma, zeta, A)

% ES_SAMPLE_ELASTICITIES - Elasticity phase
%
% result = es_sample_elasticities(N, W, ind_ext, es_constraints, es_options, c0, c, u, J, Keq, mu, K, Kma, zeta, A)
%
% Elasticity sampling: run elasticity sampling phase and MCA phase 
%
%
% Inputs (= Outputs of es_sample_steady_state; with nm: # metabolites; nr: # reactions)
%
%   N         - Stoichiometric matrix (nm x nr)
%   W         - Allosteric regulation matrix (nr x nm)
%   ind_ext   - indices of external metabolites
%   c0        - Vector of chemical equilibrium concentrations (optional)
%   c         - Vector of concentrations
%   u         - Vector of enzyme levels
%   J         - Flux vector
%   Keq       - Vector of equilibrium constants
%   mu        - Vector of chemical potentials
%   K         - Kernel matrix 
%   Kma       - Vector of mass-action ratios
%   zeta      - zeta vector computed from reaction affinities
%   A         - Vector of reaction affinities (thermodynamic driving forces)
%
% For the inputs 'es_constraints' and 'es_options', see 'es_default_options'
%
% Fields in 'es_constraints' and 'es_options' relevant to this function:
%
%   es_constraints.beta_M_fix     - predefined beta values (only used if non-zero)
%   es_constraints.beta_A_fix     - predefined beta values (only used if non-zero)
%   es_constraints.beta_I_fix     - predefined beta values (only used if non-zero)
%   es_options.kinetic_law        - Rate law {'cs','ms', ...}
%   es_options.h                  - vector of Hill exponents (one for each reaction)
%   es_options.no_equilibrium     - flag: assume no equilibrium for inactive reactions
%   es_options.KV_prior_mean      - .. and KV value to be used in this case
%   es_options.flag_test          - Run tests? (Boolean)
%   es_options.flag_second_order  - Compute second-order elasticities? (Boolean)
%   es_options.zv                 - Flux weights for 2nd order resp. coeff output function
%   es_options.zc                 - Conc. weights for 2nd order resp. coeff output function
%
% All output data are stored in a structure 'result' with fields:
%   result.c            - Metabolite levels      
%   result.c0           - Equilibrium metabolite levels      
%   result.u            - Enzyme levels
%   result.A            - Reaction affinities
%   result.zeta         - zeta vector (computed from reaction affinities)
%   result.mu           - Chemical potentials
%   result.J            - Fluxes
%   result.v            - Reaction rates (computed from rate laws, v=J)
%   result.v_plus       - Microscopic forward rates
%   result.v_minus      - Microscopic reverse rates
%   result.Keq          - Equilibrium constants
%   result.Kma          - Mass-action ratios
%   result.KA           - Activation constants
%   result.KI           - Inhibition constants
%   result.KM           - Reactant constants
%   result.KV           - Velocity constants (geometric mean of catalytic constants)
%   result.Kplus        - Forward catalytic constants
%   result.Kminus       - Reverse catalytic constants
%   result.saturation   - Saturation values
%   result.elasticities - Elasticities (scaled and unscaled)
%   result.control      - Control and response coefficients (scaled and unscaled)
%   result.K            - Kernel matrix 
%   result.kinetics     - kinetic rate laws (to be inserted into 'network')
%
% See also: es_sample_model


%% Set global variables to speed up function modular_velocities
global global_structure_matrices 
global_structure_matrices = 1;
global Mplus Mminus Wplus Wminus nm nr ind_M ind_Wp ind_Wm
h = ones(nr,1);
[Mplus, Mminus, Wplus, Wminus, nm, nr, N_int,ind_M,ind_Wp,ind_Wm] = make_structure_matrices(N,W,ind_ext,h);
%% END Set global variables
 
  
% ----------------------------------------------------------------
% sample alpha, values and compute beta and gamma values

[nm,nr] = size(N);

[alpha_A, alpha_I, alpha_M] = sample_saturation_parameters(N,W,ind_ext,es_options,es_constraints);

[beta_A,gamma_A] = alpha_to_betagamma(alpha_A);
[beta_I,gamma_I] = alpha_to_betagamma(alpha_I);
[beta_M,gamma_M] = alpha_to_betagamma(alpha_M);
KA               = alpha_to_k(alpha_A,c);
KI               = alpha_to_k(alpha_I,c);
KM               = alpha_to_k(alpha_M,c);

% ----------------------------------------------------------------
% compute remaining model parameters (KV, Kplus, Kminus)

% compute KV (given sampled u and comparison between J and preliminary v

KM_full     = full(KM); KM_full(KM_full==0) = 1;
this_KV     = ones(nr,1);
v_pre       = modular_velocities(es_options.kinetic_law,N,W,ind_ext,u,c,KA,KI,KM,this_KV,Keq,es_options.h);
ind_off     = find([v_pre==0].*[A~=0]);
ind_off_A   = find([v_pre==0].*[A==0]);

%-----------------------------
% FIX, TO BE REMOVED LATER:
es_options = join_struct(struct('no_equilibrium',1),es_options);
%-----------------------------


if es_options.no_equilibrium,
  if length(ind_off),  display('Vanishing rates: setting enzyme values = 0'); end
  u(ind_off)     = 0;
  v_pre(ind_off) = 1;
  KV          = J./v_pre;
  KV(ind_off) = es_options.KV_prior_mean;
else,
  if length(ind_off),  display('Vanishing rates: setting KV values = 0'); end
  v_pre(ind_off) = 1;
  KV          = J./v_pre;
  KV(ind_off) = 0;
end

% arbitrary choice (since nothing is known...)

u(ind_off_A)  = nanmedian(u);
KV(ind_off_A) = nanmedian(KV);
KV(isnan(KV)) = nanmedian(KV);

if sum(KV<0), error('Problem with flux directions'); end

% compute reaction velocities
[v, v_plus, v_minus] = modular_velocities(es_options.kinetic_law,N,W,ind_ext,u,c,KA,KI,KM,KV,Keq,es_options.h);
[Kplus,Kminus]       = ms_compute_Kcat(N,KM,KV,Keq);

% ----------------------------------------------------------------
% Tests

if es_options.flag_test,
  [Mplus, Mminus, Wplus, Wminus, nm, nr, N_int] = make_structure_matrices(N,W,ind_ext,es_options.h);
  dc_td             = N_int * J 
  if length(c0),
    def_keq_mismatch  = [log(Keq) - N' * log(c0)]
    v_eq = modular_velocities(es_options.kinetic_law,N,W,ind_ext,u,c0,KA,KI,KM,KV,Keq,es_options.h)
    max(abs( def_keq_mismatch)) 
    max(abs( v_eq  )) 
  end
  sign_J_A_mismatch = [sign(J) - sign(A)];
  max(abs( dc_td)) 
  max(abs( sign_J_A_mismatch )) 
end


% ----------------------------------------------------------------
% Compute elasticities:
%  o scaled elasticities (single reaction directions) only require alpha values
%  o unscaled and total elasticities also require velocities and concentrations

v_plus_fallback  = ones(size(v));
v_minus_fallback = ones(size(v));

E = compute_modular_elasticities(es_options.kinetic_law, N, W, ind_ext, alpha_A, alpha_I, alpha_M, v, A, u, c, es_options.h, v_plus_fallback, v_minus_fallback, es_options.flag_second_order);

if find(~isfinite(E.un_E_c)), 
  warning('Elasticity matrix contains non-finite values; I replace them by zeros'); 
  E.un_E_c(isfinite(E.un_E_c)==0) = 0;
  E.un_E_u(isfinite(E.un_E_u)==0) = 0;
end


% -------------------------------------------------------------
% Compute unscaled and scaled response coefficients 

Ec   = E.un_E_c;
Eu   = E.un_E_u;
Es   = E.un_E_c(:,ind_ext);
Ep   = [Eu Es];

external  = zeros(nm,1); external(ind_ext) = 1;

if prod(size(N)) > 10^4, % compute approximation of NR_int and L_int 
  epsilon_nullspace = 10^-20;
  N_int = N(find(external==0),:);
  %%[U,S,V]  = svds(sparse(N_int),100);
  [U,S,V]  = svd(full(N_int));
  ii =  find(abs(diag(S))>epsilon_nullspace); 
  L_int  = U(:,ii) * S(ii,ii);
  NR_int = V(:,ii)';
  [CJ, CS, L_int, NR_int, M, M_adj] = control_coefficients(N, Ec, external, [], NR_int, L_int);
else,
  [CJ, CS, L_int, NR_int, M, M_adj] = control_coefficients(N, Ec, external);
end


% -------------------------------------------------------------
% second order

if es_options.flag_second_order * [length(es_options.zc) + length(es_options.zv) == 0],

   %% compute full second-order response matrices
   
   Ecc  = E.un_E_cc;
   Ecu  = E.un_E_cu;
   Euu  = E.un_E_uu;
   Ecs  = E.un_E_cc(:,:,ind_ext);
   Ess  = E.un_E_cc(:,ind_ext,ind_ext);
   Esu  = E.un_E_cu(:,ind_ext,:);
  
   np  = size(Ep,2);
   Ecp = sptensor([nr,nm,np]);
   Epp = sptensor([nr,np,np]);
   Ecp(1:nr,1:nm,1:nr)       = Ecu;
   Ecp(1:nr,1:nm,nr+1:np)    = Ecs;
   Epp(1:nr,1:nr,1:nr)       = Euu;
   Epp(1:nr,1:nr,nr+1:np)    = permute(Esu,[1 3 2]);
   Epp(1:nr,nr+1:np,1:nr)    = Esu;
   Epp(1:nr,nr+1:np,nr+1:np) = Ess;
   
   [RSp,RJp,RSpp,RJpp] = response_coefficients_sparse(CS,Ec,Ep,Ecc,Ecp,Epp);

   %% "response" effect of external metabolite to itself
   RSp(ind_ext,nr+1:np) = eye(length(ind_ext));
   control              = es_control_analysis(RSp, RJp, RSpp, RJpp, c, v, u, ind_ext);

else
  
  [RSp_un,RJp_un] = response_coefficients(CS,Ec,Ep);
  %% "response" effect of external metabolite to itself
  np  = size(Ep,2);
  RSp_un(ind_ext,nr+1:np) = eye(length(ind_ext));
  [RSp_sc, RJp_sc] = norm_response_coefficients(c, v, [u; c(ind_ext)], RSp_un, RJp_un);

  control.CS  = CS;
  control.CJ  = CJ;
  control.RSp = RSp_un;
  control.RJp = RJp_un;
  control.RSp_sc = RSp_sc;
  control.RJp_sc = RJp_sc;
  control.CS_sc  = diag(1./(c)) * control.CS  * diag(v);
  control.CJ_sc  = diag(1./(v)) * control.CJ  * diag(v);
  n_ext = length(ind_ext);
  control.RSu_un = control.RSp(:,1:nr);
  control.RJu_un = control.RJp(:,1:nr);
  control.RSu_sc = diag(1./(c)) * control.RSp(:,1:nr) * diag(u);
  control.RJu_sc = diag(1./(v)) * control.RJp(:,1:nr) * diag(u);
  control.RSs_sc = diag(1./(c)) * control.RSp(:,end-n_ext+1:end) * diag(c(ind_ext));
  control.RJs_sc = diag(1./(v)) * control.RJp(:,end-n_ext+1:end) * diag(c(ind_ext));
  control.RSs_un = control.RSp(:,end-n_ext+1:end);
  control.RJs_un = control.RJp(:,end-n_ext+1:end);

  %% Second order effects on specified target function 
  %% (some linear function of concentrations and fluxes)
  if [length(es_options.zc) + length(es_options.zv)],    
    [control.Rtarget_sc_u, control.Rtarget_sc_uu] = compute_modular_response_second(...
        es_options.zc, es_options.zv, es_options.kinetic_law, N, W, ...
        ind_ext, alpha_A, alpha_I, alpha_M, v, A, u, c, es_options.h, ...
        v_plus_fallback, v_minus_fallback, CS, CJ);
  end
  
end
  
%% fix response coefficient of external concentration w.r.t. itself
for ittt = 1:length(ind_ext),
  control.RSs_sc(ind_ext(ittt),ittt) = 1; 
end

if length(es_options.zc),
  [control.Rtarget_sc_u, control.Rtarget_sc_uu] = compute_modular_response_second(es_options.zc,es_options.zv,es_options.kinetic_law,N, W, ind_ext, alpha_A, alpha_I, alpha_M, v, A, u, c, es_options.h, v_plus_fallback, v_minus_fallback,CS,CJ);
end

control.M     = M;
control.M_adj = M_adj;

eigmax = max(real(eig(M)));

control.stable = [eigmax <= 0];

if control.stable,
  display('The metabolic state is stable.');
else
  display('WARNING: The metabolic state is unstable.');
end

% -------------------------------------------------------------
% write all results to data structure 'result'

result.c                  = c;
result.c0                 = c0;
result.u                  = u;
result.A                  = A;
result.zeta               = zeta;
result.mu                 = mu;
result.J                  = J; % for double checking between J and v
result.v                  = v;
result.v_plus             = v_plus;
result.v_minus            = v_minus;
result.Keq                = Keq;
result.Kma                = Kma;
result.KA                 = KA;
result.KI                 = KI;
result.KM                 = KM;
result.KV                 = KV;
result.Kplus              = Kplus;
result.Kminus             = Kminus;
result.saturation.alpha_A = alpha_A;
result.saturation.alpha_I = alpha_I;
result.saturation.alpha_M = alpha_M;
result.saturation.beta_A  = beta_A;
result.saturation.beta_I  = beta_I;
result.saturation.beta_M  = beta_M;
result.saturation.gamma_A = gamma_A;
result.saturation.gamma_I = gamma_I;
result.saturation.gamma_M = gamma_M;
result.elasticities       = E;
result.control            = control;
result.control.CS         = CS;
result.control.CJ         = CJ;
result.K                  = K;


% -------------------------------------------------------------
% kinetic parameters ready to use

result.kinetics.type = es_options.kinetic_law;     
result.kinetics.u    = result.u;     
result.kinetics.c    = result.c;     
result.kinetics.KA   = result.KA;    
result.kinetics.KI   = result.KI;    
result.kinetics.KM   = result.KM;    
result.kinetics.KV   = result.KV;    
result.kinetics.Keq  = result.Keq;
result.kinetics.h    = es_options.h;


% ----------------------------------------------------------------
% test: do the first elasticities coincide with numerically computed elasticities?

if 0,
  network          = network_construct(N,ones(size(v)),ind_ext);
  network.kinetics = result.kinetics;
  Ec               = elasticities(network,result.c);
  
  if [norm(full(Ec)) - norm(full(E.un_E_c))] / norm(full(Ec)) > 10^-3,
    warning('Error in elasticities');
    [norm(full(Ec)) - norm(full(E.un_E_c))]
  end
end

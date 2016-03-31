function [c,dmu] = es_fix_constraint_violations(v,c,dmu,constraints,N)

%ES_FIX_CONSTRAINT_VIOLATIONS Change metabolic state variables to enforce a feasible state
%
% [c,dmu] = es_fix_constraint_violations(v,c,dmu,constraints,N)
%
% 1. Threshold small and large concentrations
% 2. Recalculate reaction affinities accordingly
% 3. Threshold the resulting reaction affinities, giving up their proper 
%    relationship to concentrations

eval(default('constraints','[]'));

if isempty(constraints), constraints = struct; end 

constraints_default.A_limit_min  = 2;
constraints_default.A_limit_max  = 200;
constraints_default.c_limit_min  = 10^-6;
constraints_default.c_limit_max  = 100;

constraints = join_struct(constraints_default,constraints);

% fix very small concentrations
c_new            = c;
ind_small        = find([c<constraints.c_limit_min]);
c_new(ind_small) = constraints.c_limit_min;

% fix very large concentrations
ind_large        = find([c>constraints.c_limit_max]);
c_new(ind_large) = constraints.c_limit_max;

dmu = dmu + RT * N' * [log(c_new) - log(c)]; 

ind_violated = find( [[v~=0] .* [sign(v.*dmu)~=-1]]  + isnan(dmu));

if length(ind_violated),
  warning(sprintf('%d sign constraints violated; changing some of the reaction affinities',length(ind_violated)));
end

% Fix wrong reaction affinities
dmu_geometric_median = exp(nanmean(log(abs(dmu(find(dmu))))));
dmu(ind_violated) = -sign(v(ind_violated)) .* dmu_geometric_median;
plot(v,dmu,'.');

display('Adjusting small and large reaction affinities');

% Fix very small or large reaction affinities

ind_small      = find(abs(dmu)<constraints.A_limit_min);
dmu(ind_small) = sign(dmu(ind_small)) * constraints.A_limit_min;
ind_large      = find(abs(dmu)>constraints.A_limit_max);
dmu(ind_large) = sign(dmu(ind_large)) * constraints.A_limit_max;

c = c_new;
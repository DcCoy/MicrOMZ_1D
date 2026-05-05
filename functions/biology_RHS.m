function [R,diag] = biology_RHS(z, X, T_degC, T, params, ftypes, use_MLR)
% Placeholder biology RHS (mmol m^-3 s^-1). Return zeros by default.
% Access to T(z), S(z) is provided for your rate scaling if needed.
% Replace contents with your ecosystem tendencies.

% List all possible ftypes
all_ftypes = {'aer','aoa','nob','nar','nai','nao','nir','nio','nos'};

% Temperature function %
Tfunc = params.temp.Tfunc(T_degC); 

% Boundary mask for biology
mask = true(length(z),1);
mask([1 end]) = false;

% --------------- %
% Extract tracers %
% --------------- %
ff = fields(T);
for i = 1:length(ff)
	tmp.(ff{i}) = X(:,T.(ff{i}));
end

% Set diagnostics to zero for all ftypes
for i = 1:length(all_ftypes)
	diag.([all_ftypes{i},'_L']) = zeros(size(tmp.O2));      % per-capita loss rate
	diag.([all_ftypes{i},'_Lb']) = zeros(size(tmp.O2));     % per-capita bioloss rate
	diag.([all_ftypes{i},'_Lg']) = zeros(size(tmp.O2));     % per-capita graze rate
	diag.([all_ftypes{i},'_mu']) = zeros(size(tmp.O2));     % per-capita growth rate
	diag.([all_ftypes{i},'_bio']) = zeros(size(tmp.O2));    % growth rate
	diag.([all_ftypes{i},'_loss']) = zeros(size(tmp.O2));   % loss rate
	if use_MLR == 0
		diag.([all_ftypes{i},'_limit']) = zeros(size(tmp.O2));  % limiting resource idx
		diag.([all_ftypes{i},'_mu_oxi']) = zeros(size(tmp.O2)); % growth on oxidant
		diag.([all_ftypes{i},'_mu_red']) = zeros(size(tmp.O2)); % growth on reductant
	end
end

% Cycle through active metabolisms and update diagnostics
for i = 1:length(ftypes)
	% Get ftype-specific oxidants/reductants
	this_oxi = tmp.(params.(ftypes{i}).oxi); % mmol m^-3
	this_red = tmp.(params.(ftypes{i}).red); % mmol m^-3
	% Get M-M limitation
	oxi_lim  = this_oxi ./ (this_oxi + params.(ftypes{i}).K_oxi); % 0 - 1
	red_lim  = this_red ./ (this_red + params.(ftypes{i}).K_red); % 0 - 1
	% Choose growth mode
	if use_MLR == 0
		% Get oxidant and reductant-limited growth
		oxi_up = (params.(ftypes{i}).Vmax_oxi .* oxi_lim) ./ params.(ftypes{i}).y_oxi; % s^-1 
		red_up = (params.(ftypes{i}).Vmax_red .* red_lim) ./ params.(ftypes{i}).y_red; % s^-1
		% Get growth via Liebig's law
		mu = min(oxi_up, red_up) .* mask;   % s^-1 
	else
		% Growth using multiple-resource-limitation
		% For now, use maximum growth on reductants
		mu = (params.(ftypes{i}).Vmax_red / params.(ftypes{i}).y_red) .* oxi_lim .* red_lim; 
	end
	bio      = tmp.(ftypes{i}) .* mu; % mmol C m^-3 s^-1    
	loss     = tmp.(ftypes{i}) .* (params.m_l + tmp.(ftypes{i}) .* params.m_q); % mmol C m^-3 s^-1
	Lb       = params.m_l + tmp.(ftypes{i}) .* params.m_q; % s^-1
	% Save diagnostics
	diag.([ftypes{i},'_Lb']) = Lb .* Tfunc;		    % 1 / s
	diag.([ftypes{i},'_mu']) = mu .* Tfunc;		    % 1 / s
	diag.([ftypes{i},'_bio']) = bio .* Tfunc;	    % mmol C / m3 s
	diag.([ftypes{i},'_loss']) = loss .* Tfunc;	    % mmol C / m3 s
	if use_MLR == 0
		diag.([ftypes{i},'_limit']) = oxi_up < red_up;  % logical
		diag.([ftypes{i},'_mu_oxi']) = oxi_up .* Tfunc; % 1 / s
		diag.([ftypes{i},'_mu_red']) = red_up .* Tfunc; % 1 / s
	end
end

% ------- %
% Grazing % 
% ------- %
% Assemble total prey carbon
totalPrey = 0;
for i = 1:length(ftypes)
	% Get total consumable biomass
	totalPrey = totalPrey + max(tmp.(ftypes{i}) - params.bmin,0);
end

% Calculate o2 limitation
z_o2lim = exp(-tmp.O2 ./ params.zoo.o2lim);

% Get max grazing rate (modified by temp and O2)
zoo_umax = params.zoo.u_max .* Tfunc .* (1 - z_o2lim);

% Consume prey
grazing = zeros(length(z),length(ftypes));
for i = 1:length(ftypes)
	consumablePrey = max(tmp.(ftypes{i}) - params.bmin,0);
	grazing(:,i) = zoo_umax .* tmp.zoo .* ...
		(consumablePrey ./ (totalPrey + params.zoo.K_B));
	Lg = zoo_umax .* tmp.zoo .* ...
		(1 ./ (totalPrey + params.zoo.K_B)); % s^-1
	diag.([ftypes{i},'_graze']) = grazing(:,i); % mmol C / m3 s
	diag.([ftypes{i},'_Lg']) = Lg; % s^-1
	diag.([ftypes{i},'_L']) = Lg + diag.([ftypes{i},'_Lb']); % add to bioloss rate, s^-1
end

% Get carbon assimilated into zoo growth
grazeC_total = sum(grazing,2);                       % total C grazed
grazeC_zoo   = grazeC_total .* params.zoo.G_Z;       % fraction routed to zoo
grazeC_poc   = grazeC_total .* params.graze_poc;	 % fraction routed to poc (rest to DIC) 

% Nitrogen routing
grazeN_total  = grazeC_total ./ params.CN_bio;        % total N ingested
N_need_Z      = grazeC_zoo ./ params.CN_det;          % zoo N requirement 
N_need_det    = grazeC_poc ./ params.CN_det;          % N from non-assimilated C
grazeN_to_NH4 = grazeN_total - N_need_Z - N_need_det; % Excretion due to N-rich prey
grazeN_to_NH4 = max(grazeN_to_NH4, 0);                % Guard against negatives

% Zooplankton loss and growth
loss = tmp.zoo .* (params.zoo.m_l + tmp.zoo .* params.zoo.m_q);
diag.zoo_loss  = loss .* Tfunc; % mmol C / m3 s
diag.zoo_bio   = grazeC_zoo;    % mmol C / m3 s
diag.zoo_mu    = zoo_umax;      % 1 / s
diag.zoo_graze = grazeC_total;  % mmol C / m3 s 

% Other diagnostics
diag.graze_poc = grazeC_poc;

% ------------- %
% Stoichiometry %
% ------------- %
% All carbon from mortality becomes POC
lossC_total = zeros(size(tmp.O2));
for i = 1:length(ftypes)
	lossC_total = lossC_total + diag.([ftypes{i},'_loss']); % mmol C / m3 s
end
% Save diagnostic

% Biomass N associated with that carbon 
N_from_biomass = lossC_total / params.CN_bio; % mmol PON / m3 s 

% N that should remain with newly created detritus (based on Anderson)
N_needed_for_POM = lossC_total / params.CN_det; % mmol PON / m3 s 

% Route excess N to NH4
lossN_to_NH4 = max(N_from_biomass - N_needed_for_POM,0); % mmol NH4 / m3 s 

% Other diagnostics
diag.loss_poc = lossC_total .* params.loss_poc;

% --------------------------- %
% Sources & Sinks (Chemistry) %
% --------------------------- %
% POC
dt.POC = 0 ...
		 - diag.aer_bio .* params.aer.y_red ... % aer uptake (mmol POC / m3 s)
		 - diag.nar_bio .* params.nar.y_red ... % nar uptake (mmol POC / m3 s)
		 - diag.nai_bio .* params.nai.y_red ... % nai uptake (mmol POC / m3 s)
		 - diag.nao_bio .* params.nao.y_red ... % nao uptake (mmol POC / m3 s)
		 - diag.nir_bio .* params.nir.y_red ... % nir uptake (mmol POC / m3 s)
		 - diag.nio_bio .* params.nio.y_red ... % nio uptake (mmol POC / m3 s)
		 - diag.nos_bio .* params.nos.y_red ... % nos uptake (mmol POC / m3 s)
		 + diag.loss_poc ...                    % mortality  (mmol POC / m3 s)
		 + diag.zoo_loss ...                    % zoo loss   (mmol POC / m3 s)
		 + diag.graze_poc;                      % grazing    (mmol POC / m3 s)

% O2
dt.O2 = 0 ...
		- diag.aer_bio .* params.aer.y_oxi ... % aer uptake (mmol O2 / m3 s)
	    - diag.aoa_bio .* params.aoa.y_oxi ... % aoa uptake (mmol O2 / m3 s)
		- diag.nob_bio .* params.nob.y_oxi;    % nob uptake (mmol O2 / m3 s)

% NH4
dt.NH4 = 0 ...
		 + diag.aer_bio .* params.aer.e_nh4 ... % aer excretion (mmol NH4 / m3 s) 
		 + diag.nar_bio .* params.nar.e_nh4 ... % nar excretion (mmol NH4 / m3 s) 
		 + diag.nai_bio .* params.nai.e_nh4 ... % nai excretion (mmol NH4 / m3 s) 
		 + diag.nao_bio .* params.nao.e_nh4 ... % nao excretion (mmol NH4 / m3 s) 
		 + diag.nir_bio .* params.nir.e_nh4 ... % nir excretion (mmol NH4 / m3 s) 
		 + diag.nio_bio .* params.nio.e_nh4 ... % nio excretion (mmol NH4 / m3 s) 
		 + diag.nos_bio .* params.nos.e_nh4 ... % nos excretion (mmol NH4 / m3 s) 
		 - diag.aoa_bio .* params.aoa.y_red ... % aoa uptake    (mmol NH4 / m3 s)
		 - diag.aox_bio .* params.aox.y_red ... % aox uptake    (mmol NH4 / m3 s)
	     + lossN_to_NH4 ...                     % mortality     (mmol NH4 / m3 s)
	     + diag.zoo_loss ./ params.CN_det ...   % zoo loss      (mmol NH4 / m3 s)
	     + grazeN_to_NH4;                       % grazing       (mmol NH4 / m3 s)

% NO2
dt.NO2 = 0 ...
		 + diag.aoa_bio .* params.aoa.e_no2 ... % aoa excretion (mmol NO2 / m3 s) 
		 + diag.nar_bio .* params.nar.e_no2 ... % nar excretion (mmol NO2 / m3 s) 
		 - diag.nob_bio .* params.nob.y_red ... % nob uptake    (mmol NO2 / m3 s)
		 - diag.nir_bio .* params.nir.y_oxi ... % nir uptake    (mmol NO2 / m3 s)
		 - diag.nio_bio .* params.nio.y_oxi ... % nio uptake    (mmol NO2 / m3 s)
		 - diag.aox_bio .* params.aox.y_oxi;    % aox uptake    (mmol NO2 / m3 s)

% NO3
dt.NO3 = 0 ...
         - diag.nar_bio .* params.nar.y_oxi ... % nar uptake    (mmol NO3 / m3 s)
         - diag.nai_bio .* params.nai.y_oxi ... % nai uptake    (mmol NO3 / m3 s)
         - diag.nao_bio .* params.nao.y_oxi ... % nao uptake    (mmol NO3 / m3 s)
		 + diag.nob_bio .* params.nob.e_no3 ... % nob excretion (mmol NO3 / m3 s)
		 + diag.aox_bio .* params.aox.e_no3;    % aox excretion (mmol NO3 / m3 s)

% N2O
dt.N2O = 0 ...
		 + diag.nai_bio .* params.nai.e_n2o ... % nai excretion (mmol N2O / m3 s)
		 + diag.nir_bio .* params.nir.e_n2o ... % nai excretion (mmol N2O / m3 s)
		 - diag.nos_bio .* params.nos.y_oxi;    % nos uptake    (mmol N2O / m3 s)

% N2
dt.N2  = 0 ...
		 + diag.nao_bio .* params.nao.e_n2 ...  % nao excretion (mmol N2 / m3 s)
		 + diag.nio_bio .* params.nio.e_n2 ...  % nio excretion (mmol N2 / m3 s)
		 + diag.aox_bio .* params.aox.e_n2;     % aox excretion (mmol N2 / m3 s)

% ------------------------- %
% Sources & Sinks (Biology) %
% ------------------------- %
for i = 1:length(ftypes)
	dt.(ftypes{i}) = 0 ...
		+ diag.([ftypes{i},'_bio']) ...  % biomass growth  (mmol C / m3 s)
		- diag.([ftypes{i},'_loss']) ... % biomass loss    (mmol C / m3 s)
		- diag.([ftypes{i},'_graze']);   % biomass grazing (mmol C / m3 s)
end

% zoo
dt.zoo = 0 ...
		 + diag.zoo_bio ... % biomass growth (mmol C / m3 s)
		 - diag.zoo_loss;   % biomass loss   (mmol C / m3 s)

% --------------------%
% N cycle diagnostics %
% --------------------%
diag.ammox     = diag.aoa_bio .* params.aoa.e_no2; % NO2 from ammox 
diag.nitrox    = diag.nob_bio .* params.nob.e_no3; % NO3 from nitrox
diag.denitrif1 = diag.nar_bio .* params.nar.e_no2; % NO2 from NO3-reduction
diag.denitrif2 = diag.nir_bio .* params.nir.e_n2o; % N2O from NO2-reduction
diag.denitrif3 = diag.nos_bio .* params.nos.e_n2;  % N2  from N2O-reduction
diag.denitrif4 = diag.nai_bio .* params.nai.e_n2o; % N2O from NO3-reduction
diag.denitrif5 = diag.nio_bio .* params.nio.e_n2;  % N2  from NO2-reduction
diag.denitrif6 = diag.nao_bio .* params.nao.e_n2;  % N2  from NO3-reduction
diag.anammox   = diag.aox_bio .* params.aox.e_n2;  % Anammox to N2

% ------------------- %
% Reassemble dt array %
% ------------------- %
R = zeros(size(X));
for i = 1:length(ff)
	R(:,T.(ff{i})) = dt.(ff{i});
end

% Make diagnostics single
ff = fields(diag);
for i = 1:length(ff)
	diag.(ff{i}) = single(diag.(ff{i}));
end


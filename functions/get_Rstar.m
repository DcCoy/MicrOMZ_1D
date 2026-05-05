function Rstar = get_Rstar(model)
% ---------------------------------------------------------
% ---------------------------------------------------------
% Rationale
% ---------------------------------------------------------
% dB/dt = uB - LB
% where:
% u = growth rate (1/s)
% L = loss rate (1/s)
% B = biomass (mmol/m3)
% ---------------------------------------------------------
% ---------------------------------------------------------
% At steady state:
% u = L
% where:
% u = ((Vmax*Tfunc)/Y)*(R/R+Km)
% ---------------------------------------------------------
% ---------------------------------------------------------
% Solving for R (Michaelis-Menten uptake dynamics)
% R = (L*Y*Km) / (Tfunc*Vmax - L*Y)
% ---------------------------------------------------------

% Unpack input model
unpackStruct(model);
diag = model.diag; % since 'diag' is a matlab function so this throws errors

% Temperature function
Tfunc = params.temp.Tfunc(grd.T_degC); 

% Calculate R* 
for i = 1:length(opt.ftypes);
	% Now saving per-capita loss rate (Lg + Lb)
	L = diag.([opt.ftypes{i},'_L']); 

	% Extract parameters
	unpackStruct(params.(opt.ftypes{i}));

	% Solve for reductant* and oxidant*
	Rstar.(red).(opt.ftypes{i}) = (L.*y_red.*K_red)./(Tfunc.*Vmax_red - L.*y_red);
	Rstar.(oxi).(opt.ftypes{i}) = (L.*y_oxi.*K_oxi)./(Tfunc.*Vmax_oxi - L.*y_oxi);
end


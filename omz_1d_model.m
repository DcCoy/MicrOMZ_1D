function [model] = omz_1d_model(opt)
% 1-D vertical column for ETSP-OMZ style tracers with strict BCs.
% Tracers (all mmol m^-3): O2, NO3, NO2, NH4, N2, POC
% Physics: vertical diffusion (variable Kz), optional advection per tracer
% Special: POC sinks (downward) and is forced at the top; outflows at bottom.
% Biology: hook via RHS function; enabled/disabled with a switch below.
% ------------------------------------------------------------------------

warning off
addpath functions
mkdir output

% --- POC parameters (flux boundary at top)
% Default settings
if nargin<1
	opt.F_poc_top_mean   = 20;     % mmol C m^-2 d^-1 (mean export)
	opt.F_poc_top_amp    = 0;      % relative amplitude (0.5 => ±50%)
	opt.F_poc_top_period = 30;     % days (period if varying)
	opt.ecosystem        = 'simple_no2'; % full or simple_no2
	opt.fname            = ['omz_1d_model_',opt.ecosystem]; 
end

%% ===================== USER SETTINGS ======================================
% --- domain & numerics
z_top   = 30;             % m
z_bot   = 1030;           % m
dz      = 10;             % m
ndays   = 365*20;         % model days to run

% --- switches
make_plots      = 0; % make plots
use_biology     = 1; % set true to include biology RHS
use_advection   = 1; % physics advection (POC sinking) on/off
use_nudging     = 1; % nudge to farfield profiles (O2, NO3) 
use_farfield    = 1; % farfield observations from ETSP 
use_MLR         = 0; % use multiple-resource-limitation vs. Liebig's law
use_restart     = 1; % load previous restart
save_restart    = 0; % overwrite previous restart
save_output     = 1; % Save daily-averaged output
if use_restart
	% Get restart file
	restart_file = ['output/omz_1d_model_restart_',opt.ecosystem,'.mat'];
	disp(['Loading restart file: ',restart_file]);
end

% Tracer registry
opt.tracers = {...
	'O2','NO3','NO2','NH4',...  % chemistry
	'N2O','N2','POC',...        % chemistry
	'aer','aoa','nob',...       % aerobic ftypes
	'nar','nai','nao',...       % anaerobic ftypes
	'nir','nio','nos','aox',... % anaerobic ftypes
	'zoo'};                     % grazers
opt.tnames = {...
	'O$_2$','NO$_3^-$','NO$_2^-$','NH$_4^+$',...
	'N$_2$O','N$_2$','OrgC',...
	'$B_{\mathrm{aer}}$','$B_{\mathrm{aoa}}$','$B_{\mathrm{nob}}$',...
	'$B_{\mathrm{nar}}$','$B_{\mathrm{nai}}$','$B_{\mathrm{nao}}$',...
	'$B_{\mathrm{nir}}$','$B_{\mathrm{nio}}$','$B_{\mathrm{nos{}$','$B_{\mathrm{aox}}$',...
	'Zoo'};
if strcmp(opt.ecosystem,'full');
	% do nothing
elseif strcmp(opt.ecosystem,'simple_no2')
	% remove unneeded ftypes
	ind = find(ismember(opt.tracers,{'nai','nao','nir','nos'})==1);
	opt.tracers(ind) = [];
	opt.tnames(ind)  = [];
end
% Set up registry
for i = 1:length(opt.tracers)
	T.(opt.tracers{i}) = i;
end
ntr = numel(opt.tracers);
% Save ftypes
opt.ftypes = opt.tracers(8:end-1);      % ignore zoo
opt.fnames = opt.tnames(8:end-1);       % ignore zoo
ftypes_idx = 8:(length(opt.tracers)-1); % ignore zoo
% Get distinct colors for each ftype
opt.fclrs  = colormix(length(ftypes_idx),{'w','k','y'});

% Load parameters 
get_params

% --- physics profiles
z_ref   = 150;        % m
K_top   = 5e-5;       % m^2/s
K_bot   = 1e-7;       % m^2/s
T_top   = 20;         % deg C
T_bot   = 4;          % deg C
tau_top = 30*86400;   % seconds 
tau_bot = 300*86400;  % seconds

% --- BC values (Dirichlet)
for i = 1:length(opt.tracers)
	BC.(opt.tracers{i}).top = params.bmin; % mmol m^-3
	BC.(opt.tracers{i}).bot = params.bmin; % mmol m^-3
end
% override for chemistry
%BC.O2.top  = 187.000; BC.O2.bot  = 55.000; % mmol m^-3
BC.O2.top  = 200.000; BC.O2.bot  = 55.000; % mmol m^-3
BC.NO3.top =   8.000; BC.NO3.bot = 43.000; % mmol m^-3
BC.NO2.top =   0.000; BC.NO2.bot =  0.000; % mmol m^-3
BC.NH4.top =   0.000; BC.NH4.bot =  0.000; % mmol m^-3
BC.N2O.top =   0.000; BC.N2O.bot =  0.000; % mmol m^-3
BC.N2.top  =   0.000; BC.N2.bot  =  0.000; % mmol m^-3
%BC.N2O.top =   0.013; BC.N2O.bot =  0.035; % mmol m^-3

%% ===================== SETUP ==============================================
z   = (z_top:dz:z_bot)';  nz = numel(z);

% Kv, temperature. nudging tau
Kz        = K_bot   + (K_top   - K_bot)   .* exp(-(z - z(1))/z_ref);
T_degC    = T_bot   + (T_top   - T_bot)   .* exp(-(z - z(1))/z_ref);
tau_nudge = tau_bot - (tau_bot - tau_top) .* exp(-(z - z(1))/z_ref);

% Override Kz?
if (0)
	%% --- Kz: strong surface mixing -> deep background (no bottom enhancement)
	K_bg_deep     = 5e-6;   % m^2/s  interior/deep-ocean background
	K_surf_target = 5e-4;   % m^2/s  near-surface target max
	mlz           = 20;     % m      e-fold depth scale for surface decay
	A_top = max(K_surf_target - K_bg_deep, 0);
	Kz = K_bg_deep + A_top .* exp(-(z - z_top) ./ mlz);   % monotone ↓ with depth
elseif (0)
	%% NitrOMZ-style Kz: increases with depth (m^2/s), below mixed layer
	% Magnitudes from their code
	K_top_int = 0.70 * 2.0 * 1.701e-5;   % ≈ 2.3814e-5 m^2/s (just below ML)
	K_bot_int = 1.00 * 2.0 * 1.701e-5;   % ≈ 3.4020e-5 m^2/s (deep end)

	% Shape (their flex=250 m below surface; your top is at 30 m)
	z_flex   = z_top + 250;   % = 280 m in your domain
	z_width  = 300;           % m (thickness of the transition)

	% Sigmoid ramp (positive-down depths)
	Kz = 0.5*(K_top_int + K_bot_int) ...
	   + 0.5*(K_bot_int - K_top_int) .* tanh((z - z_flex) ./ (0.5*z_width));
elseif (1)
	% --- choose magnitudes (m^2 s^-1) ---
	K_avg   = 1e-5;      % mean Kz
	K_top   = K_avg*1.5; % top Kz
	K_deep  = K_avg;     % deep Kz
	K_mid   = K_avg*2.5; % mid Kz

	% --- target the mid-depth bump ---
	zc      = 500;       % bump center depth (m, absolute)
	sigma   = 150;       % bump half-width (m) ~ controls thickness of ventilated band

	% --- 1) linear background from top->bottom ---
	z0 = z(1); z1 = z(end);
	K_bg = K_top + (K_deep - K_top) * (z - z0) / (z1 - z0);

	% --- 2) Gaussian bump (auto amplitude so it peaks at K_mid) ---
	A    = max(0, K_mid - interp1(z, K_bg, zc, 'linear', 'extrap'));   % amplitude
	bump = A * exp(-0.5*((z - zc)/sigma).^2);

	% --- total profile ---
	Kz = K_bg + bump;

	% (optional safety: keep within reasonable bounds)
	Kz = max(Kz, min(K_top, K_deep));    % never below K_min or above K_deep if you want
end

% --- far-field target profiles for nudging (consistent with Kz & BCs) ---
if use_farfield	
	fname = ['/data/project1/demccoy/iNitrOMZ/Data/farfield_ETSP_gridded.mat'];
	load(fname);
	o2    = farfield_ETSP_gridded.o2; % slighty nudge
	no3   = farfield_ETSP_gridded.no3;
	dep   = -farfield_ETSP_gridded.zgrid;
	farfield.O2  = interp1(dep,o2,z);
	farfield.NO3 = interp1(dep,no3,z); 
else
	farfield.O2  = BC.O2.bot  + (BC.O2.top  - BC.O2.bot)  .* exp(-(z - z(1))/z_ref);
	farfield.NO3 = BC.NO3.bot + (BC.NO3.top - BC.NO3.bot) .* exp(-(z - z(1))/z_ref);
end

% Tracer-specific advection (downward positive). Only POC sinks.
w_up_m_per_y  = -10;            % Upwelling speed [m yr^-1];
w_poc_m_per_d = 20;             % POC sinking speed [m d^-1]
w_poc = w_poc_m_per_d/86400;    % m s^-1 (downward positive)
w_up  = w_up_m_per_y/86400/365; % m s^-1 (downward negative)
w = zeros(nz,ntr);
if use_advection
	w(:,:) = w_up;      % upwelling for all tracers 
    w(:,T.POC) = w_poc; % POC sinking overrides upwelling 
end

% Choose dt based on CFL condition
CFL_max = 0.5;
safety = 0.9;
dt_cfl = CFL_max * dz / w_poc;
dt     = safety*dt_cfl;
nsteps_per_day = floor(86400 / dt);
nsteps_per_day = max(nsteps_per_day,1);     % at least one step/day
dt = 86400 / nsteps_per_day;

% Number of timsteps
nt  = ceil((ndays*86400)/dt);
t   = (0:nt)*dt;

% Initial conditions (simple shapes; consistent units)
X = zeros(nz,ntr);
if ~use_restart
	% Simple linear profile
	for i = 1:length(opt.tracers)
		if strcmp(opt.tracers{i},'POC');
			X(:,T.POC) = eps;
		else
			X(:,T.(opt.tracers{i}))  = lin_profile(z, z_top, z_bot, ...
				BC.(opt.tracers{i}).top, BC.(opt.tracers{i}).bot);
		end
	end
	% Set O2, NO3 to farfield
	X(:,T.O2)  = farfield.O2;
	X(:,T.NO3) = farfield.NO3;
else
	% Use final timestep of restart file
	load([restart_file]);
	% Check that restart resolution matches requested resolution
	if dz == model.grd.dz; 
		for i = 1:length(opt.tracers)
			X(:,T.(opt.tracers{i})) = model.out.(opt.tracers{i})(:,end);
		end
	else 
		% Get restart grid
		rst_depth  = model.grd.z_r;
		for i = 1:length(opt.tracers)
			% Interpolate to new grid
			rst_tracer = model.out.(opt.tracers{i})(:,end);
			X(:,T.(opt.tracers{i})) = interp1(rst_depth,rst_tracer,z);
		end
	end
	clear model
end

% Pre-build BE operators per tracer (without Dirichlet—BCs applied each step)
A_BE  = cell(ntr,1);
Mrhs  = cell(ntr,1);
w_center = zeros(nz, ntr); % no implicit advection inside A_BE
for k = 1:ntr
	% No implicit advection inside A_BE
    [A_BE{k}, Mrhs{k}] = build_BE_matrix(z, dz, Kz, zeros(nz,1), dt);
end

% Which tracers have strict Dirichlet BCs?
hasDirichletTop    = false(1,ntr);
hasDirichletBottom = false(1,ntr);
DirichletTopVal    = nan(1,ntr);
DirichletBottomVal = nan(1,ntr);

% O2 and NO3 always Dirichlet at top & bottom:
hasDirichletTop([T.O2 T.NO3 T.NO2 T.NH4 T.N2O T.N2])    = true;
hasDirichletBottom([T.O2 T.NO3 T.NO2 T.NH4 T.N2O T.N2]) = true;
DirichletTopVal([T.O2 T.NO3 T.NO2 T.NH4 T.N2O T.N2])    = ...
	[BC.O2.top, BC.NO3.top, BC.NO2.top, BC.NH4.top, BC.N2O.top, BC.N2.top];
DirichletBottomVal([T.O2 T.NO3 T.NO2 T.NH4 T.N2O T.N2]) = ...
	[BC.O2.bot, BC.NO3.bot, BC.NO2.bot, BC.NH4.bot, BC.N2O.bot, BC.N2.bot];

% POC: flux BC at top, open at bottom 
hasDirichletTop(T.POC)    = false;  DirichletTopVal(T.POC)    = NaN;  % unused
hasDirichletBottom(T.POC) = false;  DirichletBottomVal(T.POC) = NaN;  % unused

% Storage for daily-averaged tracer fields
out = struct();
for i = 1:ntr
	out.(opt.tracers{i}) = single(zeros(nz,ndays));
end
out.time  = 1:ndays;

% Save grid
grd.z_top  = z_top;
grd.z_bot  = z_bot;
grd.dz     = dz;
grd.z_r    = z;
grd.Kz     = Kz;
grd.T_degC = T_degC; 
grd.tau    = tau_nudge;

% Storage for daily-averaged rates
[~,diag] = biology_RHS(z, X, T_degC, T, params, opt.ftypes, use_MLR);   % mmol m^-3 s^-1
diag_fields = fields(diag);
if save_output
	for i = 1:length(diag_fields)
		diag_out.(diag_fields{i}) = single(zeros(nz,ndays));
	end
	% Also save advection, diffusion, and net sources-minus-sinks of tracers
	for i = 1:ntr
		diag_out.([opt.tracers{i},'_adv']) = single(zeros(nz,ndays));
		diag_out.([opt.tracers{i},'_dfz']) = single(zeros(nz,ndays));
		diag_out.([opt.tracers{i},'_sms']) = single(zeros(nz,ndays));
	end
end

% Daily accumulators
sumTracerDay = zeros(nz, ntr);
sumAdvDay    = zeros(nz, ntr);
sumDfzDay    = zeros(nz, ntr);
sumSmsDay    = zeros(nz, ntr);
sumDiagDay   = cell(length(diag_fields),1);
for i = 1:length(diag_fields)
	sumDiagDay{i} = zeros(nz,1);
end
count   = 0;
day_idx = 1;

% Display progress
prog_check = floor(linspace(1,nt,101));
prog_check = prog_check(2:end);
prog_cnt   = 1;

%% ===================== TIME LOOP ==========================================
disp(['Running ',num2str(nt),' steps (',num2str(t(end)/86400),' days -- ',num2str(t(end)/86400/365),' years)']);
for n = 1:nt
    tn = t(n+1);

    % ---- display progress
    if prog_cnt <= numel(prog_check) && n > prog_check(prog_cnt)
        disp([num2str(prog_cnt),'% complete']);
        prog_cnt = prog_cnt + 1;
    end

    % ---- Reaction tendencies (biology hook) ----
    if use_biology
        % If biology_RHS expects tracers cell array instead of T, call accordingly.
        [R,diag] = biology_RHS(z, X, T_degC, T, params, opt.ftypes, use_MLR); % mmol m^-3 s^-1
        if n == 1, ff = fields(diag); end
    else
        R    = zeros(size(X));
        diag = struct(); if n == 1, ff = {}; end
    end

    % ---- Solve each tracer (diffusion implicit; advection explicit) ----
    for k = 1:ntr
        % ---------------- Far-field nudging (O2, NO3 only) ----------------
        S_nudge = 0;
        if use_nudging
            % if k==T.O2,  
				%S_nudge = (farfield.O2  - X(:,k)) ./ tau_nudge; 
			%end
            if k==T.NO3, 
				S_nudge = (farfield.NO3 - X(:,k)) ./ tau_nudge; 
			end
        end

		% --- Top flux boundary for POC ---
		S_flux = zeros(nz,1);   % mmol m^-3 s^-1
		if k == T.POC
			% time-varying export flux (downward positive)
			tdays = n*dt/86400;
			Ftop = opt.F_poc_top_mean * (1 + opt.F_poc_top_amp * ...
				   sin(2*pi*tdays/opt.F_poc_top_period));
			Ftop = Ftop / 86400;      % convert mmol m^-2 d^-1 -> mmol m^-2 s^-1
			S_flux(1) = Ftop / dz;    % apply flux as source to top cell (mmol m^-3 s^-1)
		end

        % ---------------- Explicit advection (conservative) ----------------
        % Build face velocities (nz+1): top face = 1, bottom face = end
        wi = zeros(nz+1,1);
        if use_advection
            wc     = w(:,k);                  % cell-centred velocity (m s^-1), + = downward
            wi(2:nz) = 0.5*(wc(1:nz-1) + wc(2:nz));
            wi(1)    = wc(1);                 % top boundary face
            wi(end)  = wc(end);               % bottom boundary face
        end

        % Inflow concentrations at boundaries come from Dirichlet BCs
        Cin_top = DirichletTopVal(k);
        Cin_bot = DirichletBottomVal(k);
        if isnan(Cin_top), Cin_top = 0; end   % e.g., non-Dirichlet species
        if isnan(Cin_bot), Cin_bot = 0; end   % POC bottom is open

        % Upwind value at each face
        Cup        = zeros(nz+1,1);
        % TOP face: inflow if wi(1) > 0 (downwelling); otherwise outflow uses cell 1
        if wi(1) > 0,  Cup(1) = Cin_top; else, Cup(1) = X(1,k); end
        % INTERIOR faces
        for f = 2:nz
            if wi(f) >= 0
                Cup(f) = X(f-1,k);           % flow downward: take cell above
            else
                Cup(f) = X(f,k);             % flow upward: take cell below
            end
        end
        % BOTTOM face: inflow if wi(end) < 0 (upwelling); otherwise outflow uses cell nz
        if wi(end) < 0, Cup(end) = Cin_bot; else, Cup(end) = X(end,k); end

        % ---------------- Diagnostic vertical advection ----------------
        % Face fluxes and divergence (positive downward)
        F     = wi .* Cup;                    % mmol m^-2 s^-1
        S_adv = -(F(2:end) - F(1:end-1)) / dz; % mmol m^-3 s^-1

        % ---------------- Diagnostic vertical diffusion ----------------
        % Compute face diffusivities and gradients explicitly for diagnostics.
        Ki = zeros(nz+1,1);
        Ki(2:nz) = 0.5*(Kz(1:nz-1) + Kz(2:nz));   % interface Kz
        Ki(1)    = Kz(1);
        Ki(end)  = Kz(end);
        F_diff = -Ki .* [ (X(1,k) - DirichletTopVal(k)) / (0.5*dz); ...
                          diff(X(:,k)) / dz; ...
                          (DirichletBottomVal(k) - X(end,k)) / (0.5*dz) ];   % mmol m^-2 s^-1

        % Divergence of diffusive flux (positive when net flux converges downward)
        S_diff = -(F_diff(2:end) - F_diff(1:end-1)) / dz;   % mmol m^-3 s^-1

        % ---------------- Assemble RHS and solve (diffusion implicit) -----
        b = X(:,k) + dt * ( R(:,k) + S_adv + S_nudge + S_flux);

        % Strongly enforce Dirichlet rows on A and b
        [Amod, bmod] = apply_dirichlet_rows(A_BE{k}, b, ...
            hasDirichletTop(k),    DirichletTopVal(k), ...
            hasDirichletBottom(k), DirichletBottomVal(k));

        % Backward-Euler solve for diffusion
        X(:,k) = Amod \ bmod;

        % Belt & suspenders: enforce Dirichlet after solve
        if hasDirichletTop(k),    X(1,  k) = DirichletTopVal(k);    end
        if hasDirichletBottom(k), X(end,k) = DirichletBottomVal(k); end

        % Non-negativity (skip for N2xs if you allow negatives)
        if k ~= T.N2
            X(:,k) = max(X(:,k), 0);
        end

		% Save advection and diffusion
		diag.adv(:,k) = S_adv;  % mmol m^-3 s^-1 
		diag.dfz(:,k) = S_diff; % mmol m^-3 s^-1
    end

    % ---- Prevent extinction for biology pools (small positive floor) ----
    tmp = X(:,[ftypes_idx]);
    tmp(tmp < BC.(opt.ftypes{1}).top) = BC.(opt.ftypes{1}).top;
    X(:,[ftypes_idx]) = tmp;

	% Accumulate this step into today's running mean
	if save_output
		sumTracerDay = sumTracerDay + X;  % tracer concentration
		sumSmsDay = sumSmsDay + R;        % sources-minus-sinks
		sumAdvDay = sumAdvDay + diag.adv; % advection 
		sumDfzDay = sumDfzDay + diag.dfz; % diffusion
		for k = 1:length(diag_fields)
			sumDiagDay{k} = sumDiagDay{k} + diag.(diag_fields{k});
		end
		count = count + 1;
	end

	% Save daily-averaged tracers and diagnostics
	if save_output
		if mod(n, nsteps_per_day) == 0
			% Write daily-averaged tracers, and budget terms
			for k = 1:ntr
				out.(opt.tracers{k})(:, day_idx) = single(sumTracerDay(:,k) / max(1, count));
				diag_out.([opt.tracers{k},'_adv'])(:, day_idx) = single(sumAdvDay(:,k)./max(1,count));
				diag_out.([opt.tracers{k},'_dfz'])(:, day_idx) = single(sumDfzDay(:,k)./max(1,count));
				diag_out.([opt.tracers{k},'_sms'])(:, day_idx) = single(sumSmsDay(:,k)./max(1,count));
			end

			% Write daily-averaged diagnostics
			for f = 1:length(diag_fields)
				diag_out.(diag_fields{f})(:, day_idx) = single(sumDiagDay{f} / max(1, count));
				sumDiagDay{f}(:) = 0;
			end

			% Reset daily accumulators for physics too
			sumTracerDay(:,:) = 0;
			sumAdvDay(:,:)    = 0;
			sumDfzDay(:,:)    = 0;
			sumSmsDay(:,:)    = 0;

			% Advance day
			count = 0;
			day_idx = day_idx + 1;
		end
	% Only save final timestep, for spinup purposes
	else
		if n == nt; 
            for k = 1:ntr
                out.(opt.tracers{k}) = X(:,k);
            end
			diag_out = diag;
		end
	end
end

% Make plots?
if make_plots
	plot_results
end

% Save results into structure 
model.out      = out;      % tracers
model.diag     = diag_out; % diagnostics
model.grd      = grd;      % grid
model.params   = params;   % parameters
model.opt      = opt;      % options

% Save
if save_output
	disp(['Saving output/',opt.fname,'.mat']);
	save(['output/',opt.fname,'.mat'],'model');
end
if save_restart
	disp('Overwrite restart file');
	% Copy model
	MODEL = model;	
	% Only save last output, for restart conditions
	for i = 1:length(opt.tracers)		
		model.out.(opt.tracers{i}) = model.out.(opt.tracers{i})(:,end);
	end
	save(['output/omz_1d_model_restart_',opt.ecosystem,'.mat'],'model');
	% Reset
	model = MODEL;
end

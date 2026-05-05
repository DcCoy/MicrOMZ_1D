% Script to call omz_1d_model with different OM flux settings
% Model will run 20 years from previous restart with default settings

% Default settings
opt.F_poc_top_mean   = 20; % mmol C m^-2 d^-1 (mean export)
opt.F_poc_top_amp    = 0;  % relative amplitude (0.5 => ±50%)
opt.F_poc_top_period = 30; % days (period if varying)
opt.ecosystem        = 'simple_no2'; % 

% Run a spectrum of amplitude settings
if (1)
    % Good mix of lower (e.g., 120+) and higher (e.g., 14-) frequencies 
    period_range = [1 2 3 5 7 10 14 21 30 45 60 75 90 120 150 180 240 270 365];

	for i = 1:length(period_range)
		% Track progress
		disp(' ');
		disp([num2str(i),'/',num2str(length(period_range))]);
		disp(' ');
		
		% Update options
		opt.F_poc_top_amp = 0.5; % default of 50% amplitude
		opt.F_poc_top_period = period_range(i);
		if i < 10
			opt.fname = ['omz_1d_model_period_range_0',num2str(i)];
		else
			opt.fname = ['omz_1d_model_period_range_',num2str(i)];
		end

		% Run
		[~] = omz_1d_model(opt); 
	end
end

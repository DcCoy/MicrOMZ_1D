% Script to organize parameters for functional types
days2secs = 1/86400;

% ------------------------- %
% aerobic heterotroph (aer) %
% ------------------------- %
params.aer.red      = 'POC';      % reductant
params.aer.oxi      = 'O2';       % oxidant
params.aer.K_red    = 10;         % uM reductant
params.aer.K_oxi    = 0.2;        % uM oxidant
params.aer.Vmax_red = 2.4756e-05; % mol R / mol B s
params.aer.Vmax_oxi = 2.3148e-05; % mol R / mol B s
params.aer.y_red    = 4.2779;     % mol R / mol B
params.aer.y_oxi    = 4.0000;     % mol R / mol B
params.aer.e_nh4    = 0.3850;     % mol R / mol B s

% -------------------------- %
% aerobic nh4-oxidizer (aoa) %
% -------------------------- %
params.aoa.red      = 'NH4';      % reductant
params.aoa.oxi      = 'O2';       % oxidant
params.aoa.K_red    = 0.1340;     % uM reductant
params.aoa.K_oxi    = 0.3330;     % uM oxidant
params.aoa.Vmax_red = 9.4482e-05; % mol R / mol B s
params.aoa.Vmax_oxi = 1.2668e-04; % mol R / mol B s
params.aoa.y_red    = 8.1633;     % mol R / mol B
params.aoa.y_oxi    = 10.9449;    % mol R / mol B
params.aoa.e_no2    = 7.9633;     % mol R / mol B s

% -------------------------- %
% aerobic no2-oxidizer (nob) %
% -------------------------- %
params.nob.red      = 'NO2';      % reductant
params.nob.oxi      = 'O2';       % oxidant
params.nob.K_red    = 0.7780;     % uM reductant
params.nob.K_oxi    = 0.2540;     % uM oxidant
params.nob.Vmax_red = 2.7557e-04; % mol R / mol B s
params.nob.Vmax_oxi = 1.2043e-04; % mol R / mol B s
params.nob.y_red    = 15.8730;    % mol R / mol B
params.nob.y_oxi    = 6.9365;     % mol R / mol B
params.nob.e_no3    = 15.8730;    % mol R / mol B s

% ------------------------ %
% no3-reducer to no2 (nar) %
% ------------------------% 
params.nar.red      = 'POC';      % reductant
params.nar.oxi      = 'NO3';      % oxidant
params.nar.K_red    = 10;         % uM reductant
params.nar.K_oxi    = 10;         % uM oxidant
params.nar.Vmax_red = 2.4756e-05; % mol R / mol B s
params.nar.Vmax_oxi = 2.3148e-05; % mol R / mol B s
params.nar.y_red    = 6.0104;     % mol R / mol B
params.nar.y_oxi    = 12.0498;    % mol R / mol B
params.nar.e_nh4    = 0.6219;     % mol R / mol B s
params.nar.e_no2    = 12.0498;    % mol R / mol B s

% ------------------------ %
% no3-reducer to n2o (nai) %
% ------------------------% 
params.nai.red      = 'POC';      % reductant
params.nai.oxi      = 'NO3';      % oxidant
params.nai.K_red    = 10;         % uM reductant
params.nai.K_oxi    = 10;         % uM oxidant
params.nai.Vmax_red = 2.4756e-05; % mol R / mol B s
params.nai.Vmax_oxi = 2.3148e-05; % mol R / mol B s
params.nai.y_red    = 6.1822;     % mol R / mol B
params.nai.y_oxi    = 6.2258;     % mol R / mol B
params.nai.e_nh4    = 0.6454;     % mol R / mol B s
params.nai.e_n2o    = 3.1129;     % mol R / mol B s

% ----------------------- %
% no3-reducer to n2 (nao) %
% ----------------------- % 
params.nao.red      = 'POC';      % reductant
params.nao.oxi      = 'NO3';      % oxidant
params.nao.K_red    = 10;         % uM reductant
params.nao.K_oxi    = 10;         % uM oxidant
params.nao.Vmax_red = 2.4756e-05; % mol R / mol B s
params.nao.Vmax_oxi = 2.3148e-05; % mol R / mol B s
params.nao.y_red    = 6.7279;     % mol R / mol B
params.nao.y_oxi    = 5.4909;     % mol R / mol B
params.nao.e_nh4    = 0.7201;     % mol R / mol B s
params.nao.e_n2     = 2.7454;     % mol R / mol B s

% ------------------------ %
% no2-reducer to n2o (nir) %
% ------------------------ %
params.nir.red      = 'POC';      % reductant
params.nir.oxi      = 'NO2';      % oxidant
params.nir.K_red    = 10;         % uM reductant
params.nir.K_oxi    = 10;         % uM oxidant
params.nir.Vmax_red = 2.4756e-05; % mol R / mol B s
params.nir.Vmax_oxi = 2.3148e-05; % mol R / mol B s
params.nir.y_red    = 4.5994;     % mol R / mol B
params.nir.y_oxi    = 8.7517;     % mol R / mol B
params.nir.e_nh4    = 0.4290;     % mol R / mol B s
params.nir.e_n2o    = 4.3758;     % mol R / mol B s

% ----------------------- %
% no2-reducer to n2 (nio) %
% ----------------------- %
params.nio.red      = 'POC';      % reductant
params.nio.oxi      = 'NO2';      % oxidant
params.nio.K_red    = 10;         % uM reductant
params.nio.K_oxi    = 10;         % uM oxidant
params.nio.Vmax_red = 2.4756e-05; % mol R / mol B s
params.nio.Vmax_oxi = 2.3148e-05; % mol R / mol B s
params.nio.y_red    = 4.7531;     % mol R / mol B
params.nio.y_oxi    = 6.0740;     % mol R / mol B
params.nio.e_nh4    = 0.4500;     % mol R / mol B s
params.nio.e_n2     = 3.0370;     % mol R / mol B s

% ----------------------- %
% n2o-reducer to n2 (nos) %
% ----------------------- %
params.nos.red      = 'POC';      % reductant
params.nos.oxi      = 'N2O';      % oxidant
params.nos.K_red    = 10;         % uM reductant
params.nos.K_oxi    = 0.4;        % uM oxidant (NOTE: LOWER)
params.nos.Vmax_red = 2.4756e-05; % mol R / mol B s
params.nos.Vmax_oxi = 2.3148e-05; % mol R / mol B s
params.nos.y_red    = 3.2247;     % mol R / mol B
params.nos.y_oxi    = 5.5380;     % mol R / mol B
params.nos.e_nh4    = 0.2410;     % mol R / mol B s
params.nos.e_n2     = 5.5380;     % mol R / mol B s

% ---------------------------- %
% anaerobic nh4-oxidizer (aox) %
% ---------------------------- %
params.aox.red      = 'NH4';      % reductant
params.aox.oxi      = 'NO2';      % oxidant
params.aox.K_red    = 0.45;       % uM reductant
params.aox.K_oxi    = 0.45;       % uM oxidant
params.aox.Vmax_red = 4.3293e-05; % mol R / mol B s
params.aox.Vmax_oxi = 5.1257e-05; % mol R / mol B s
params.aox.y_red    = 17.7143;    % mol R / mol B
params.aox.y_oxi    = 14.9619;    % mol R / mol B
params.aox.e_no3    = 2.5714;     % mol R / mol B s
params.aox.e_n2     = 14.9524;    % mol R / mol B s

% -------------------%
% zooplankton grazer %
% -------------------%
params.zoo.u_max = 1.00*days2secs; % 1 / s (Zakem 2018)
params.zoo.m_l   = 0.00*days2secs; % 1 / s (Zakem 2018)
params.zoo.m_q   = 0.7*days2secs;  % 1 / uM C s (Zakem 2018)
params.zoo.K_B   = 1.0;            % mmol C / m3 (Zakem 2018) 
params.zoo.G_Z   = 0.5;            % percent prey to zoo biomass (Zakem 2018)
params.zoo.o2lim = 10;             % mmol O2 / m3 (Zakem 2018)

% --------------------- %
% loss parameters (all) %
% --------------------- %
params.m_l       = 0.01*days2secs; % 1 / s
params.m_q       = 0.10*days2secs; % 1 / uM C s
params.bmin      = 1e-10;          % mmol C / m3
params.loss_poc  = 0.06;           % percent losses to POC
params.graze_poc = 0;              % percent grazing to POC

% ------------------------ %
% stoichiometry parameters %
% ------------------------ %
params.CN_bio = 5;      % biomass C:N                (mol C / mol N)
params.CN_det = 117/16; % Anderson 1995 detrital C:N (mol C / mol N)

% ---------------------- % 
% Temperature dependency %
% ---------------------- % 
params.temp.TempCoeffArr = 1;
params.temp.TempAeArr    = -4000;
params.temp.TemprefArr   = 293.15;   
params.temp.Tkel         = 273.15;
params.temp.Tfunc = @(temp) params.temp.TempCoeffArr .* ...
    exp(params.temp.TempAeArr .* (1 ./ (temp + params.temp.Tkel) - 1 ./ params.temp.TemprefArr));

% OVERRIDES
hetero = {'aer','nar','nai','nao','nir','nio','nos'};
for i = 1:length(hetero)
	params.(hetero{i}).K_red = 5; % from 10
	params.(hetero{i}).K_oxi = 5; % from 10;
end
params.nos.K_oxi = 0.4; % match O2
params.m_l = 0.02*days2secs; % raise linear loss

% Test
% N2O tests
%params.nos.y_oxi = params.nai.e_n2o; % same yields in P = C equation
%params.nir.Vmax_red = 0; % kill nir
%params.nir.Vmax_oxi = 0; % kill nir
%params.nos.K_oxi = 10; % same kinetics as NO3-to-N2O

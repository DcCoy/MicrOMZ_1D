function [A_BE, rhs_mass] = build_BE_matrix(z, dz, Kz, w, dt)
% Backward-Euler matrix for 1-D diffusion + advection.
% Diffusion: Neumann (zero diffusive flux) at both ends via Ki(1)=Ki(end)=0.
% Advection: upwind in flux form at interior interfaces.
%            At the bottom face, include an explicit outflow sink when w>0
%            so advected tracers (e.g., POC) can EXIT the domain.

nz = numel(z);

% ---- Diffusivity at interfaces (harmonic mean); zero at boundaries (Neumann)
Ki = zeros(nz+1,1);
Ki(2:nz) = 2*(Kz(1:nz-1).*Kz(2:nz))./(Kz(1:nz-1)+Kz(2:nz) + eps);
Ki(1)    = 0;     % top Neumann (no diffusive flux across boundary)
Ki(end)  = 0;     % bottom Neumann

% ---- Advection velocity at interfaces (simple average to faces)
wi = zeros(nz+1,1);
wi(2:nz) = 0.5*(w(1:nz-1) + w(2:nz));
wi(1)    = w(1);
wi(end)  = w(end);

% ---- Tridiagonal coefficients for L (dC/dt = L*C + ...)
main = zeros(nz,1);
low  = zeros(nz-1,1);
upp  = zeros(nz-1,1);

for i = 1:nz
    im = i;         % index for lower interface (i-1/2)
    ip = i+1;       % index for upper interface (i+1/2)

    % Diffusive couplings (centered)
    Dm = Ki(im)/dz^2;    % to i-1
    Dp = Ki(ip)/dz^2;    % to i+1

    % -------- Upwind advection contributions --------
    % Lower interface (i-1/2): contributes to row i and (maybe) i-1
    if i>1
        wl = wi(im);
        if wl >= 0
            % inflow from below: flux uses C_{i-1}; adds to neighbor coupling
            Al =  wl/dz;   % multiplies C_{i-1}
            Bl =  0;       % multiplies C_i
        else
            % inflow from above (into lower face): flux uses C_i
            Al =  0;
            Bl = -wl/dz;   % multiplies C_i (note wl<0)
        end
    else
        Al = 0; Bl = 0;    % top boundary handled by BC
		% --- NEW: top outflow for upwelling (w < 0 at top face) ---
		wtop = wi(1);            % velocity at the TOP boundary face
		if wtop < 0              % outflow upward through the top
			% Divergence includes -F_top/dz = -(wtop*C_i)/dz
			Bl = Bl - wtop/dz;   % acts on C_1 (note: wtop < 0, so -wtop/dz > 0)
		end
    end

    % Upper interface (i+1/2): interior OR bottom face
    if i < nz
        wu = wi(ip);
        if wu >= 0
            % outflow upward from cell i: flux uses C_i
            Au = -wu/dz;   % multiplies C_i
            Bu =  0;       % multiplies C_{i+1}
        else
            % inflow from above into cell i: flux uses C_{i+1}
            Au =  0;
            Bu =  wu/dz;   % multiplies C_{i+1} (wu<0)
        end
    else
        % ---- Bottom cell: include bottom-face outflow sink when w>0 ----
        Au = 0; Bu = 0;
        wbot = wi(end);          % velocity at bottom boundary face
        if wbot > 0
            % Divergence includes -F_bottom/dz = -(wbot*C_i)/dz  -> add to main
            Au = Au - wbot/dz;   % acts on C_i (pure sink; no inflow from below)
        end
        % (If wbot<0, that would imply inflow from below; we assume none.)
    end

    % Assemble row i of L
    main(i) = -(Dp + Dm) + (Au + Bl);
    if i < nz,  upp(i)   =  Dp + Bu; end
    if i > 1,   low(i-1) =  Dm + Al; end
end

% ---- Backward-Euler system: A = I - dt * L
rhs_mass = ones(nz,1);

B = zeros(nz,3);                 % diags [-1, 0, +1]
B(:,1) = -dt * [low; 0];         % sub-diagonal
B(:,2) =  1  - dt * main;        % main diagonal
B(:,3) = -dt * [0;  upp];        % super-diagonal

A_BE = spdiags(B, -1:1, nz, nz);
end

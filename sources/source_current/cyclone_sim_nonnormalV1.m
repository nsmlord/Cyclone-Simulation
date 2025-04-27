function shallow_water_cyclone_demo_with_diagnostics()
    % SHALLOW_WATER_CYCLONE_DEMO_WITH_DIAGNOSTICS
    %
    % This code solves a shallow-water "cyclone" system on a sphere, with
    % the option to turn OFF non-conservative forcing/damping for validation
    % of approximate conservation of mass, energy, enstrophy, and angular momentum.
    %
    % To VALIDATE (and check conservation):
    %  - Set validation_mode = true, which zeroes out friction, forcing, etc.
    %  - Then the shallow-water equations are almost ideal/inviscid.
    %  - Check mass, energy, enstrophy, momentum plots after the run.

    clc; clear;

    %% USER OPTION: Toggle validation mode
    validation_mode = false;
    % - If true, the code turns OFF friction, solar heating Q, and height relaxation
    %   so that mass, energy, etc. are closer to being conserved.
    % - If false, runs the full code with friction, forcing, etc.

    %% === 1) Load Earth Globe (for display) ===
    [globe, lat_rad, lon_rad, x, y, z] = interactive_globe();

    %% === 2) Basic Parameters ===
    R        = 1.0;         % nondimensional Earth radius
    g        = 9.81;        % gravitational acceleration (m/s^2)
    Omega    = 7.2921e-5;   % Earth's rotation rate (s^-1)
    lat_res  = 2;           % grid resolution (degrees)
    lon_res  = 2;           % grid resolution (degrees)
    dLat     = deg2rad(lat_res);
    dLon     = deg2rad(lon_res);

    % Time parameters
    dt       = 0.000001;       % time step (nondimensional)
    nSteps   = 1000;       % number of simulation steps (make smaller for quick tests)
    skipFrame = 10;         % update display every few steps

    % Define the simulation grid
    latVec = deg2rad(-90 : lat_res : 90);
    lonVec = deg2rad(0   : lon_res : 360 - lon_res);
    [LAM, PHI] = meshgrid(lonVec, latVec);

    nLat = size(PHI,1);
    nLon = size(PHI,2);

    % Precompute trigonometric factors
    cosPHI = cos(PHI);
    sinPHI = sin(PHI);
    tanPHI = tan(PHI);

    %% === 3) Topography (Gaussian "mountain") ===
    lat0 = deg2rad(30);
    lon0 = deg2rad(90);
    b = 0.1 * exp(-((PHI - lat0).^2 + (LAM - lon0).^2) / 0.1);

    %% === 4) Balanced Cyclone Initial Conditions ===
    h0      = 1;              % background fluid depth
    delta_h = 0.05;           % maximum depression at cyclone center
    sigma   = 0.3;            % radial scale of the depression
    lat_c   = deg2rad(35);    % cyclone center latitude
    lon_c   = deg2rad(100);   % cyclone center longitude

    [u, v, h] = balanced_cyclone_initial_condition(PHI, LAM, h0, delta_h, sigma, lat_c, lon_c, Omega, g);

    %% === 5) Coriolis Parameter ===
    % Amplify slightly for visualization
    f = 2 * Omega * sinPHI * 3;

    %% === 6) Derivative Operators ===
    d_dlambda = @(A) (circshift(A, [0, -1]) - circshift(A, [0, 1])) ./ (2 * dLon);
    d_dphi    = @(A) (circshift(A, [-1, 0]) - circshift(A, [1, 0])) ./ (2 * dLat);

    %% === 7) Quiver3 Setup for Display ===
    [coarseLon, coarseLat] = meshgrid(lonVec, latVec);
    x_coarse = R * cos(coarseLat) .* cos(coarseLon);
    y_coarse = R * cos(coarseLat) .* sin(coarseLon);
    z_coarse = R * sin(coarseLat);

    figure(gcf);
    hold on;

    step_equator = 3;
    step_midlat  = 4;
    step_poles   = 8;

    eq_indices = find(abs(PHI(:,1)) < deg2rad(20));
    midlat_indices = find(abs(PHI(:,1)) >= deg2rad(20) & abs(PHI(:,1)) < deg2rad(60));
    pole_indices   = find(abs(PHI(:,1)) >= deg2rad(60));

    r_indices = sort([eq_indices(1:step_equator:end); 
                     midlat_indices(1:step_midlat:end);
                     pole_indices(1:step_poles:end)]);
    c_indices = 1:6:nLon;

    px = x_coarse(r_indices, c_indices);
    py = y_coarse(r_indices, c_indices);
    pz = z_coarse(r_indices, c_indices);

    px_orig = px;
    py_orig = py;
    pz_orig = pz;

    [wx, wy, wz] = velocity_sphere_to_cart(u, v, PHI, LAM);
    uq = wx(r_indices, c_indices);
    vq = wy(r_indices, c_indices);
    wq = wz(r_indices, c_indices);

    h_arrows = quiver3(px, py, pz, uq, vq, wq, 0.7, 'w', 'LineWidth', 0.8);

    %% === 8) Model Constants ===
    % Original code used friction, height relaxation, heating
    gamma       = 0.008;   % friction
    relax_rate  = 0.005;   % height relaxation
    sigma_lat   = deg2rad(15);
    sigma_lon   = deg2rad(15);

    if validation_mode
        % Turn off friction & forcing for approximate conservation checks
        gamma      = 0;    % no friction
        relax_rate = 0;    % no relaxation
        % We'll still define Q, but set it to zero each step below
    end

    % For arrow smoothing
    prev_uq = uq; 
    prev_vq = vq; 
    prev_wq = wq;
    smooth_factor = 0.85;

    %% === Storage for Diagnostics ===
    % We'll store total mass, total energy, enstrophy, momentum each step
    diag_mass      = zeros(nSteps, 1);
    diag_energy    = zeros(nSteps, 1);
    diag_enstrophy = zeros(nSteps, 1);
    diag_momentum  = zeros(nSteps, 1);

    % Start timer
    tStart = tic;

    %% === 9) Main Simulation Loop ===
    for it = 1:nSteps
        % Geopotential
        Phi_field = g * (h + b);

        % Momentum eq terms
        gradPhi_lambda = d_dlambda(Phi_field);
        gradPhi_phi    = d_dphi(Phi_field);

        damp_factor = 0.8;
        gradPhi_lambda = damp_factor * gradPhi_lambda;
        gradPhi_phi    = damp_factor * gradPhi_phi;

        adv_u = u .* d_dlambda(u) ./ (R * cosPHI) + v .* d_dphi(u) ./ R;
        adv_v = u .* d_dlambda(v) ./ (R * cosPHI) + v .* d_dphi(v) ./ R;

        geo_u =  u .* v .* tanPHI ./ R;
        geo_v = -u.^2 .* tanPHI ./ R;

        dudt = -adv_u - (1./(R*cosPHI)) .* gradPhi_lambda + f .* v - geo_u;
        dvdt = -adv_v - (1./R)          .* gradPhi_phi    - f .* u - geo_v;

        % friction
        dudt = dudt - gamma .* u;
        dvdt = dvdt - gamma .* v;

        % forcing
        solar_speed = 2*pi / 5000;
        solar_lon   = mod(solar_speed * it, 2*pi);
        solar_lat   = deg2rad(10) * cos(2*pi*it / 20000);

        Q_forcing = 0.0005 * exp(-(((PHI - solar_lat).^2)/sigma_lat^2 + ...
                                   ((LAM - solar_lon).^2)/sigma_lon^2));
        if validation_mode
            % Turn off Q in validation mode
            Q_forcing = 0;
        end

        % Symplectic Euler update for momentum
        u_new = u + dt .* dudt;
        v_new = v + dt .* dvdt;

        % velocity clamp
        max_speed = 0.2;
        speed = sqrt(u_new.^2 + v_new.^2);
        too_fast = speed > max_speed;
        if any(too_fast(:))
            scale_factor = max_speed ./ speed;
            u_new(too_fast) = u_new(too_fast) .* scale_factor(too_fast);
            v_new(too_fast) = v_new(too_fast) .* scale_factor(too_fast);
        end

        % Continuity (height)
        flux_lambda_new = h .* u_new .* cosPHI;
        flux_phi_new    = h .* v_new;
        dhdt_new = -(1./(R*cosPHI)) .* d_dlambda(flux_lambda_new) ...
                   - (1./R)         .* d_dphi(flux_phi_new);
        dhdt_new = dhdt_new + Q_forcing - relax_rate .* (h - h0);
        h_new = h + dt .* dhdt_new;

        % Minimal filtering
        polar_filter = abs(PHI) > deg2rad(60);
        if any(polar_filter(:))
            h_smooth = smoothdata(h_new, 'gaussian', 5);
            u_smooth = smoothdata(u_new, 'gaussian', 5);
            v_smooth = smoothdata(v_new, 'gaussian', 5);
            h_new(polar_filter) = h_smooth(polar_filter);
            u_new(polar_filter) = u_smooth(polar_filter);
            v_new(polar_filter) = v_smooth(polar_filter);
        end

        midlat_filter = abs(PHI) > deg2rad(30) & abs(PHI) <= deg2rad(60);
        if any(midlat_filter(:))
            h_smooth = smoothdata(h_new, 'gaussian', 3);
            u_smooth = smoothdata(u_new, 'gaussian', 3);
            v_smooth = smoothdata(v_new, 'gaussian', 3);
            blend = 0.2;
            h_new(midlat_filter) = blend*h_smooth(midlat_filter) + ...
                                   (1-blend)*h_new(midlat_filter);
            u_new(midlat_filter) = blend*u_smooth(midlat_filter) + ...
                                   (1-blend)*u_new(midlat_filter);
            v_new(midlat_filter) = blend*v_smooth(midlat_filter) + ...
                                   (1-blend)*v_new(midlat_filter);
        end

        % Update fields
        h = h_new;
        u = u_new;
        v = v_new;

        %% === DIAGNOSTICS: Mass, Energy, Enstrophy, Momentum ===
        diag_mass(it)      = compute_total_mass(h, PHI, dLat, dLon);
        diag_energy(it)    = compute_total_energy(h, u, v, g, PHI, dLat, dLon, b);
        diag_enstrophy(it) = compute_potential_enstrophy(h, u, v, f, PHI, dLat, dLon);
        diag_momentum(it)  = compute_angular_momentum(h, u, v, PHI, LAM, dLat, dLon);

        %% === Update Quiver Arrows ===
        [wx, wy, wz] = velocity_sphere_to_cart(u, v, PHI, LAM);
        uq = wx(r_indices, c_indices);
        vq = wy(r_indices, c_indices);
        wq = wz(r_indices, c_indices);

        % Smooth arrow velocities for display
        smooth_uq = smooth_factor * prev_uq + (1 - smooth_factor)*uq;
        smooth_vq = smooth_factor * prev_vq + (1 - smooth_factor)*vq;
        smooth_wq = smooth_factor * prev_wq + (1 - smooth_factor)*wq;
        prev_uq = smooth_uq;
        prev_vq = smooth_vq;
        prev_wq = smooth_wq;

        move_factor = 0.0001;
        px = px + move_factor * smooth_uq;
        py = py + move_factor * smooth_vq;
        pz = pz + move_factor * smooth_wq;

        % Project back on sphere
        mag = sqrt(px.^2 + py.^2 + pz.^2);
        px = px ./ mag;
        py = py ./ mag;
        pz = pz ./ mag;

        % Reset arrows that wander too far
        reset_dist = 0.4;
        for idx = 1:numel(px)
            dist_arrow = sqrt((px(idx) - px_orig(idx))^2 + ...
                              (py(idx) - py_orig(idx))^2 + ...
                              (pz(idx) - pz_orig(idx))^2);
            if dist_arrow > reset_dist
                px(idx) = px_orig(idx);
                py(idx) = py_orig(idx);
                pz(idx) = pz_orig(idx);
            end
        end

        % Periodic partial reset
        if mod(it, 200) == 0
            reset_percent = 0.03;
            num_to_reset = ceil(numel(px) * reset_percent);
            reset_indices = randperm(numel(px), num_to_reset);
            for idx = reset_indices
                px(idx) = px_orig(idx);
                py(idx) = py_orig(idx);
                pz(idx) = pz_orig(idx);
            end
        end

        if mod(it, skipFrame) == 0
            scale = 0.7;
            set(h_arrows, 'XData', px, 'YData', py, 'ZData', pz, ...
                'UData', scale*smooth_uq, 'VData', scale*smooth_vq, 'WData', scale*smooth_wq);

            sim_time = it * dt * 3600;
            real_time = toc(tStart);
            hrs_per_sec = sim_time / real_time;
            title(sprintf('Shallow-Water Cyclone Demo — t = %.1f hrs (%.1f sim hrs/sec)', ...
                          sim_time, hrs_per_sec), 'Color', 'w', 'FontSize', 14);
            drawnow limitrate;
        end
    end

    %% === AFTER SIM: Plot Diagnostics ===
    figure('Color','w');
    times = (1:nSteps)*dt;  % "model time" in nondim (could multiply by real scale for hours, etc.)

    % 1) Mass
    subplot(2,2,1);
    plot(times, diag_mass, 'LineWidth',1.5);
    xlabel('Time (nondim)'); ylabel('Total Mass');
    title('Total Mass vs. Time');
    grid on;

    % 2) Energy
    subplot(2,2,2);
    plot(times, diag_energy, 'LineWidth',1.5);
    xlabel('Time (nondim)'); ylabel('Total Energy');
    title('Total Energy vs. Time');
    grid on;

    % 3) Potential Enstrophy
    subplot(2,2,3);
    plot(times, diag_enstrophy, 'LineWidth',1.5);
    xlabel('Time (nondim)'); ylabel('Potential Enstrophy');
    title('Potential Enstrophy vs. Time');
    grid on;

    % 4) Angular Momentum
    subplot(2,2,4);
    plot(times, diag_momentum, 'LineWidth',1.5);
    xlabel('Time (nondim)'); ylabel('Angular Momentum (approx)');
    title('Angular Momentum vs. Time');
    grid on;

    sgtitle('Diagnostics: Mass, Energy, Enstrophy, Momentum','FontSize',14);

    % If friction or forcing is on, expect non-conservation in these plots.
    % If validation_mode = true, friction & forcing are off, so these should
    % be more conserved (though still not perfect due to numeric discretization).
end

%% --------------------------------------------------------------------- %%
function M = compute_total_mass(h, PHI, dLat, dLon)
    % Integrate h over the sphere with lat-lon weighting:
    % M = \int h * dA, with dA ~ cos(phi) dphi dlon.
    % Summation method: sum( h(i,j) * cosPHI(i,j) ) * dLat * dLon
    cosPHI = cos(PHI);
    area_elem = dLat * dLon; 
    M = sum(sum(h .* cosPHI)) * area_elem;
end

function E = compute_total_energy(h, u, v, g, PHI, dLat, dLon, b)
    % Kinetic + potential energy. In shallow-water:
    % E = ∫ [1/2 * h * (u^2+v^2) + g * h * (h/2 + b) ] dA
    % Some references do potential as 1/2*g*h^2 + g*h*b, etc. We’ll do:
    %   KE = 1/2 * h(u^2+v^2)
    %   PE = g * h(b + 1/2*h)
    % and integrate over sphere.
    cosPHI = cos(PHI);
    area_elem = dLat * dLon;

    KE = 0.5 * h .* (u.^2 + v.^2);
    PE = g * h .* (b + 0.5*h);

    density = KE + PE;
    E = sum(sum(density .* cosPHI)) * area_elem;
end

function Z = compute_potential_enstrophy(h, u, v, f, PHI, dLat, dLon)
    % A rough measure of potential enstrophy = ∫ (ζ + f)^2 / (2*h0) * h ...
    % or a simpler version: ∫ (ζ+f)^2 / (2*h) maybe. There's some variation
    % in definitions. We do a basic approach:
    %
    %   ζ = relative vorticity (on sphere approx).
    %   q = (ζ + f) / h => potential vorticity.
    %   "potential enstrophy" ~ ∫ q^2 * h dA
    %
    % We'll do a naive finite-diff for ζ. Real codes use more elaborate methods.
    %
    % Note: If friction, forcing, or smoothing is on, enstrophy won't be conserved.

    % Finite diff for vorticity in spherical coords (rough):
    %   ζ ~ (1/(R cosφ)) (∂v/∂λ) - (1/R) (∂(u cosφ)/∂φ) / cosφ
    % We skip some constants, or set R=1 for dimensionless. 
    % Let’s do a simpler approach:
    R = 1;  % from code
    dv_dlambda = (circshift(v, [0 -1]) - circshift(v, [0 1])) ./ (2 * dLon);
    % for the other part: we want ∂(u cosφ)/∂φ
    ucos = u .* cos(PHI);
    ducos_dphi = (circshift(ucos, [-1 0]) - circshift(ucos, [1 0])) ./ (2 * dLat);

    zeta = (1./(R.*cos(PHI))) .* dv_dlambda - (1./(R*cos(PHI))) .* (ducos_dphi);

    % Potential vorticity q = (ζ + f) / h
    q = (zeta + f) ./ h;
    % Potential enstrophy ~ ∫ q^2 * h dA
    ps_density = q.^2 .* h;

    cosPHI_ = cos(PHI);
    area_elem = dLat * dLon;
    Z = sum(sum(ps_density .* cosPHI_)) * area_elem;
end

function L = compute_angular_momentum(h, u, v, PHI, LAM, dLat, dLon)
    % Approximate global angular momentum
    %   \mathbf{L} = ∫ ρ h (r x u) dA
    % We'll just compute the magnitude in a rough way. Or do an "equatorial" approach.
    % A simpler approach: L_z ~ ∫ ρ h * (x u_y - y u_x) dA, etc.
    %
    % Here we do a direct "z-component" approximation: 
    %   Lz ~ ∫ (X * momentum in east-west direction, etc.)
    %
    % Because we have a 1-layer fluid, we do a naive version:
    %   Lz ~ ∫ h(φ,λ) * ( x * v - y * u ) dA
    % (This is reminiscent of ∫ (r × ρ v).z ) for R=1 sphere. 
    %
    % Real codes might do a 3D cross product, but let's keep it simple.

    x = cos(PHI).*cos(LAM);
    y = cos(PHI).*sin(LAM);

    % Momentum in the "z sense": x*v - y*u
    momZ = x .* v - y .* u;
    integrand = h .* momZ;

    cosPHI_ = cos(PHI);
    area_elem = dLat * dLon;
    L = sum(sum(integrand .* cosPHI_)) * area_elem;
end

%% ===================================================================== %%
function [u, v, h] = balanced_cyclone_initial_condition(PHI, LAM, h0, delta_h, sigma, lat_c, lon_c, Omega, g)
    f_c = 2 * Omega * sin(lat_c);  
    [nLat, nLon] = size(PHI);

    u = zeros(nLat, nLon);
    v = zeros(nLat, nLon);
    h = zeros(nLat, nLon);

    for i = 1:nLat
        for j = 1:nLon
            dlat = PHI(i,j) - lat_c;
            dlon = LAM(i,j) - lon_c;
            dlon = mod(dlon + pi, 2*pi) - pi;
            r = sqrt(dlat^2 + dlon^2);

            h(i,j) = h0 - delta_h * exp(-(r/sigma)^2);

            if r < 1e-6
                u(i,j) = 0; 
                v(i,j) = 0;
            else
                dhdr = (2 * delta_h / sigma^2) * r * exp(-(r/sigma)^2);
                % Gradient wind balance
                V = (-f_c * r + sqrt((f_c*r)^2 + 4*g*r*dhdr)) / 2;
                % Tangential direction in lat-lon coords
                u(i,j) = -V * (dlat / r);
                v(i,j) =  V * (dlon / r);
            end
        end
    end
end

%% ===================================================================== %%
function [globe, lat_rad, lon_rad, x, y, z] = interactive_globe(texture_path)
    lat_res = 0.1;
    lon_res = 0.1;
    lat = -90:lat_res:90;
    lon = 0:lon_res:360 - lon_res;
    [lon_rad, lat_rad] = meshgrid(deg2rad(lon), deg2rad(lat));

    R = 1;
    x = R * cos(lat_rad) .* cos(lon_rad);
    y = R * cos(lat_rad) .* sin(lon_rad);
    z = R * sin(lat_rad);

    if nargin < 1 || isempty(texture_path)
        texture_path = 'C:\\Users\\nisha\\OneDrive\\Desktop\\matlab\\intro_to_comp_sim\\Simulating Mid-Latitude Cyclones with the Shallow Water Model\\textures\\earth_reg_10k.jpg';
    end

    try
        earth_img = imread(texture_path);
        if size(earth_img, 3) == 1
            earth_img = repmat(earth_img, 1, 1, 3);
        end
        earth_img = imresize(earth_img, [length(lat), length(lon)]);

        figure('Color', 'k');
        globe = surf(x, y, z, flipud(earth_img), ...
                     'EdgeColor', 'none', ...
                     'FaceColor', 'texturemap');
    catch
        figure('Color', 'k');
        globe = surf(x, y, z, ...
                     'EdgeColor', 'none', ...
                     'FaceColor', [0.3, 0.5, 0.9]);
    end

    shading interp;
    axis equal off;
    light('Position', [1 0 1], 'Style', 'infinite');
    lighting gouraud;
    material dull;
    view(45, 25);
    rotate3d on;
    title('Interactive Earth Globe', 'Color', 'w', 'FontSize', 14);
end

%% ===================================================================== %%
function [vx, vy, vz] = velocity_sphere_to_cart(u, v, lat, lon)
    sinLat = sin(lat);
    cosLat = cos(lat);
    sinLon = sin(lon);
    cosLon = cos(lon);

    % Eastward unit vector
    ex = -sinLon;
    ey =  cosLon;
    ez =  0 * lat;

    % Northward unit vector
    nx = -cosLon .* sinLat;
    ny = -sinLon .* sinLat;
    nz =  cosLat;

    vx = u .* ex + v .* nx;
    vy = u .* ey + v .* ny;
    vz = u .* ez + v .* nz;
end

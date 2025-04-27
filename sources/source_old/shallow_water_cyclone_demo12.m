function shallow_water_cyclone_demo()
    clc; clear;
    
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
    dt       = 0.001;       % time step (nondimensional)
    nSteps   = 100000;      % number of simulation steps
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
    % Set background fluid depth and cyclone parameters:
    h0      = 1;              % background fluid depth
    delta_h = 0.05;           % maximum depression at cyclone center
    sigma   = 0.3;            % radial scale of the depression
    lat_c   = deg2rad(35);    % cyclone center latitude
    lon_c   = deg2rad(100);   % cyclone center longitude

    % Compute balanced initial conditions (u,v,h) using gradient wind balance.
    [u, v, h] = balanced_cyclone_initial_condition(PHI, LAM, h0, delta_h, sigma, lat_c, lon_c, Omega, g);

    %% === 5) Coriolis Parameter ===
    % Amplify slightly for visualization purposes.
    f = 2 * Omega * sinPHI * 3; 

    %% === 6) Derivative Operators in Spherical Coordinates ===
    d_dlambda = @(A) (circshift(A, [0, -1]) - circshift(A, [0, 1])) ./ (2 * dLon);
    d_dphi    = @(A) (circshift(A, [-1, 0]) - circshift(A, [1, 0])) ./ (2 * dLat);

    %% === 7) Quiver3 Setup for Display ===
    [coarseLon, coarseLat] = meshgrid(lonVec, latVec);
    x_coarse = R * cos(coarseLat) .* cos(coarseLon);
    y_coarse = R * cos(coarseLat) .* sin(coarseLon);
    z_coarse = R * sin(coarseLat);

    figure(gcf);
    hold on;
    
    % Define arrow spacing (varying by latitude for clarity)
    step_equator = 3;
    step_midlat = 4;
    step_poles = 8;
    
    eq_indices = find(abs(PHI(:,1)) < deg2rad(20));
    midlat_indices = find(abs(PHI(:,1)) >= deg2rad(20) & abs(PHI(:,1)) < deg2rad(60));
    pole_indices = find(abs(PHI(:,1)) >= deg2rad(60));
    r_indices = sort([eq_indices(1:step_equator:end); 
                     midlat_indices(1:step_midlat:end);
                     pole_indices(1:step_poles:end)]);
    c_indices = 1:6:nLon;
    
    % Get coordinates for arrow positions
    px = x_coarse(r_indices, c_indices);
    py = y_coarse(r_indices, c_indices);
    pz = z_coarse(r_indices, c_indices);
    
    % Store original positions for later resets
    px_orig = px;
    py_orig = py;
    pz_orig = pz;
    
    % Compute initial velocity vectors for arrows
    [wx, wy, wz] = velocity_sphere_to_cart(u, v, PHI, LAM);
    uq = wx(r_indices, c_indices);
    vq = wy(r_indices, c_indices);
    wq = wz(r_indices, c_indices);
    
    % Create quiver plot
    h_arrows = quiver3(px, py, pz, uq, vq, wq, 0.7, 'w', 'LineWidth', 0.8);

    %% === 8) Model Constants ===
    gamma = 0.008;       % friction (linear drag)
    relax_rate = 0.005;  % height relaxation rate
    sigma_lat = deg2rad(15); % heating zone scale in latitude
    sigma_lon = deg2rad(15); % heating zone scale in longitude

    %% === Initialize Variables for Temporal Smoothing ===
    prev_uq = uq;
    prev_vq = vq;
    prev_wq = wq;
    smooth_factor = 0.85;  % smoothing factor for arrow updates
    
    % Start the timer for simulation speed calculation
    tStart = tic;

    %% === 9) Main Simulation Loop ===
    for it = 1:nSteps
        % -- Geopotential (including topography) --
        Phi_field = g * (h + b);

        % -- Momentum Equations --
        gradPhi_lambda = d_dlambda(Phi_field);
        gradPhi_phi = d_dphi(Phi_field);
        
        % Damping on pressure gradients to control spurious features
        damp_factor = 0.8;
        gradPhi_lambda = damp_factor * gradPhi_lambda;
        gradPhi_phi = damp_factor * gradPhi_phi;
        
        % Advection terms (simplified)
        adv_u = u .* d_dlambda(u) ./ (R * cosPHI) + v .* d_dphi(u) ./ R;
        adv_v = u .* d_dlambda(v) ./ (R * cosPHI) + v .* d_dphi(v) ./ R;
        
        % Geometric terms
        geo_u = u .* v .* tanPHI ./ R;
        geo_v = - (u.^2) .* tanPHI ./ R;
        
        % Full momentum equations
        dudt = -adv_u - (1./(R * cosPHI)) .* gradPhi_lambda + f .* v - geo_u;
        dvdt = -adv_v - (1./R) .* gradPhi_phi - f .* u - geo_v;
        
        % Apply frictional damping
        dudt = dudt - gamma .* u;
        dvdt = dvdt - gamma .* v;
        
        % -- Forcing: Solar Heating & Height Relaxation --
        solar_speed = 2*pi / 5000;  % slow solar motion
        solar_lon = mod(solar_speed * it, 2*pi);
        solar_lat = deg2rad(10) * cos(2*pi*it / 20000);  
        Q = 0.0005 * exp(-(((PHI - solar_lat).^2) / sigma_lat^2 + ((LAM - solar_lon).^2) / sigma_lon^2));
        % Height relaxation term toward h0 is included in the h-update
        
        % --- Symplectic Euler Update ---
        % 1. Update momentum (u and v)
        u_new = u + dt .* dudt;
        v_new = v + dt .* dvdt;
        
        % Apply stability limit to velocity
        max_speed = 0.2;
        speed = sqrt(u_new.^2 + v_new.^2);
        too_fast = speed > max_speed;
        if any(too_fast(:))
            scale_factor = max_speed ./ speed;
            u_new(too_fast) = u_new(too_fast) .* scale_factor(too_fast);
            v_new(too_fast) = v_new(too_fast) .* scale_factor(too_fast);
        end
        
        % 2. Update fluid depth using continuity
        flux_lambda_new = h .* u_new .* cosPHI;
        flux_phi_new = h .* v_new;
        dhdt_new = -(1./(R * cosPHI)) .* d_dlambda(flux_lambda_new) - (1/R) .* d_dphi(flux_phi_new);
        dhdt_new = dhdt_new + Q - relax_rate .* (h - h0);
        h_new = h + dt .* dhdt_new;
        
        % Global smoothing every 50 steps
        if mod(it, 50) == 0
            h_new = smoothdata(h_new, 'gaussian', 3);
            u_new = smoothdata(u_new, 'gaussian', 3);
            v_new = smoothdata(v_new, 'gaussian', 3);
        end

        % Extended polar filter
        polar_filter = abs(PHI) > deg2rad(60);
        if any(polar_filter(:))
            h_smooth = smoothdata(h_new, 'gaussian', 9);
            u_smooth = smoothdata(u_new, 'gaussian', 9);
            v_smooth = smoothdata(v_new, 'gaussian', 9);
            h_new(polar_filter) = h_smooth(polar_filter);
            u_new(polar_filter) = u_smooth(polar_filter);
            v_new(polar_filter) = v_smooth(polar_filter);
        end

        % Mid-latitude smoother
        midlat_filter = abs(PHI) > deg2rad(30) & abs(PHI) <= deg2rad(60);
        if any(midlat_filter(:))
            h_smooth = smoothdata(h_new, 'gaussian', 5);
            u_smooth = smoothdata(u_new, 'gaussian', 5);
            v_smooth = smoothdata(v_new, 'gaussian', 5);
            blend = 0.3;
            h_new(midlat_filter) = blend * h_smooth(midlat_filter) + (1-blend) * h_new(midlat_filter);
            u_new(midlat_filter) = blend * u_smooth(midlat_filter) + (1-blend) * u_new(midlat_filter);
            v_new(midlat_filter) = blend * v_smooth(midlat_filter) + (1-blend) * v_new(midlat_filter);
        end
        
        % Update fields for next iteration
        h = h_new;
        u = u_new;
        v = v_new;

        % -- Update Quiver Arrows for Display --
        [wx, wy, wz] = velocity_sphere_to_cart(u, v, PHI, LAM);
        uq = wx(r_indices, c_indices);
        vq = wy(r_indices, c_indices);
        wq = wz(r_indices, c_indices);
        
        move_factor = 0.0001;
        smooth_uq = smooth_factor * prev_uq + (1-smooth_factor) * uq;
        smooth_vq = smooth_factor * prev_vq + (1-smooth_factor) * vq;
        smooth_wq = smooth_factor * prev_wq + (1-smooth_factor) * wq;
        prev_uq = smooth_uq;
        prev_vq = smooth_vq;
        prev_wq = smooth_wq;
        
        % Move arrows
        px = px + move_factor * smooth_uq;
        py = py + move_factor * smooth_vq;
        pz = pz + move_factor * smooth_wq;
        
        % Project back onto sphere
        mag = sqrt(px.^2 + py.^2 + pz.^2);
        px = px ./ mag;
        py = py ./ mag;
        pz = pz ./ mag;
        
        % Reset arrows that stray too far from original positions
        reset_dist = 0.4;
        for idx = 1:numel(px)
            dist_arrow = sqrt((px(idx) - px_orig(idx))^2 + (py(idx) - py_orig(idx))^2 + (pz(idx) - pz_orig(idx))^2);
            if dist_arrow > reset_dist
                px(idx) = px_orig(idx);
                py(idx) = py_orig(idx);
                pz(idx) = pz_orig(idx);
            end
        end
        
        % Periodically reset a small percentage of arrows
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
                'UData', scale * smooth_uq, 'VData', scale * smooth_vq, 'WData', scale * smooth_wq);
            
            % Display simulation speed in real hours per second
            sim_time = it * dt * 3600;
            real_time = toc(tStart);
            hours_per_sec = sim_time / real_time;
            title(sprintf('Shallow-Water Cyclone Demo — t = %.1f hours (%.1f sim hours/sec)', sim_time, hours_per_sec), 'Color', 'w', 'FontSize', 14);
            drawnow limitrate;
        end
    end
end

%% ===================================================================== %%
function [u, v, h] = balanced_cyclone_initial_condition(PHI, LAM, h0, delta_h, sigma, lat_c, lon_c, Omega, g)
    % Computes balanced initial fields for a cyclone based on gradient wind balance.
    %
    % Inputs:
    %   PHI, LAM  - 2D arrays of latitudes and longitudes (in radians).
    %   h0        - Background fluid depth.
    %   delta_h   - Maximum depression (positive number) at the cyclone center.
    %   sigma     - Scale of the cyclone (controls the radial extent).
    %   lat_c     - Cyclone center latitude (in radians).
    %   lon_c     - Cyclone center longitude (in radians).
    %   Omega     - Earth’s rotation rate.
    %   g         - Gravitational acceleration.
    %
    % Outputs:
    %   u, v      - Zonal and meridional velocity components.
    %   h         - Fluid depth with a Gaussian depression.
    
    % Use the Coriolis parameter at the cyclone center for balance.
    f_c = 2 * Omega * sin(lat_c);
    
    [nLat, nLon] = size(PHI);
    u = zeros(nLat, nLon);
    v = zeros(nLat, nLon);
    h = zeros(nLat, nLon);
    
    for i = 1:nLat
        for j = 1:nLon
            % Compute distance from cyclone center.
            dlat = PHI(i,j) - lat_c;
            dlon = LAM(i,j) - lon_c;
            % Adjust dlon to [-pi, pi]
            dlon = mod(dlon + pi, 2*pi) - pi;
            r = sqrt(dlat^2 + dlon^2);
            
            % Set fluid depth with a Gaussian depression:
            % h(r) = h0 - delta_h * exp[-(r/sigma)^2]
            h(i,j) = h0 - delta_h * exp(-(r/sigma)^2);
            
            if r < 1e-6
                % Avoid singularity at the center.
                u(i,j) = 0;
                v(i,j) = 0;
            else
                % Compute the radial gradient of h.
                dhdr = (2 * delta_h / sigma^2) * r * exp(-(r/sigma)^2);
                
                % Gradient wind balance:
                % V^2/r + f_c * V = g * (dh/dr)
                % Multiply through by r: V^2 + (f_c * r) * V - g*r*dhdr = 0.
                % Solve quadratic for V:
                V = (-f_c * r + sqrt((f_c * r)^2 + 4 * g * r * dhdr)) / 2;
                
                % Project the tangential velocity into u and v.
                % The counterclockwise unit tangent vector is (-dlat/r, dlon/r).
                u(i,j) = -V * (dlat / r);
                v(i,j) =  V * (dlon / r);
            end
        end
    end
end

%% ===================================================================== %%
function [globe, lat_rad, lon_rad, x, y, z] = interactive_globe(texture_path)
    % Creates a rotatable Earth globe with a texture.
    lat_res = 0.5;
    lon_res = 0.5;
    lat = -90:lat_res:90;
    lon = 0:lon_res:360 - lon_res;
    [lon_rad, lat_rad] = meshgrid(deg2rad(lon), deg2rad(lat));

    R = 1;
    x = R * cos(lat_rad) .* cos(lon_rad);
    y = R * cos(lat_rad) .* sin(lon_rad);
    z = R * sin(lat_rad);

    if nargin < 1 || isempty(texture_path)
        texture_path = 'textures/earth_reg_10k.jpg';
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
                     'FaceColor', [0.3 0.5 0.9]);
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
    % Converts (u,v) in spherical coordinates to 3D Cartesian velocities on a unit sphere.
    sinLat = sin(lat);
    cosLat = cos(lat);
    sinLon = sin(lon);
    cosLon = cos(lon);

    % Eastward unit vector
    ex = -sinLon;
    ey = cosLon;
    ez = zeros(size(lat));

    % Northward unit vector
    nx = -cosLon .* sinLat;
    ny = -sinLon .* sinLat;
    nz = cosLat;

    vx = u .* ex + v .* nx;
    vy = u .* ey + v .* ny;
    vz = u .* ez + v .* nz;
end

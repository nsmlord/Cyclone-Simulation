function shallow_water_cyclone_demo()
    clc; clear;

    %% === 1) Load Earth Globe (for display) ===
    [globe, lat_rad, lon_rad, x, y, z] = interactive_globe();

    %% === 2) Basic Parameters ===
    R        = 1.0;         % radius of Earth (nondimensional here)
    g        = 9.81;        % gravitational acceleration
    Omega    = 7.2921e-5;   % Earth's rotation rate (s^-1)
    lat_res  = 2;           % grid resolution (degrees) for simulation
    lon_res  = 2;           % grid resolution (degrees) for simulation
    dLat     = deg2rad(lat_res);
    dLon     = deg2rad(lon_res);

    % Time parameters
    dt       = 0.001;       % time step (nondimensional)
    nSteps   = 100000;      % increased for longer simulation
    skipFrame = 10;         % increased to reduce visual update frequency

    % Define the simulation grid (coarse grid)
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
    % CHANGED: Reduced mountain height and made it wider/smoother
    lat0 = deg2rad(30);
    lon0 = deg2rad(90);
    b = 0.1 * exp(-((PHI - lat0).^2 + (LAM - lon0).^2) / 0.1);  % smaller, smoother bump

    %% === 4) Initial Conditions (height & velocity) ===
    h0 = 1.0;                  
    h  = h0 * ones(nLat, nLon); % uniform fluid depth
    
    % CHANGED: Further reduced baseline flow
    u = 0.01 * ones(nLat, nLon);  % smaller zonal wind
    v = zeros(nLat, nLon);        % zero meridional wind
    
    % Create an initial cyclonic disturbance
    disturbance_lat = deg2rad(35); % Northern hemisphere
    disturbance_lon = deg2rad(100);
    radius = 0.3; % Size of disturbance
    
    for i = 1:nLat
        for j = 1:nLon
            % Calculate distance to center of disturbance
            dlat = PHI(i,j) - disturbance_lat;
            dlon = min(abs(LAM(i,j) - disturbance_lon), 2*pi - abs(LAM(i,j) - disturbance_lon));
            dist = sqrt(dlat^2 + dlon^2);
            
            if dist < radius
                % CHANGED: Further reduced strength of initial disturbance
                strength = 0.05 * (1 - dist/radius); % Weaker in center
                angle = atan2(dlat, dlon);
                
                u(i,j) = u(i,j) + strength * sin(angle);  % Add to baseline flow
                v(i,j) = v(i,j) - strength * cos(angle); 
                
                % CHANGED: Reduced height depression
                h(i,j) = h0 - 0.05 * (1 - dist/radius)^2;
            end
        end
    end

    %% === 5) Coriolis Parameter ===
    % CHANGED: Reduced amplification even more
    f = 2 * Omega * sinPHI * 3; % Less amplified for visualization

    %% === 6) Derivative Operators in Spherical Coordinates ===
    d_dlambda = @(A) (circshift(A, [0, -1]) - circshift(A, [0, 1])) ./ (2 * dLon);
    d_dphi    = @(A) (circshift(A, [-1, 0]) - circshift(A, [1, 0])) ./ (2 * dLat);

    %% === 7) Quiver3 Setup for Display ===
    % Compute sphere coordinates for the coarse simulation grid
    [coarseLon, coarseLat] = meshgrid(lonVec, latVec);
    x_coarse = R * cos(coarseLat) .* cos(coarseLon);
    y_coarse = R * cos(coarseLat) .* sin(coarseLon);
    z_coarse = R * sin(coarseLat);

    figure(gcf);  % reuse globe figure
    hold on;
    
    % Use different arrow spacing in different regions for better visualization
    step_equator = 3;   % More dense near equator
    step_midlat = 4;    % Medium density in mid-latitudes
    step_poles = 8;     % Sparse near poles
    
    % Calculate indices for each region
    eq_indices = find(abs(PHI(:,1)) < deg2rad(20));
    midlat_indices = find(abs(PHI(:,1)) >= deg2rad(20) & abs(PHI(:,1)) < deg2rad(60));
    pole_indices = find(abs(PHI(:,1)) >= deg2rad(60));
    
    % Combine all selected indices for rows
    r_indices = sort([eq_indices(1:step_equator:end); 
                     midlat_indices(1:step_midlat:end);
                     pole_indices(1:step_poles:end)]);
    
    % Use regular spacing for columns
    c_indices = 1:6:nLon;
    
    % Get coordinates for arrow positions
    px = x_coarse(r_indices, c_indices);
    py = y_coarse(r_indices, c_indices);
    pz = z_coarse(r_indices, c_indices);
    
    % Store original positions for reset
    px_orig = px;
    py_orig = py;
    pz_orig = pz;

    % Compute initial velocity vectors
    [wx, wy, wz] = velocity_sphere_to_cart(u, v, PHI, LAM);
    uq = wx(r_indices, c_indices);
    vq = wy(r_indices, c_indices);
    wq = wz(r_indices, c_indices);

    % Create quiver plot with shorter arrows and thinner lines
    h_arrows = quiver3(px, py, pz, uq, vq, wq, 0.7, 'w', 'LineWidth', 0.8);

    %% === 8) Model Constants ===
    % CHANGED: Increased stability factors
    gamma = 0.008;       % Increased friction for more stability
    relax_rate = 0.005;  % Faster relaxation to equilibrium
    sigma_lat = deg2rad(15); % Wider heating zone (smoother)
    sigma_lon = deg2rad(15); % Wider heating zone (smoother)

    %% === Initialize Variables for Temporal Smoothing ===
    prev_uq = uq;
    prev_vq = vq;
    prev_wq = wq;
    smooth_factor = 0.85;  % Increased temporal smoothing factor
    
    % Start the timer for simulation speed calculation
    tStart = tic;

    %% === 9) Main Simulation Loop ===
    for it = 1:nSteps
        % -- Geopotential --
        Phi = g * (h + b);

        % -- 9a) Continuity Equation --
        flux_lambda = h .* u .* cosPHI;
        flux_phi = h .* v;
        dhdt = -(1 ./ (R * cosPHI)) .* d_dlambda(flux_lambda) ...
               -(1 / R) .* d_dphi(flux_phi);

        % -- 9b) Momentum Equations --
        % Pressure gradient force
        gradPhi_lambda = d_dlambda(Phi);
        gradPhi_phi = d_dphi(Phi);
        
        % CHANGED: Added damping to pressure gradients to reduce mini-cyclones
        damp_factor = 0.8;
        gradPhi_lambda = damp_factor * gradPhi_lambda;
        gradPhi_phi = damp_factor * gradPhi_phi;
        
        % Advection terms (simplified for stability)
        adv_u = u .* d_dlambda(u) ./ (R * cosPHI) + v .* d_dphi(u) ./ R;
        adv_v = u .* d_dlambda(v) ./ (R * cosPHI) + v .* d_dphi(v) ./ R;
        
        % Geometric terms
        geo_u = u .* v .* tanPHI ./ R;
        geo_v = -u .* u .* tanPHI ./ R;
        
        % Full momentum equations
        dudt = -adv_u - (1./(R * cosPHI)) .* gradPhi_lambda + f .* v - geo_u;
        dvdt = -adv_v - (1./R) .* gradPhi_phi - f .* u - geo_v;

        % -- 9c) Friction (linear drag) --
        dudt = dudt - gamma .* u;
        dvdt = dvdt - gamma .* v;

        % -- 9d) Forcing: Solar Heating & Height Relaxation --
        % CHANGED: Even slower and weaker solar heating
        solar_speed = 2*pi / 5000;  % Much slower solar motion 
        solar_lon = mod(solar_speed * it, 2*pi);
        solar_lat = deg2rad(10) * cos(2*pi*it / 20000);  
        
        % CHANGED: Significantly reduced heating rate
        Q = 0.0005 * exp(-((PHI - solar_lat).^2 / sigma_lat^2 + ...
                         (LAM - solar_lon).^2 / sigma_lon^2));
        dhdt = dhdt + Q - relax_rate .* (h - h0);

        % -- 9e) Update the Fields (Forward Euler) --
        h_new = h + dt .* dhdt;
        u_new = u + dt .* dudt;
        v_new = v + dt .* dvdt;
        
        % CHANGED: Added global smoothing every 50 steps to reduce instabilities
        if mod(it, 50) == 0
            h_new = smoothdata(h_new, 'gaussian', 3);
            u_new = smoothdata(u_new, 'gaussian', 3);
            v_new = smoothdata(v_new, 'gaussian', 3);
        end
        
        % Apply stability limit to velocity (prevent extreme values)
        % CHANGED: Further reduced maximum speed
        max_speed = 0.2;
        speed = sqrt(u_new.^2 + v_new.^2);
        too_fast = speed > max_speed;
        
        if any(too_fast(:))
            scale_factor = max_speed ./ speed;
            u_new(too_fast) = u_new(too_fast) .* scale_factor(too_fast);
            v_new(too_fast) = v_new(too_fast) .* scale_factor(too_fast);
        end
        
        % CHANGED: Extended and stronger polar filter
        polar_filter = abs(PHI) > deg2rad(60); % Expanded region
        if any(polar_filter(:))
            h_smooth = smoothdata(h_new, 'gaussian', 9); % Increased smoothing
            u_smooth = smoothdata(u_new, 'gaussian', 9);
            v_smooth = smoothdata(v_new, 'gaussian', 9);
            
            h_new(polar_filter) = h_smooth(polar_filter);
            u_new(polar_filter) = u_smooth(polar_filter);
            v_new(polar_filter) = v_smooth(polar_filter);
        end
        
        % CHANGED: Added mid-latitude smoother (milder)
        midlat_filter = abs(PHI) > deg2rad(30) & abs(PHI) <= deg2rad(60);
        if any(midlat_filter(:))
            h_smooth = smoothdata(h_new, 'gaussian', 5);
            u_smooth = smoothdata(u_new, 'gaussian', 5);
            v_smooth = smoothdata(v_new, 'gaussian', 5);
            
            % Partial application (blend)
            blend = 0.3; % 30% smoothed, 70% original
            h_new(midlat_filter) = blend * h_smooth(midlat_filter) + (1-blend) * h_new(midlat_filter);
            u_new(midlat_filter) = blend * u_smooth(midlat_filter) + (1-blend) * u_new(midlat_filter);
            v_new(midlat_filter) = blend * v_smooth(midlat_filter) + (1-blend) * v_new(midlat_filter);
        end
        
        % Update fields
        h = h_new;
        u = u_new;
        v = v_new;

        % -- 9f) Update Quiver Arrows (Display) --
        % Compute 3D velocity vectors
        [wx, wy, wz] = velocity_sphere_to_cart(u, v, PHI, LAM);
        uq = wx(r_indices, c_indices);
        vq = wy(r_indices, c_indices);
        wq = wz(r_indices, c_indices);
        
        % Move arrow positions with flow (even slower movement)
        move_factor = 0.0001; % Much slower movement of particles
        
        % Apply temporal smoothing (exponential moving average)
        smooth_uq = smooth_factor * prev_uq + (1-smooth_factor) * uq;
        smooth_vq = smooth_factor * prev_vq + (1-smooth_factor) * vq;
        smooth_wq = smooth_factor * prev_wq + (1-smooth_factor) * wq;
        
        % Save current values for next iteration
        prev_uq = smooth_uq;
        prev_vq = smooth_vq;
        prev_wq = smooth_wq;
        
        % Move arrows with smoothed velocity
        px = px + move_factor * smooth_uq;
        py = py + move_factor * smooth_vq;
        pz = pz + move_factor * smooth_wq;
        
        % Project back to sphere surface
        mag = sqrt(px.^2 + py.^2 + pz.^2);
        px = px ./ mag;
        py = py ./ mag;
        pz = pz ./ mag;
        
        % Reset arrows that move too far from their original positions
        reset_dist = 0.4; % maximum distance threshold
        for i = 1:length(px(:))
            dist = sqrt((px(i) - px_orig(i))^2 + (py(i) - py_orig(i))^2 + (pz(i) - pz_orig(i))^2);
            if dist > reset_dist
                px(i) = px_orig(i);
                py(i) = py_orig(i);
                pz(i) = pz_orig(i);
            end
        end
        
        % More gradual arrow resets
        if mod(it, 200) == 0
            % Randomly reset some arrows to original positions - more gradually
            reset_percent = 0.03; % Reset just 3% of arrows
            num_to_reset = ceil(length(px(:)) * reset_percent);
            reset_indices = randperm(length(px(:)), num_to_reset);
            
            for idx = reset_indices
                px(idx) = px_orig(idx);
                py(idx) = py_orig(idx);
                pz(idx) = pz_orig(idx);
            end
        end

        if mod(it, skipFrame) == 0
            % Scale factor for velocity display
            scale = 0.7; % Scale down arrows for display
            
            set(h_arrows, 'XData', px, 'YData', py, 'ZData', pz, ...
                         'UData', scale*smooth_uq, 'VData', scale*smooth_vq, 'WData', scale*smooth_wq);
            
            % Show simulation speed in real hours per second of simulation
            sim_time = it*dt*3600; % Simulation time in hours
            real_time = toc(tStart); % Elapsed time since first tic
            
            hours_per_sec = sim_time / real_time;
            
            title(sprintf('Shallow-Water Cyclone Demo — t = %.1f hours (%.1f sim hours/sec)', ...
                 sim_time, hours_per_sec), 'Color','w','FontSize',14);
                 
            drawnow limitrate;
        end

    end % end simulation loop

end

%% ===================================================================== %%
function [globe, lat_rad, lon_rad, x, y, z] = interactive_globe(texture_path)
% INTERACTIVE_GLOBE Creates a rotatable textured Earth globe.

    % 1. Define resolution for display
    lat_res = 0.5;
    lon_res = 0.5;
    lat = -90:lat_res:90;
    lon = 0:lon_res:360 - lon_res;
    [lon_rad, lat_rad] = meshgrid(deg2rad(lon), deg2rad(lat));

    % 2. Compute sphere coordinates
    R = 1;
    x = R * cos(lat_rad) .* cos(lon_rad);
    y = R * cos(lat_rad) .* sin(lon_rad);
    z = R * sin(lat_rad);

    % 3. Load texture image
    if nargin < 1 || isempty(texture_path)
        texture_path = 'textures/earth_reg_10k.jpg';
    end
    
    try
        earth_img = imread(texture_path);
        if size(earth_img, 3) == 1
            earth_img = repmat(earth_img, 1, 1, 3);
        end
        earth_img = imresize(earth_img, [length(lat), length(lon)]);
        
        % 4. Plot globe with texture
        figure('Color', 'k');
        globe = surf(x, y, z, flipud(earth_img), ...
                     'EdgeColor', 'none', ...
                     'FaceColor', 'texturemap');
    catch
        % Fallback if texture not found
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
% Convert shallow-water (u,v) in lat-lon coordinates to 3D Cartesian (vx,vy,vz)
% for a unit sphere.
%
%   u = eastward velocity (increasing longitude)
%   v = northward velocity (increasing latitude)

    sinLat = sin(lat);   
    cosLat = cos(lat);
    sinLon = sin(lon);   
    cosLon = cos(lon);

    % Eastward unit vector
    ex = -sinLon;
    ey = cosLon;
    ez = zeros(size(lat));

    % Northward unit vector
    nx = -cosLon.*sinLat;
    ny = -sinLon.*sinLat;
    nz = cosLat;

    % Combine into 3D velocities
    vx = u .* ex + v .* nx;
    vy = u .* ey + v .* ny;
    vz = u .* ez + v .* nz;
end
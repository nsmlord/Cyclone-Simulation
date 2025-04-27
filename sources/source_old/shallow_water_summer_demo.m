function shallow_water_day_demo()
% shallow_water_day_demo simulates a regular day on Earth using a shallow‐water
% model with a diurnal solar forcing cycle. The simulation starts with a near‐balanced 
% state (nearly uniform fluid depth and weak zonal flow), and the solar heating 
% pattern moves across the globe on a 24‑hour cycle.
%
% Adjust parameters (grid resolution, friction, relaxation, and forcing amplitudes)
% to best suit your needs.

clc; clear;

%% === 1) Load Earth Globe (for display) ===
[globe, lat_rad, lon_rad, x, y, z] = interactive_globe();

%% === 2) Basic Parameters ===
R        = 1.0;         % nondimensional Earth radius
g        = 9.81;        % gravitational acceleration (m/s^2)
Omega    = 7.2921e-5;   % Earth's rotation rate (s^-1)
lat_res  = 0.5;           % grid resolution (degrees)
lon_res  = 0.5;           % grid resolution (degrees)
dLat     = deg2rad(lat_res);
dLon     = deg2rad(lon_res);

% Time parameters
dt       = 0.001;       % time step (nondimensional)
nSteps   = 36000;       % total number of simulation steps (≈ 1.5 days)
skipFrame = 10;         % update frequency for visualization

% Define simulation grid
latVec = deg2rad(-90:lat_res:90);
lonVec = deg2rad(0:lon_res:360 - lon_res);
[LAM, PHI] = meshgrid(lonVec, latVec);
nLat = size(PHI, 1);
nLon = size(PHI, 2);

% Precompute trigonometric factors
cosPHI = cos(PHI);
sinPHI = sin(PHI);
tanPHI = tan(PHI);

%% === 3) Topography (Mild Gaussian "bump") ===
% We include a mild topographic feature.
lat0 = deg2rad(30);
lon0 = deg2rad(90);
b = 0.05 * exp(-((PHI - lat0).^2 + (LAM - lon0).^2)/0.2);

%% === 4) Initial Conditions (Regular Day State) ===
% Nearly uniform fluid depth with a gentle eastward (zonal) wind.
h0 = 1.0;
h = h0 * ones(nLat, nLon);
u = 0.002 * ones(nLat, nLon);  % weak eastward flow
v = zeros(nLat, nLon);

%% === 5) Coriolis Parameter ===
% Standard Coriolis on a sphere.
f = 2 * Omega * sinPHI;

%% === 6) Derivative Operators in Spherical Coordinates ===
d_dlambda = @(A) (circshift(A, [0, -1]) - circshift(A, [0, 1])) / (2 * dLon);
d_dphi    = @(A) (circshift(A, [-1, 0]) - circshift(A, [1, 0])) / (2 * dLat);

%% === 7) Quiver3 Setup for Display ===
[coarseLon, coarseLat] = meshgrid(lonVec, latVec);
x_coarse = R * cos(coarseLat) .* cos(coarseLon);
y_coarse = R * cos(coarseLat) .* sin(coarseLon);
z_coarse = R * sin(coarseLat);

figure(gcf);
hold on;

% Define arrow spacing (varying by latitude)
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

px = x_coarse(r_indices, c_indices);
py = y_coarse(r_indices, c_indices);
pz = z_coarse(r_indices, c_indices);

% Save original arrow positions for occasional reset.
px_orig = px;
py_orig = py;
pz_orig = pz;

[wx, wy, wz] = velocity_sphere_to_cart(u, v, PHI, LAM);
uq = wx(r_indices, c_indices);
vq = wy(r_indices, c_indices);
wq = wz(r_indices, c_indices);

% Use a larger arrow scaling for visibility (adjust as needed)
h_arrows = quiver3(px, py, pz, uq, vq, wq, 2, 'r', 'LineWidth', 0.8);

%% === 8) Model Constants and Diurnal Solar Forcing ===
gamma = 0;       % friction coefficient
relax_rate = 0;  % relaxation rate for h

% Diurnal solar forcing parameters
Q0 = 0.0005;               % solar heating amplitude
% For a diurnal cycle, the subsolar point moves in longitude.
% We choose a fixed subsolar latitude (e.g., 20°N) and let the longitude vary.
solar_lat = deg2rad(20);   
sigma_lat = deg2rad(30);   % broad heating zone in latitude
sigma_lon = deg2rad(60);   % broad heating zone in longitude

% For a 24-hour cycle: real simulation time (in hours) = it*dt*3600.
% When sim_time = 24 hrs, solar_lon should complete one full revolution.
% Since 24 hrs correspond to 24 hours = 24, we have: solar_lon = mod(2*pi*(sim_time/24), 2*pi)
% Note: sim_time = it * dt * 3600, so:
%    solar_lon = mod(2*pi * it * dt * 3600/24, 2*pi) = mod(2*pi * it * dt * 150, 2*pi)
    
%% === 9) Initialize Variables for Temporal Smoothing (Arrows) ===
prev_uq = uq;
prev_vq = vq;
prev_wq = wq;
smooth_factor = 0.85;

% Start timer for simulation speed calculation
tStart = tic;

%% === 10) Main Simulation Loop ===
for it = 1:nSteps
    % Geopotential including topography
    Phi_field = g * (h + b);
    
    % Compute gradients for momentum equations
    gradPhi_lambda = d_dlambda(Phi_field);
    gradPhi_phi = d_dphi(Phi_field);
    
    % Optionally damp pressure gradients
    damp_factor = 0.9;
    gradPhi_lambda = damp_factor * gradPhi_lambda;
    gradPhi_phi = damp_factor * gradPhi_phi;
    
    % Advection terms
    adv_u = u .* d_dlambda(u) ./ (R * cosPHI) + v .* d_dphi(u) ./ R;
    adv_v = u .* d_dlambda(v) ./ (R * cosPHI) + v .* d_dphi(v) ./ R;
    
    % Geometric terms (curvature)
    geo_u = u .* v .* tanPHI ./ R;
    geo_v = -(u.^2) .* tanPHI ./ R;
    
    % Full momentum equations
    dudt = -adv_u - (1./(R * cosPHI)) .* gradPhi_lambda + f .* v - geo_u;
    dvdt = -adv_v - (1./R) .* gradPhi_phi - f .* u - geo_v;
    
    % Apply friction
    dudt = dudt - gamma .* u;
    dvdt = dvdt - gamma .* v;
    
    % --- Diurnal Solar Forcing ---
    % Compute simulation time in hours
    sim_time = it * dt * 3600;
    % Subsolar longitude: completes one full rotation (2*pi) every 24 hours.
    solar_lon = mod(2*pi * sim_time/24, 2*pi);
    % Solar heating: peak at (solar_lat, solar_lon)
    Q = Q0 * exp(-(((PHI - solar_lat).^2)/(sigma_lat^2) + ((LAM - solar_lon).^2)/(sigma_lon^2)));
    
    % Update fluid depth using continuity equation
    flux_lambda_new = h .* u .* cosPHI;
    flux_phi_new = h .* v;
    dhdt = -(1./(R * cosPHI)) .* d_dlambda(flux_lambda_new) - (1/R) .* d_dphi(flux_phi_new);
    dhdt = dhdt + Q - relax_rate .* (h - h0);
    h_new = h + dt .* dhdt;
    
    % Optional smoothing every 50 iterations to reduce noise
    if mod(it, 50) == 0
        h_new = smoothdata(h_new, 'gaussian', 3);
        u = smoothdata(u, 'gaussian', 3);
        v = smoothdata(v, 'gaussian', 3);
    end
    
    % Update fields
    h = h_new;
    u = u + dt .* dudt;
    v = v + dt .* dvdt;
    
    % Update arrow velocities for display
    [wx, wy, wz] = velocity_sphere_to_cart(u, v, PHI, LAM);
    uq = wx(r_indices, c_indices);
    vq = wy(r_indices, c_indices);
    wq = wz(r_indices, c_indices);
    
    % Temporal smoothing of arrow velocities
    smooth_uq = smooth_factor * prev_uq + (1 - smooth_factor) * uq;
    smooth_vq = smooth_factor * prev_vq + (1 - smooth_factor) * vq;
    smooth_wq = smooth_factor * prev_wq + (1 - smooth_factor) * wq;
    prev_uq = smooth_uq;
    prev_vq = smooth_vq;
    prev_wq = smooth_wq;
    
    % Update arrow positions for visualization
    move_factor = 0.0001;
    px = px + move_factor * smooth_uq;
    py = py + move_factor * smooth_vq;
    pz = pz + move_factor * smooth_wq;
    
    % Project arrow positions back onto the sphere
    mag = sqrt(px.^2 + py.^2 + pz.^2);
    px = px ./ mag;
    py = py ./ mag;
    pz = pz ./ mag;
    
    % Reset arrows that stray too far from their original positions
    reset_dist = 0.4;
    for idx = 1:numel(px)
        dist_arrow = sqrt((px(idx)-px_orig(idx))^2 + (py(idx)-py_orig(idx))^2 + (pz(idx)-pz_orig(idx))^2);
        if dist_arrow > reset_dist
            px(idx) = px_orig(idx);
            py(idx) = py_orig(idx);
            pz(idx) = pz_orig(idx);
        end
    end
    
    % Update quiver arrows every skipFrame iterations if the object still exists
    if ishandle(h_arrows) && mod(it, skipFrame) == 0
        scale = 10000;  % use the same scaling factor as initial quiver
        set(h_arrows, 'XData', px, 'YData', py, 'ZData', pz, ...
            'UData', scale * smooth_uq, 'VData', scale * smooth_vq, 'WData', scale * smooth_wq);
        
        % Display simulation speed and simulation time in hours
        real_time = toc(tStart);
        title(sprintf('Shallow-Water Day Demo — t = %.1f hours (%.1f sim hours/sec)', sim_time, sim_time/real_time), 'Color', 'w', 'FontSize', 14);
        drawnow limitrate;
    end
end

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
    nz = cos(lat);
    
    vx = u .* ex + v .* nx;
    vy = u .* ey + v .* ny;
    vz = u .* ez + v .* nz;
end

%% ===================================================================== %%
function [globe, lat_rad, lon_rad, x, y, z] = interactive_globe(texture_path)
% Creates a rotatable Earth globe with texture.
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
        if size(earth_img,3)==1
            earth_img = repmat(earth_img,1,1,3);
        end
        earth_img = imresize(earth_img, [length(lat) length(lon)]);
        figure('Color', 'k');
        globe = surf(x, y, z, flipud(earth_img), 'EdgeColor','none','FaceColor','texturemap');
    catch
        figure('Color', 'k');
        globe = surf(x, y, z, 'EdgeColor','none','FaceColor',[0.3 0.5 0.9]);
    end
    
    shading interp;
    axis equal off;
    light('Position',[1 0 1],'Style','infinite');
    lighting gouraud;
    material dull;
    view(45,25);
    rotate3d on;
    title('Interactive Earth Globe', 'Color','w','FontSize',14);
end

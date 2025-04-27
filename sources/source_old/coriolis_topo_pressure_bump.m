clc; clear;

%% === Load Earth Globe with RGB Texture ===
% This function must be in your path. It returns a textured Earth globe.
[globe, lat_rad, lon_rad, x, y, z] = interactive_globe();

%% === Parameters ===
lat_res = 2;          % degrees resolution for simulation grid
lon_res = 2;
R = 1;                % normalized radius (for the simulation, not real Earth scale)
g = 9.81;             % gravitational acceleration

% Compute approximate grid spacing (in radians * R)
dx = deg2rad(lon_res) * R;
dy = deg2rad(lat_res) * R;
dt = 0.005;           % time step (adjust for stability/visual speed)
nSteps = 5000;         % total number of time steps

%% === Topography (static terrain height) ===
% Here we create a Gaussian bump to simulate a mountain.
% Centered near 30°N, 90°E.
lat0 = deg2rad(30);
lon0 = deg2rad(90);
% Adjust the exponent denominator to set the width of the bump.
b = 0.2 * exp(-((lat_rad - lat0).^2 + (lon_rad - lon0).^2) / 0.05);

%% === Initial Conditions ===
h0 = 1.0;                     % Base fluid height
h = h0 * ones(size(lat_rad)); % Fluid height above terrain (without topography)
% You can think of the free surface as "h + b"
u = 0.1 * ones(size(lat_rad));  % Initial eastward wind (zonal)
v = zeros(size(lat_rad));       % No initial meridional wind

%% === Coriolis Force ===
Omega = 7.2921e-5; 
% We use a negative sign so that in the Northern Hemisphere an eastward flow is deflected left.
f = -2 * Omega * sin(lat_rad) * 1000;  % scaled for visibility

%% === Local Unit Vectors (for projecting 2D wind to 3D) ===
% Local east unit vector (tangent to circles of constant latitude)
ex = -sin(lon_rad); 
ey = cos(lon_rad);  
ez = zeros(size(lat_rad));
% Local north unit vector (points toward increasing latitude)
nx = -cos(lon_rad) .* sin(lat_rad);
ny = -sin(lon_rad) .* sin(lat_rad);
nz = cos(lat_rad);

%% === Derivative Functions (centered finite differences) ===
d_dx = @(A) (circshift(A, [0, -1]) - circshift(A, [0, 1])) / (2*dx);
d_dy = @(A) (circshift(A, [-1, 0]) - circshift(A, [1, 0])) / (2*dy);

%% === Downsample Grid for Visualization ===
step = 10;
r = 1:step:size(lat_rad, 1);
c = 1:step:size(lat_rad, 2);
xq = x(r, c); 
yq = y(r, c); 
zq = z(r, c);
% Initialize arrow base positions for visualization:
px = xq; 
py = yq; 
pz = zq;

% Compute initial 3D wind vectors:
wx = u .* ex + v .* nx;
wy = u .* ey + v .* ny;
wz = u .* ez + v .* nz;
uq = wx(r, c); 
vq = wy(r, c); 
wq = wz(r, c);

hold on;
% Plot initial wind arrows
h_arrows = quiver3(px, py, pz, uq, vq, wq, 0.5, 'w', 'LineWidth', 1.2);

%% === Simulation Loop ===
skipFrame = 4;  % update display every skipFrame frames
for t = 1:nSteps
    % Compute total free-surface height (fluid height + topography)
    eta = h + b;
    
    % Compute spatial gradients of eta (free surface)
    eta_x = d_dx(eta);
    eta_y = d_dy(eta);
    
    % Compute simple advection terms (centered differences)
    u_adv = u .* d_dx(u) + v .* d_dy(u);
    v_adv = u .* d_dx(v) + v .* d_dy(v);
    
    % Symplectic Euler update for u, v (momentum equations)
    u_new = u + dt * (-u_adv - g * eta_x + f .* v);
    v_new = v + dt * (-v_adv - g * eta_y - f .* u_new);
    
    % Continuity equation: update fluid depth h
    div = d_dx(h .* u_new) + d_dy(h .* v_new);
    h_new = h - dt * div;
    
    % Update fields for next time step
    u = u_new;
    v = v_new;
    h = h_new;
    
    % Project updated 2D wind to 3D using local unit vectors
    wx = u .* ex + v .* nx;
    wy = u .* ey + v .* ny;
    wz = u .* ez + v .* nz;
    uq = wx(r, c); 
    vq = wy(r, c); 
    wq = wz(r, c);
    
    % Advect arrow positions along the wind field
    px = px + 0.002 * uq;
    py = py + 0.002 * vq;
    pz = pz + 0.002 * wq;
    mag = sqrt(px.^2 + py.^2 + pz.^2);
    px = px ./ mag;
    py = py ./ mag;
    pz = pz ./ mag;
    
    % Update display every skipFrame iterations
    if mod(t, skipFrame) == 0
        set(h_arrows, 'XData', px, 'YData', py, 'ZData', pz, ...
                      'UData', uq, 'VData', vq, 'WData', wq);
        drawnow limitrate;
    end
end

%% === Topography Definition Function ===
function elevation = define_topography(lat, lon)
    % Initialize elevation field to zero.
    elevation = zeros(size(lat));
    
    % Define major mountain ranges as Gaussian bumps.
    % Format for each: {Name, Center Lat (deg), Center Lon (deg), Std Dev Lat (deg), Std Dev Lon (deg), Peak Elevation}
    ranges = {
        'Himalayas',       30,         85,          3,          10,        8849;
        'Andes',          -20,        -70,         15,           5,        6961;
        'Rockies',         45,       -110,          7,           5,        4401;
        'Alps',            46,         10,          2,           5,        4809;
        'Urals',           60,         60,          5,           2,        1895;
        'Atlas',           32,          0,          2,           5,        4167;
        'Great Dividing', -30,        150,         10,           5,        2228;
    };
    
    % Loop over each range and add its Gaussian contribution.
    for i = 1:size(ranges, 1)
        lat0 = deg2rad(ranges{i, 2});
        lon0 = deg2rad(ranges{i, 3});
        sigmaLat = deg2rad(ranges{i, 4});
        sigmaLon = deg2rad(ranges{i, 5});
        peak = ranges{i, 6};
        elevation = elevation + peak * exp(-(((lat - lat0)/sigmaLat).^2 + ((lon - lon0)/sigmaLon).^2));
    end
end

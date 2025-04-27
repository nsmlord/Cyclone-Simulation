clc; clear;

%% === Load Earth Globe with RGB Texture ===
[globe, lat_rad, lon_rad, x, y, z] = interactive_globe();

%% === Grid and Physical Parameters ===
lat_res = 3;      % degrees (should match interactive_globe)
lon_res = 3;
R = 1;            % normalized Earth radius
g = 9.81;         % gravitational acceleration

% Approximate grid spacing (note: dx ideally varies with latitude)
dx = deg2rad(lon_res) * R;
dy = deg2rad(lat_res) * R;

%% === Time Step and Skip Frame ===
dt = 0.005;       % time step (tweak as needed)
skipFrame = 1;    % update display every 'skipFrame' iterations

%% === Initial Fields ===
% Fluid depth (h): base value plus a Gaussian pressure bump.
h0 = 1;  % base fluid depth
h = h0 + 0.3 * exp( - ( (lat_rad).^2 + (lon_rad - pi).^2 )/0.2 );

% Initial wind: start with a small eastward flow.
u = 0.1 * ones(size(lat_rad));  % zonal component (eastward)
v = zeros(size(lat_rad));         % meridional component

%% === Define Coriolis Parameter ===
Omega = 7.2921e-5;              % Earth's rotation rate (rad/s)
phi = lat_rad;                  % latitude in radians
f = -2 * Omega * sin(phi);      % negative sign for leftward deflection in NH

% Scale f for visible simulation dynamics.
scale_f = 1000;  
f = f * scale_f;

%% === Finite Difference Operators (Centered Differences) ---
d_dx = @(A) (circshift(A, [0, -1]) - circshift(A, [0, 1]))/(2*dx);
d_dy = @(A) (circshift(A, [-1, 0]) - circshift(A, [1, 0]))/(2*dy);

%% === Define Local Unit Vectors (East & North) ===
% Local east unit vector (tangent to constant latitude circles)
ex = -sin(lon_rad);
ey = cos(lon_rad);
ez = zeros(size(lat_rad));

% Local north unit vector (points toward increasing latitude)
nx = -cos(lon_rad) .* sin(lat_rad);
ny = -sin(lon_rad) .* sin(lat_rad);
nz = cos(lat_rad);

%% === Project Initial Wind Field to 3D ---
% Wind3 = u*east + v*north.
wx = u .* ex + v .* nx;
wy = u .* ey + v .* ny;
wz = u .* ez + v .* nz;

%% === Downsample for Display ---
step = 5;  % adjust for arrow density
row_idx = 1:step:size(lat_rad, 1);
col_idx = 1:step:size(lat_rad, 2);

% 3D positions for arrows (initial arrow base positions)
xq = x(row_idx, col_idx);
yq = y(row_idx, col_idx);
zq = z(row_idx, col_idx);
px = xq; py = yq; pz = zq;

% Downsample initial wind vectors for arrows
uq = wx(row_idx, col_idx);
vq = wy(row_idx, col_idx);
wq = wz(row_idx, col_idx);

% Plot initial arrows
hold on;
h_arrows = quiver3(px, py, pz, uq, vq, wq, 0.5, 'w', 'LineWidth', 1.2);

%% === Main Simulation Loop ===
nSteps = 20000;
for t = 1:nSteps
    % --- Compute Pressure Gradient Terms ---
    h_x = d_dx(h);
    h_y = d_dy(h);
    
    % --- Compute Advection Terms (using centered differences) ---
    u_adv = u .* d_dx(u) + v .* d_dy(u);
    v_adv = u .* d_dx(v) + v .* d_dy(v);
    
    % --- Update Momentum Equations (Symplectic Euler) ---
    % Update u using old v and h:
    u_new = u + dt * ( - u_adv - g * h_x + f .* v );
    % Then update v using new u:
    v_new = v + dt * ( - v_adv - g * h_y - f .* u_new );
    
    % --- Update Continuity Equation (Fluid Depth h) ---
    % Compute divergence of the mass flux:
    div = d_dx(h .* u_new) + d_dy(h .* v_new);
    h_new = h - dt * div;
    
    % Update the fields
    u = u_new;
    v = v_new;
    h = h_new;
    
    % --- Project Updated 2D Wind to 3D for Visualization ---
    wx = 2*(u .* ex + v .* nx);
    wy = 2*(u .* ey + v .* ny);
    wz = 2*(u .* ez + v .* nz);
    
    % Downsample updated wind vectors:
    uq = wx(row_idx, col_idx);
    vq = wy(row_idx, col_idx);
    wq = wz(row_idx, col_idx);
    
    % --- Advect Arrow Positions (simulate particle motion) ---
    px = px + 0.002 * uq;
    py = py + 0.002 * vq;
    pz = pz + 0.002 * wq;
    
    % Reproject arrow positions onto the sphere (R = 1)
    mag = sqrt(px.^2 + py.^2 + pz.^2);
    px = px ./ mag;
    py = py ./ mag;
    pz = pz ./ mag;
    
    % --- Update Display Every 'skipFrame' Iterations ---
    if mod(t, skipFrame) == 0
        set(h_arrows, 'XData', px, 'YData', py, 'ZData', pz, ...
                      'UData', uq, 'VData', vq, 'WData', wq);
        drawnow limitrate;
    end
end

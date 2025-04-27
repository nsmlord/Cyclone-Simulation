clc; clear;

%% === Load Earth Globe with RGB Texture ===
[globe, lat_rad, lon_rad, x, y, z] = interactive_globe();

%% === Initial Wind Field ===
% Define an eastward wind (u positive means eastward).
u = sin(lat_rad);            % Eastward wind
v = zeros(size(lat_rad));    % No initial meridional wind

% Define a polar mask (|lat| >= 80°)
mask_polar = abs(lat_rad) >= deg2rad(90);

% Delete all wind at the poles
u(mask_polar) = 0;
v(mask_polar) = 0;

%% === Coriolis Parameter ===
Omega = 7.2921e-5;    % Earth's rotation rate (rad/s)
phi = lat_rad;        % latitude (radians)

% Use a Coriolis parameter that yields leftward deflection for eastward flow in NH.
% Then, zero the parameter in the polar regions.
f = -2 * Omega * sin(phi);
f(mask_polar) = 0;

% Scale f for visible simulation dynamics.
scale_f = 1000;
f = f * scale_f;

%% === Time Step ===
dt = 0.005;

%% === Define Local Unit Vectors (East & North) ===
% Local east unit vector (tangent to constant latitude circles)
ex = -sin(lon_rad);
ey = cos(lon_rad);
ez = zeros(size(lat_rad));

% Local north unit vector (points toward increasing latitude)
nx = -cos(lon_rad) .* sin(lat_rad);
ny = -sin(lon_rad) .* sin(lat_rad);
nz = cos(lat_rad);

%% === Project Wind Vectors to 3D ===
% 3D wind vector = u*east + v*north.
wx = u .* ex + v .* nx;
wy = u .* ey + v .* ny;
wz = u .* ez + v .* nz;  % initially, since v = 0, wz is zero

%% === Downsample for Display ===
step = 10;
row_idx = 1:step:size(lat_rad, 1);
col_idx = 1:step:size(lat_rad, 2);

xq = x(row_idx, col_idx);
yq = y(row_idx, col_idx);
zq = z(row_idx, col_idx);

uq = wx(row_idx, col_idx);
vq = wy(row_idx, col_idx);
wq = wz(row_idx, col_idx);

% Set initial arrow base positions.
px = xq;
py = yq;
pz = zq;

% Create a display mask so arrows are not drawn for latitudes >= 80°.
mask_display = abs(lat_rad(row_idx, col_idx)) < deg2rad(80);

%% === Plot Initial Arrows ===
hold on;
h_arrows = quiver3(px, py, pz, uq, vq, wq, 0.5, 'w', 'LineWidth', 1.2);

%% === Main Simulation Loop ===
for t = 1:300
    % --- Update wind field using symplectic Euler ---
    u_new = u + dt * (f .* v);
    v_new = v - dt * (f .* u_new);
    
    % Delete wind in polar regions after update.
    u_new(mask_polar) = 0;
    v_new(mask_polar) = 0;
    
    u = u_new;
    v = v_new;
    
    % --- Recompute 3D wind vectors ---
    wx = u .* ex + v .* nx;
    wy = u .* ey + v .* ny;
    wz = u .* ez + v .* nz;
    
    uq = wx(row_idx, col_idx);
    vq = wy(row_idx, col_idx);
    wq = wz(row_idx, col_idx);
    
    % --- Apply display mask to avoid drawing arrows in polar zones ---
    uq(~mask_display) = 0;
    vq(~mask_display) = 0;
    wq(~mask_display) = 0;
    
    % --- Advect arrow positions along the wind field ---
    px = px + 0.002 * uq;
    py = py + 0.002 * vq;
    pz = pz + 0.002 * wq;
    
    % Normalize positions to keep them on the sphere
    mag = sqrt(px.^2 + py.^2 + pz.^2);
    px = px ./ mag;
    py = py ./ mag;
    pz = pz ./ mag;
    
    % --- Update arrow graphics ---
    set(h_arrows, 'XData', px, 'YData', py, 'ZData', pz, ...
                  'UData', uq, 'VData', vq, 'WData', wq);
    
    drawnow limitrate;
end

clc; clear;

%% === Load Earth Globe with RGB Texture ===
[globe, lat_rad, lon_rad, x, y, z] = interactive_globe();

%% === Parameters ===
lat_res = 2; lon_res = 2;
R = 1; g = 9.81;
dx = deg2rad(lon_res) * R;
dy = deg2rad(lat_res) * R;
dt = 0.005;
nSteps = 3000;

%% === Topography (static terrain height) ===
lat0 = deg2rad(30);
lon0 = deg2rad(90);
b = 0.2 * exp(-((lat_rad - lat0).^2 + (lon_rad - lon0).^2) / 0.05);

%% === Initial Conditions ===
h0 = 1.0;
h = h0 * ones(size(lat_rad));
u = 0.1 * ones(size(lat_rad));  % eastward wind
v = zeros(size(lat_rad));       % no northward wind

%% === Coriolis Force ===
Omega = 7.2921e-5;
f = -2 * Omega * sin(lat_rad) * 1000;  % scaled

%% === Local Unit Vectors ===
ex = -sin(lon_rad); ey = cos(lon_rad); ez = zeros(size(lat_rad));
nx = -cos(lon_rad).*sin(lat_rad); ny = -sin(lon_rad).*sin(lat_rad); nz = cos(lat_rad);

%% === Derivative Functions ===
d_dx = @(A) (circshift(A, [0, -1]) - circshift(A, [0, 1])) / (2*dx);
d_dy = @(A) (circshift(A, [-1, 0]) - circshift(A, [1, 0])) / (2*dy);

%% === Downsample Grid for Visualization ===
step = 10;
r = 1:step:size(lat_rad, 1);
c = 1:step:size(lat_rad, 2);
xq = x(r, c); yq = y(r, c); zq = z(r, c);
px = xq; py = yq; pz = zq;

wx = u .* ex + v .* nx;
wy = u .* ey + v .* ny;
wz = u .* ez + v .* nz;
uq = wx(r, c); vq = wy(r, c); wq = wz(r, c);

hold on;
h_arrows = quiver3(px, py, pz, uq, vq, wq, 0.5, 'w', 'LineWidth', 1.2);

%% === Simulation Loop ===
skipFrame = 4;
gamma = 0.01;  % Friction coefficient

for t = 1:nSteps
    % Free surface height
    eta = h + b;
    eta_x = d_dx(eta);
    eta_y = d_dy(eta);

    % Advection terms
    u_adv = u .* d_dx(u) + v .* d_dy(u);
    v_adv = u .* d_dx(v) + v .* d_dy(v);

    % Symplectic Euler: update u, v
    u_new = u + dt * (-u_adv - g * eta_x + f .* v);
    v_new = v + dt * (-v_adv - g * eta_y - f .* u_new);

    % === Apply linear friction ===
    u_new = u_new - dt * gamma * u_new;
    v_new = v_new - dt * gamma * v_new;

    % Continuity equation
    div = d_dx(h .* u_new) + d_dy(h .* v_new);
    h_new = h - dt * div;

    % Update fields
    u = u_new;
    v = v_new;
    h = h_new;

    % Project wind to 3D
    wx = u .* ex + v .* nx;
    wy = u .* ey + v .* ny;
    wz = u .* ez + v .* nz;
    uq = wx(r, c); vq = wy(r, c); wq = wz(r, c);

    % Move arrows along flow
    px = px + 0.002 * uq;
    py = py + 0.002 * vq;
    pz = pz + 0.002 * wq;
    mag = sqrt(px.^2 + py.^2 + pz.^2);
    px = px ./ mag; py = py ./ mag; pz = pz ./ mag;

    % Update arrows every skipFrame frames
    if mod(t, skipFrame) == 0
        set(h_arrows, 'XData', px, 'YData', py, 'ZData', pz, ...
                      'UData', uq, 'VData', vq, 'WData', wq);
        drawnow limitrate;
    end
end

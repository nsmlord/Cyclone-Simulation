clc; clear;

%% === Load Earth Globe with RGB Texture ===
[globe, lat_rad, lon_rad, x, y, z] = interactive_globe();

%% === Parameters ===
lat_res = 2; lon_res = 2;
R = 1; g = 9.81;
dx = deg2rad(lon_res) * R;
dy = deg2rad(lat_res) * R;
dt = 0.005;
nSteps = 40000;

%% === Topography (one Gaussian mountain) ===
lat0 = deg2rad(30); lon0 = deg2rad(90);
b = 0.2 * exp(-((lat_rad - lat0).^2 + (lon_rad - lon0).^2) / 0.05);

%% === Initial Conditions ===
h0 = 1.0;
h = h0 * ones(size(lat_rad));
u = 0.05 * ones(size(lat_rad));  % slight eastward base wind
v = zeros(size(lat_rad));

%% === Coriolis Force ===
Omega = 7.2921e-5;
f = -2 * Omega * sin(lat_rad) * 1000;

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

%% === Simulation Settings ===
skipFrame = 4;
gamma = 0.01;       % friction
relax_rate = 0.01;  % fluid height relaxation

%% === Define Multiple Heat Sources ===
% Each row: [lat_deg, lon_deg, strength]
heat_zones = [
    15, 300, 0.005;   % Atlantic MDR
    12, 250, 0.004;   % Eastern Pacific
    15, 135, 0.004;   % Western Pacific
    12,  90, 0.003;   % Bay of Bengal
   -10,  70, 0.003;   % Southwest Indian
   -15, 160, 0.003;   % Coral Sea
];

sigma_lat = deg2rad(10);
sigma_lon = deg2rad(10);

%% === Simulation Loop ===
for t = 1:nSteps
    % --- Pressure Gradient ---
    pressure = g * (h + b);
    dp_dx = d_dx(pressure);
    dp_dy = d_dy(pressure);

    % --- Advection ---
    u_adv = u .* d_dx(u) + v .* d_dy(u);
    v_adv = u .* d_dx(v) + v .* d_dy(v);

    % --- Momentum update ---
    u_new = u + dt * (-u_adv - dp_dx + f .* v);
    v_new = v + dt * (-v_adv - dp_dy - f .* u_new);

    % --- Friction ---
    u_new = u_new - dt * gamma * u_new;
    v_new = v_new - dt * gamma * v_new;

    % --- Continuity ---
    div = d_dx(h .* u_new) + d_dy(h .* v_new);
    h_new = h - dt * div;

    % --- Thermal Forcing from All Heat Zones ---
    Q = zeros(size(h));
    for i = 1:size(heat_zones, 1)
        lat_c = deg2rad(heat_zones(i, 1));
        lon_c = deg2rad(mod(heat_zones(i, 2), 360));  % ensure 0–2π
        amp = heat_zones(i, 3);
        Q = Q + amp * exp(-((lat_rad - lat_c).^2 / sigma_lat^2 + ...
                            (lon_rad - lon_c).^2 / sigma_lon^2)) ...
              .* (1 + 0.5 * sin(2 * pi * t / 100));
    end
    h_new = h_new + dt * Q;

    % --- Relax toward base height ---
    h_new = h_new + dt * relax_rate * (h0 - h_new);

    % --- Update fields ---
    u = u_new;
    v = v_new;
    h = h_new;

    % --- Project to 3D ---
    wx = u .* ex + v .* nx;
    wy = u .* ey + v .* ny;
    wz = u .* ez + v .* nz;
    uq = wx(r, c); vq = wy(r, c); wq = wz(r, c);

    % --- Move arrows ---
    px = px + 0.002 * uq;
    py = py + 0.002 * vq;
    pz = pz + 0.002 * wq;
    mag = sqrt(px.^2 + py.^2 + pz.^2);
    px = px ./ mag; py = py ./ mag; pz = pz ./ mag;

    % --- Update Visualization ---
    if mod(t, skipFrame) == 0
        set(h_arrows, 'XData', px, 'YData', py, 'ZData', pz, ...
                      'UData', uq, 'VData', vq, 'WData', wq);
        title(sprintf('Global Cyclone Simulation — t = %.1f min', t*dt*60), ...
              'Color', 'w', 'FontSize', 14);
        drawnow limitrate;
    end
end

clc; clear;

%% === Load Earth Globe with RGB Texture ===
[globe, lat_rad, lon_rad, x, y, z] = interactive_globe();

%% === Parameters ===
lat_res = 2; lon_res = 2;
R = 1; g = 9.81;
dx = deg2rad(lon_res) * R;
dy = deg2rad(lat_res) * R;
dt = 0.000001;
nSteps = 300;

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

%% === Diagnostics ===
diag_mass = zeros(nSteps, 1);
diag_energy = zeros(nSteps, 1);
diag_enstrophy = zeros(nSteps, 1);
diag_momentum = zeros(nSteps, 1);
diag_max_speed = zeros(nSteps, 1);
dLat = deg2rad(lat_res);
dLon = deg2rad(lon_res);
area_elem = dLat * dLon;
cosphi = cos(lat_rad);

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

    % Apply linear friction
    u_new = u_new - dt * gamma * u_new;
    v_new = v_new - dt * gamma * v_new;

    % Continuity equation
    div = d_dx(h .* u_new) + d_dy(h .* v_new);
    h_new = h - dt * div;

    % Update fields
    u = u_new;
    v = v_new;
    h = h_new;

    % === DIAGNOSTICS ===
    diag_mass(t) = sum(h(:) .* cosphi(:)) * area_elem;

    KE = 0.5 * h .* (u.^2 + v.^2);
    PE = g * h .* (b + 0.5 * h);
    diag_energy(t) = sum((KE(:) + PE(:)) .* cosphi(:)) * area_elem;

    dv_dx = d_dx(v);
    ucos = u .* cosphi;
    ducos_dy = d_dy(ucos);
    zeta = (1 ./ cosphi) .* dv_dx - (1 ./ cosphi) .* ducos_dy;
    q = (zeta + f) ./ h;
    ps_density = q.^2 .* h;
    diag_enstrophy(t) = sum(ps_density(:) .* cosphi(:)) * area_elem;

    x_ = cos(lat_rad) .* cos(lon_rad);
    y_ = cos(lat_rad) .* sin(lon_rad);
    momZ = x_ .* v - y_ .* u;
    diag_momentum(t) = sum((h(:) .* momZ(:)) .* cosphi(:)) * area_elem;

    diag_max_speed(t) = max(sqrt(u.^2 + v.^2), [], 'all');

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

%% === Plot Diagnostics with Labels ===
figure;
times = (1:nSteps) * dt;  % nondimensional time

subplot(3,2,1);
plot(times, diag_mass, 'LineWidth', 1.5);
xlabel('Time (non-dim)'); ylabel('Total Mass');
title('Mass Conservation'); grid on;

subplot(3,2,2);
plot(times, diag_energy, 'LineWidth', 1.5);
xlabel('Time (non-dim)'); ylabel('Total Energy');
title('Energy (Kinetic + Potential)'); grid on;

subplot(3,2,3);
plot(times, diag_enstrophy, 'LineWidth', 1.5);
xlabel('Time (non-dim)'); ylabel('Potential Enstrophy');
title('Potential Enstrophy'); grid on;

subplot(3,2,4);
plot(times, diag_momentum, 'LineWidth', 1.5);
xlabel('Time (non-dim)'); ylabel('Angular Momentum');
title('Zonal Angular Momentum'); grid on;

subplot(3,2,5);
plot(times, diag_max_speed, 'LineWidth', 1.5);
xlabel('Time (non-dim)'); ylabel('Max Wind Speed');
title('Maximum Wind Speed'); grid on;

sgtitle('Simulation Diagnostics Over Time');

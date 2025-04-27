clc; clear;

% === Load Earth Globe with RGB texture ===
[globe, lat_rad, lon_rad, x, y, z] = interactive_globe();

% === Initial toy wind field ===
v_x = -sin(lat_rad);  % zonal
v_y = zeros(size(lat_rad));  % meridional

% Convert to 3D vectors
u3 = -v_x .* sin(lon_rad);
v3 =  v_x .* cos(lon_rad);
w3 =  zeros(size(lat_rad));  % horizontal only

% Downsample for display
step = 30;
row_idx = 1:step:size(lat_rad, 1);
col_idx = 1:step:size(lat_rad, 2);

xq = x(row_idx, col_idx);
yq = y(row_idx, col_idx);
zq = z(row_idx, col_idx);

uq = u3(row_idx, col_idx);
vq = v3(row_idx, col_idx);
wq = w3(row_idx, col_idx);

% Initialize arrow positions
px = xq;
py = yq;
pz = zq;

% Plot moving arrows
hold on;
h_arrows = quiver3(px, py, pz, uq, vq, wq, 0.5, 'w', 'LineWidth', 1.2);

% Animate flow by updating arrow origins
for t = 1:300
    % Fake wind evolution (cosine pulse)
    v_x = -sin(lat_rad) .* cos(0.01 * t);
    v_y = zeros(size(lat_rad));

    % Update 3D vector field
    u3 = -v_x .* sin(lon_rad);
    v3 =  v_x .* cos(lon_rad);

    % Resample
    uq = u3(row_idx, col_idx);
    vq = v3(row_idx, col_idx);

    % Move starting points along flow (simulate particle motion)
    px = px + 0.002 * uq;
    py = py + 0.002 * vq;
    pz = pz + 0.002 * wq;

    % Optional: wrap to surface radius
    R = 1;
    mag = sqrt(px.^2 + py.^2 + pz.^2);
    px = px ./ mag * R;
    py = py ./ mag * R;
    pz = pz ./ mag * R;

    % Update arrows
    set(h_arrows, 'XData', px, 'YData', py, 'ZData', pz, ...
                  'UData', uq, 'VData', vq, 'WData', wq);

    drawnow limitrate;
end

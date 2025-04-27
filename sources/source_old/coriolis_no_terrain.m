    clc; clear;
    
    %% === Load Earth Globe with RGB Texture ===
    [globe, lat_rad, lon_rad, x, y, z] = interactive_globe();
    
    %% === Initial Wind Field (Shallow Water Momentum with Coriolis Only) ===
    % u: zonal (east-west) component, v: meridional (north-south) component.
    % Here we initialize a simple eastward flow that will be modified by Coriolis.
    u = -sin(lat_rad);       % initial zonal wind (eastward if negative here)
    v = zeros(size(lat_rad));  % no initial meridional wind
    
    %% === Define Coriolis Parameters ===
    Omega = 7.2921e-5;          % Earth's rotation rate (rad/s)
    f = 2 * Omega * sin(lat_rad);  % Coriolis parameter (varies with latitude)
    
    % For visual demonstration, scale f so that the induced rotation is visible.
    scale_f = 1000;   % (Experiment with this value)
    f = f * scale_f;
    
    %% === Time Step for Symplectic Euler Integration ===
    dt = 0.005;  % time step (in simulation units)
    
    %% === Convert Wind Field to 3D Vectors for Visualization ===
    % Here we convert the 2D wind components (u,v) on the sphere to 3D vectors.
    u3 = -u .* sin(lon_rad);
    v3 =  u .* cos(lon_rad);
    w3 = zeros(size(lat_rad));  % purely horizontal
    
    %% === Downsample for Display (to avoid cluttering) ===
    step = 30;
    row_idx = 1:step:size(lat_rad, 1);
    col_idx = 1:step:size(lat_rad, 2);
    
    xq = x(row_idx, col_idx);
    yq = y(row_idx, col_idx);
    zq = z(row_idx, col_idx);
    
    uq = u3(row_idx, col_idx);
    vq = v3(row_idx, col_idx);
    wq = w3(row_idx, col_idx);
    
    % Initialize arrow starting positions (these will be advected by the flow)
    px = xq;
    py = yq;
    pz = zq;
    
    %% === Plot Moving Arrows ===
    hold on;
    h_arrows = quiver3(px, py, pz, uq, vq, wq, 0.5, 'w', 'LineWidth', 1.2);
    
    %% === Main Simulation Loop (Symplectic Euler with Coriolis) ===
    for t = 1:300
        % --- Symplectic Euler Update ---
        % Update u using old v:
        u_new = u + dt * (f .* v);  % u_{n+1} = u_n + dt * (f*v_n)
        % Update v using the new u:
        v_new = v - dt * (f .* u_new);  % v_{n+1} = v_n - dt * (f*u_{n+1})
        
        % Assign updated values back to u and v
        u = u_new;
        v = v_new;
        
        % --- Convert Updated 2D Wind to 3D Vectors ---
        u3 = -u .* sin(lon_rad);
        v3 =  u .* cos(lon_rad);
        % w3 remains zero (horizontal flow only)
        
        % Downsample updated wind vectors for display:
        uq = u3(row_idx, col_idx);
        vq = v3(row_idx, col_idx);
        wq = w3(row_idx, col_idx);
        
        % --- Advect Arrow Starting Positions (simulate particle motion) ---
        % Here, we "move" the arrow origins along the wind field.
        px = px + 0.002 * uq;
        py = py + 0.002 * vq;
        pz = pz + 0.002 * wq;
        
        % Reproject the arrow positions onto the sphere (R = 1)
        mag = sqrt(px.^2 + py.^2 + pz.^2);
        px = px ./ mag;
        py = py ./ mag;
        pz = pz ./ mag;
        
        % --- Update Arrow Graphics ---
        set(h_arrows, 'XData', px, 'YData', py, 'ZData', pz, ...
                      'UData', uq, 'VData', vq, 'WData', wq);
        
        drawnow limitrate;
    end

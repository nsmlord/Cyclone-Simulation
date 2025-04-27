function shallow_water_cyclone_demo_with_diagnostics()
    % SHALLOW_WATER_CYCLONE_DEMO_WITH_DIAGNOSTICS
    % Solves shallow-water cyclone with diagnostics and now energy from forcing.

    clc; clear;

    %% USER OPTIONS
    study_mode = true;  
    params.validation_mode = ~study_mode;  
    params.gamma      = 1;    % friction coefficient
    params.relax_rate = 1;    % height relaxation rate
    params.Q_amp      = 1;    % solar heating amplitude
    params.h0         = 1;    % background fluid depth

    %% 1) Load Earth Globe
    [globe, lat_rad, lon_rad, x, y, z] = interactive_globe();
    hold on;

    %% 2) Basic Parameters
    R       = 1;          
    g       = 9.81;       
    Omega   = 7.2921e-5;  
    lat_res = 2;          
    lon_res = 2;          
    dLat    = deg2rad(lat_res);
    dLon    = deg2rad(lon_res);
    area_elem = dLat * dLon;

    dt       = 1e-6;    
    nSteps   = 30000;   
    skipFrame= 1;       

    latVec = deg2rad(-90:lat_res:90);
    lonVec = deg2rad(0:lon_res:360-lon_res);
    [LAM, PHI] = meshgrid(lonVec, latVec);

    cosPHI = cos(PHI);
    tanPHI = tan(PHI);

    %% 3) Topography
    lat0 = deg2rad(30);
    lon0 = deg2rad(90);
    b = 0.1 * exp(-((PHI-lat0).^2 + (LAM-lon0).^2)/0.1);
    if params.validation_mode
        b(:) = 0;
    end

    %% 4) Initial Conditions
    [u, v, h] = balanced_cyclone_initial_condition(PHI, LAM, params.h0, 0.05, 0.3, deg2rad(35), deg2rad(100), Omega, g);

    %% 5) Coriolis Parameter
    f = 2 * Omega * sin(PHI) * 3000;

    %% 6) Derivative Operators
    d_dlambda = @(A) (circshift(A,[0,-1]) - circshift(A,[0,1]))/(2*dLon);
    d_dphi    = @(A) (circshift(A,[-1,0]) - circshift(A,[1,0]))/(2*dLat);

    %% 7) Quiver3 Setup
    [cLon,cLat] = meshgrid(lonVec, latVec);
    xq = R*cos(cLat).*cos(cLon);
    yq = R*cos(cLat).*sin(cLon);
    zq = R*sin(cLat);

    eq_idx  = find(abs(PHI(:,1))<deg2rad(20));
    mid_idx = find(abs(PHI(:,1))>=deg2rad(20)&abs(PHI(:,1))<deg2rad(60));
    po_idx  = find(abs(PHI(:,1))>=deg2rad(60));
    r_idx   = sort([eq_idx(1:3:end); mid_idx(1:4:end); po_idx(1:8:end)]);
    c_idx   = 1:6:size(LAM,2);

    px = xq(r_idx,c_idx);
    py = yq(r_idx,c_idx);
    pz = zq(r_idx,c_idx);

    [WX,WY,WZ] = velocity_sphere_to_cart(u, v, PHI, LAM);
    h_arrows = quiver3(px,py,pz,WX(r_idx,c_idx),WY(r_idx,c_idx),WZ(r_idx,c_idx),...
                       0.7,'w','LineWidth',0.8);

    %% 8) Model Constants
    if params.validation_mode
        gamma      = 0;
        relax_rate = 0;
    else
        gamma      = params.gamma;
        relax_rate = params.relax_rate;
    end

    prev_uq = WX(r_idx,c_idx);
    prev_vq = WY(r_idx,c_idx);
    prev_wq = WZ(r_idx,c_idx);
    smooth_factor = 0.85;

    %% Storage for Diagnostics
    diag_mass      = zeros(nSteps,1);
    diag_energy    = zeros(nSteps,1);
    diag_enstrophy = zeros(nSteps,1);
    diag_momentum  = zeros(nSteps,1);
    diag_max_speed = zeros(nSteps,1);

    power_solar = zeros(nSteps,1);
    power_relax = zeros(nSteps,1);
    power_drag  = zeros(nSteps,1);

    energy_solar = zeros(nSteps,1);
    energy_relax = zeros(nSteps,1);
    energy_drag  = zeros(nSteps,1);

    tStart = tic;

    %% 9) Main Simulation Loop
    for it = 1:nSteps
        Phi = g*(h + b);

        if params.validation_mode
            Qf = 0;
        else
            ss = 2*pi/5000;
            lonS = mod(ss*(it*dt),2*pi);
            latS = deg2rad(10)*cos(2*pi*(it*dt)/20000);
            Qf = params.Q_amp * exp(-(((PHI-latS).^2)/deg2rad(15)^2 + ((LAM-lonS).^2)/deg2rad(15)^2));
        end

        % Instantaneous Power
        power_solar(it) = sum(sum((Phi .* Qf) .* cosPHI)) * area_elem;
        power_relax(it) = sum(sum((Phi .* (-relax_rate.*(h-params.h0))) .* cosPHI)) * area_elem;
        power_drag(it)  = sum(sum((-gamma.*h.*(u.^2 + v.^2)) .* cosPHI)) * area_elem;

        % Accumulate into Energies
        if it == 1
            energy_solar(it) = power_solar(it) * dt;
            energy_relax(it) = power_relax(it) * dt;
            energy_drag(it)  = power_drag(it)  * dt;
        else
            energy_solar(it) = energy_solar(it-1) + power_solar(it) * dt;
            energy_relax(it) = energy_relax(it-1) + power_relax(it) * dt;
            energy_drag(it)  = energy_drag(it-1) + power_drag(it)  * dt;
        end

        % Momentum Equations
        gradLon = d_dlambda(Phi);
        gradLat = d_dphi(Phi);

        adv_u = u.*d_dlambda(u)./(R*cosPHI) + v.*d_dphi(u)./R;
        adv_v = u.*d_dlambda(v)./(R*cosPHI) + v.*d_dphi(v)./R;
        geo_u = u.*v.*tanPHI./R;
        geo_v = -u.^2.*tanPHI./R;

        dudt = -adv_u - gradLon./(R*cosPHI) + f.*v - geo_u - gamma.*u;
        dvdt = -adv_v - gradLat./R           - f.*u - geo_v - gamma.*v;

        u = u + dt .* dudt;
        v = v + dt .* dvdt;

        % Velocity Clamp
        max_speed = 0.2;
        speed = sqrt(u.^2 + v.^2);
        mask = speed > max_speed;
        u(mask) = u(mask).*max_speed./speed(mask);
        v(mask) = v(mask).*max_speed./speed(mask);

        % Continuity Equation
        fluxLon = h.*u.*cosPHI;
        fluxLat = h.*v;
        dhdt = -(1./(R*cosPHI)).*d_dlambda(fluxLon) - (1./R).*d_dphi(fluxLat);
        dhdt = dhdt + Qf - relax_rate.*(h-params.h0);
        h = h + dt .* dhdt;

        % Smoothing
        if ~params.validation_mode
            polar = abs(PHI)>deg2rad(60);
            if any(polar(:))
                hs = smoothdata(h,'gaussian',5);
                hs(polar)=hs(polar);
                h(polar)=hs(polar);
            end
            mid = abs(PHI)>deg2rad(30)&abs(PHI)<=deg2rad(60);
            if any(mid(:))
                hs = smoothdata(h,'gaussian',3);
                h(mid)=0.2*hs(mid)+0.8*h(mid);
            end
        end

        % Diagnostics
        diag_mass(it)      = compute_total_mass(h,PHI,dLat,dLon);
        diag_energy(it)    = compute_total_energy(h,u,v,g,PHI,dLat,dLon,b);
        diag_enstrophy(it) = compute_potential_enstrophy(h,u,v,f,PHI,dLat,dLon);
        diag_momentum(it)  = compute_angular_momentum(h,u,v,PHI,LAM,dLat,dLon);
        diag_max_speed(it) = max(sqrt(u.^2+v.^2),[],'all');

        % Quiver update
        if mod(it,skipFrame)==0
            [WX,WY,WZ] = velocity_sphere_to_cart(u,v,PHI,LAM);
            smooth_u = smooth_factor*prev_uq + (1-smooth_factor)*WX(r_idx,c_idx);
            smooth_v = smooth_factor*prev_vq + (1-smooth_factor)*WY(r_idx,c_idx);
            prev_uq = smooth_u; prev_vq = smooth_v;
            set(h_arrows,'UData',smooth_u,'VData',smooth_v,'WData',WZ(r_idx,c_idx));
            drawnow limitrate;
        end
    end

    %% 10) Plot Diagnostics
    times = linspace(0,0.1,nSteps);

    figure('Color','w');
    subplot(3,2,1);
    plot(times,diag_mass,'LineWidth',1.5); title('Total Mass'); grid on;
    xlabel('Time (nondim)'); xlim([0 0.1]);
    subplot(3,2,2);
    plot(times,diag_energy,'LineWidth',1.5); title('Total Energy'); grid on;
    xlabel('Time (nondim)'); xlim([0 0.1]);
    subplot(3,2,3);
    plot(times,diag_enstrophy,'LineWidth',1.5); title('Potential Enstrophy'); grid on;
    xlabel('Time (nondim)'); xlim([0 0.1]);
    subplot(3,2,4);
    plot(times,diag_momentum,'LineWidth',1.5); title('Angular Momentum'); grid on;
    xlabel('Time (nondim)'); xlim([0 0.1]);
    subplot(3,2,5);
    plot(times,diag_max_speed,'LineWidth',1.5); title('Max Wind Speed'); grid on;
    xlabel('Time (nondim)'); xlim([0 0.1]);

    figure('Color','w');
    subplot(3,1,1);
    plot(times,energy_solar,'r-','LineWidth',1.5); title('Cumulative Energy from Solar Heating'); grid on;
    xlabel('Time (nondim)'); xlim([0 0.1]);
    subplot(3,1,2);
    plot(times,energy_relax,'g-','LineWidth',1.5); title('Cumulative Energy from Relaxation'); grid on;
    xlabel('Time (nondim)'); xlim([0 0.1]);
    subplot(3,1,3);
    plot(times,energy_drag,'b-','LineWidth',1.5); title('Cumulative Energy from Drag Dissipation'); grid on;
    xlabel('Time (nondim)'); xlim([0 0.1]);
end

%% ===================================================================== %%
function M = compute_total_mass(h, PHI, dLat, dLon)
    cosPHI = cos(PHI);
    area_elem = dLat * dLon; 
    M = sum(sum(h .* cosPHI)) * area_elem;
end

function E = compute_total_energy(h, u, v, g, PHI, dLat, dLon, b)
    cosPHI = cos(PHI);
    area_elem = dLat * dLon;
    KE = 0.5 * h .* (u.^2 + v.^2);
    PE = g * h .* (b + 0.5 * h);
    density = KE + PE;
    E = sum(sum(density .* cosPHI)) * area_elem;
end

function Z = compute_potential_enstrophy(h, u, v, f, PHI, dLat, dLon)
    R = 1;  % nondimensional
    dv_dlambda = (circshift(v, [0, -1]) - circshift(v, [0, 1])) ./ (2 * dLon);
    ucos = u .* cos(PHI);
    ducos_dphi = (circshift(ucos, [-1, 0]) - circshift(ucos, [1, 0])) ./ (2 * dLat);
    zeta = (1 ./ (R .* cos(PHI))) .* dv_dlambda - (1 ./ (R * cos(PHI))) .* ducos_dphi;
    q = (zeta + f) ./ h;
    ps_density = q.^2 .* h;
    cosPHI_ = cos(PHI);
    area_elem = dLat * dLon;
    Z = sum(sum(ps_density .* cosPHI_)) * area_elem;
end

function L = compute_angular_momentum(h, u, v, PHI, LAM, dLat, dLon)
    x = cos(PHI) .* cos(LAM);
    y = cos(PHI) .* sin(LAM);
    momZ = x .* v - y .* u;
    integrand = h .* momZ;
    cosPHI_ = cos(PHI);
    area_elem = dLat * dLon;
    L = sum(sum(integrand .* cosPHI_)) * area_elem;
end

%% ===================================================================== %%
function [u, v, h] = balanced_cyclone_initial_condition(PHI, LAM, h0, delta_h, sigma, lat_c, lon_c, Omega, g)
    f_c = 2 * Omega * sin(lat_c);  
    [nLat, nLon] = size(PHI);
    u = zeros(nLat, nLon);
    v = zeros(nLat, nLon);
    h = zeros(nLat, nLon);
    for i = 1:nLat
        for j = 1:nLon
            dlat = PHI(i,j) - lat_c;
            dlon = LAM(i,j) - lon_c;
            dlon = mod(dlon + pi, 2*pi) - pi;
            r = sqrt(dlat^2 + dlon^2);
            h(i,j) = h0 - delta_h * exp(-(r/sigma)^2);
            if r < 1e-6
                u(i,j) = 0; 
                v(i,j) = 0;
            else
                dhdr = (2 * delta_h / sigma^2) * r * exp(-(r/sigma)^2);
                V = (-f_c * r + sqrt((f_c * r)^2 + 4 * g * r * dhdr)) / 2;
                u(i,j) = -V * (dlat / r);
                v(i,j) =  V * (dlon / r);
            end
        end
    end
end

%% ===================================================================== %%
function [globe, lat_rad, lon_rad, x, y, z] = interactive_globe(texture_path)
    lat_res = 0.3;
    lon_res = 0.3;
    lat = -90:lat_res:90;
    lon = 0:lon_res:360 - lon_res;
    [lon_rad, lat_rad] = meshgrid(deg2rad(lon), deg2rad(lat));
    R = 1;
    x = R * cos(lat_rad) .* cos(lon_rad);
    y = R * cos(lat_rad) .* sin(lon_rad);
    z = R * sin(lat_rad);
    if nargin < 1 || isempty(texture_path)
        texture_path = 'C:\\Users\\nisha\\OneDrive\\Desktop\\matlab\\intro_to_comp_sim\\Simulating Mid-Latitude Cyclones with the Shallow Water Model\\textures\\earth_reg_10k.jpg';
    else
        disp('fuck');
    
    end
    try
        earth_img = imread(texture_path);
        if size(earth_img, 3) == 1
            earth_img = repmat(earth_img, 1, 1, 3);
        end
        earth_img = imresize(earth_img, [length(lat), length(lon)]);
        figure('Color', 'k');
        globe = surf(x, y, z, flipud(earth_img), 'EdgeColor', 'none', 'FaceColor', 'texturemap');
    catch
        figure('Color', 'k');
        globe = surf(x, y, z, 'EdgeColor', 'none', 'FaceColor', [0.3, 0.5, 0.9]);
    end
    shading interp;
    axis equal off;
    light('Position', [1 0 1], 'Style', 'infinite');
    lighting gouraud;
    material dull;
    view(160, 20);
    rotate3d on;
    title('Interactive Earth Globe', 'Color', 'w', 'FontSize', 14);
end

%% ===================================================================== %%
function [vx, vy, vz] = velocity_sphere_to_cart(u, v, lat, lon)
    sinLat = sin(lat);
    cosLat = cos(lat);
    sinLon = sin(lon);
    cosLon = cos(lon);
    ex = -sinLon;
    ey = cosLon;
    ez = zeros(size(lat));
    nx = -cosLon .* sinLat;
    ny = -sinLon .* sinLat;
    nz = cosLat;
    vx = u .* ex + v .* nx;
    vy = u .* ey + v .* ny;
    vz = u .* ez + v .* nz;
end

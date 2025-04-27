function [globe, lat_rad, lon_rad, x, y, z] = interactive_globe(texture_path)
% INTERACTIVE_GLOBE Creates a rotatable textured Earth globe.
%
% USAGE:
%   [globe, lat_rad, lon_rad, x, y, z] = interactive_globe();
%   [globe, lat_rad, lon_rad, x, y, z] = interactive_globe('textures/earth_reg_10k.jpg');
%
% INPUT:
%   texture_path (optional): path to the Earth texture image (RGB equirectangular)
%
% OUTPUT:
%   globe      - Handle to the surf plot (for updating)
%   lat_rad    - Latitude grid (radians)
%   lon_rad    - Longitude grid (radians)
%   x, y, z    - Cartesian coordinates of the sphere

    % 1. Define resolution of the grid
    lat_res = 0.5;
    lon_res = 0.5;

    % 2. Generate lat-lon grid
    lat = -90:lat_res:90;
    lon = 0:lon_res:360 - lon_res;
    [lon_rad, lat_rad] = meshgrid(deg2rad(lon), deg2rad(lat));

    % 3. Compute sphere coordinates
    R = 1;
    x = R * cos(lat_rad) .* cos(lon_rad);
    y = R * cos(lat_rad) .* sin(lon_rad);
    z = R * sin(lat_rad);

    % 4. Load texture image
    if nargin < 1 || isempty(texture_path)
        texture_path = 'textures/earth_reg_10k.jpg';  % default path
    end
    earth_img = imread(texture_path);
    
    % If grayscale, convert to RGB
    if size(earth_img, 3) == 1
        earth_img = repmat(earth_img, 1, 1, 3);
    end

    % Resize texture to match mesh
    earth_img = imresize(earth_img, [length(lat), length(lon)]);

    % 5. Plot globe
    figure('Color', 'k');
    globe = surf(x, y, z, flipud(earth_img), ...
                 'EdgeColor', 'none', ...
                 'FaceColor', 'texturemap');

    % 6. Styling
    shading interp;
    axis equal off;
    %colormap('gray');  % can be removed if using RGB texture

    % 7. Lighting
    light('Position', [1 0 1], 'Style', 'infinite');
    lighting gouraud;
    material dull;

    % 8. View
    view(45, 25);
    rotate3d on;

    % 9. Title
    title('Interactive Earth Globe', 'Color', 'w', 'FontSize', 14);
end

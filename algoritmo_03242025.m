% === 1. Leer datos desde el archivo NetCDF ===
filename = 'A3.nc';          % Nombre del archivo
time_index = 24;             % Tiempo a analizar

% Obtener dimensiones
info = ncinfo(filename, 'U');
x_len = info.Size(1);
z_len = info.Size(3);

% Parámetros espaciales
dx = 10;                     % Espaciado en X (m)
dz = 10;                     % Espaciado en Z (m)
x = (0:x_len-2) * dx;        % Eje X (por promediado de U)
z = (0:z_len-1) * dz;        % Eje Z
[X, Z] = meshgrid(x, z);     % Malla 2D

% Inicializar campos
Q_field = zeros(z_len, x_len-1);
U_field = zeros(z_len, x_len-1);
W_field = zeros(z_len, x_len-1);

% === 2. Cálculo de Q y velocidades para cada altura ===
for z_level = 1:z_len
    % Leer U y W
    u = squeeze(ncread(filename, 'U', [1, 1, z_level, time_index], [Inf, 1, 1, 1]));
    w = squeeze(ncread(filename, 'W', [1, 1, z_level, time_index], [Inf, 1, 1, 1]));

    % Promediar U para igualar tamaño con W
    u = 0.5 * (u(1:end-1) + u(2:end));

    % Derivadas
    du_dx = diff(u) / dx;
    du_dx(end+1) = du_dx(end);

    dw_dz = diff(w) / dz;
    dw_dz(end+1) = dw_dz(end);

    % Estimar V por continuidad
    divergence = du_dx + dw_dz;
    v = -cumsum(divergence) * dx;

    dv_dx = diff(v) / dx;
    dv_dx(end+1) = dv_dx(end);
    dv_dz = zeros(size(dv_dx));  % despreciamos

    % Tensor y vorticidad
    Sxx = du_dx;
    Syy = dv_dz;
    Sxy = 0.5 * (du_dx + dv_dx);
    Omega = 0.5 * (du_dx - dv_dx);

    % Criterio Q
    Q = 0.5 * (Omega.^2 - (Sxx.^2 + 2*Sxy.^2 + Syy.^2));

    % Guardar resultados
    Q_field(z_level, :) = Q;
    U_field(z_level, :) = u;
    W_field(z_level, :) = w;
end

% === 3. Visualización ===
figure;

% Mostrar Q con más contraste
contourf(X, Z, Q_field, 100, 'LineColor', 'none');
colormap(turbo);  % Colores más nítidos
colorbar;
caxis([-1e-4 1e-4]);  % Ajustar escala para mejor visualización

hold on;

% Mostrar campo de velocidades con escala amplificada
quiver(X, Z, 10*U_field, 10*W_field, 'k');  % Multiplicamos por 10

xlabel('X (m)');
ylabel('Z (m)');
title(['Criterio Q y campo de velocidades - t = ' num2str(time_index)]);
grid on;
ylim([0 400]);  % Limitar Z para enfocarnos donde hay actividad

% === 4. Detección de centros de vórtices ===
Q_threshold = 1e-5;  % Umbral mínimo

% Verificar si se tiene la función imregionalmax
if exist('imregionalmax', 'file')
    vortex_mask = imregionalmax(Q_field) & (Q_field > Q_threshold);

    % Coordenadas
    [vortex_z_idx, vortex_x_idx] = find(vortex_mask);
    vortex_x = x(vortex_x_idx);
    vortex_z = z(vortex_z_idx);

    % Dibujar centros
    plot(vortex_x, vortex_z, 'ro', 'MarkerSize', 8, 'LineWidth', 2);

    % Etiquetas
    for i = 1:length(vortex_x)
        text(vortex_x(i)+20, vortex_z(i), sprintf('V%d', i), ...
            'Color', 'w', 'FontWeight', 'bold', 'FontSize', 10);
    end
else
    disp('Advertencia: No se encontró la función imregionalmax. No se detectarán centros de vórtices automáticamente.');
end

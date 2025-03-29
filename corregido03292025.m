% === 1. Leer datos desde NetCDF ===
filename = 'A3.nc';
time_index = 24;

% Constantes de Titán
g_titan = 1.352;  % gravedad en m/s^2

% Leer variables
u_all  = squeeze(ncread(filename, 'U',  [1, 1, 1, time_index], [Inf, 1, Inf, 1]));
w_all  = squeeze(ncread(filename, 'W',  [1, 1, 1, time_index], [Inf, 1, Inf, 1]));
ph     = squeeze(ncread(filename, 'PH',  [1, 1, 1, time_index], [Inf, 1, Inf, 1]));
phb    = squeeze(ncread(filename, 'PHB', [1, 1, 1, time_index], [Inf, 1, Inf, 1]));

% Altura real
altura = (ph + phb) / g_titan;
altura = 0.5 * (altura(:,1:end-1) + altura(:,2:end));  % promedio entre niveles

% Promediar U y W
u = 0.5 * (u_all(1:end-1,:) + u_all(2:end,:));
w = 0.5 * (w_all(:,1:end-1) + w_all(:,2:end));

% Malla espacial
[nx, nz] = size(u);
dx = 10;
x = (0:nx-1) * dx;
z = mean(altura, 1);

% === 2. Derivadas ===
[du_dx, du_dz] = gradient(u, dx, mean(diff(z)));
[dw_dx, dw_dz] = gradient(w, dx, mean(diff(z)));

divergencia = du_dx + dw_dz;
v = -cumsum(divergencia, 1) * dx;

[dv_dx, dv_dz] = gradient(v, dx, mean(diff(z)));

% === 3. Tensor y Q ===
Sxx = du_dx;
Syy = dv_dz;
Sxy = 0.5 * (du_dz + dv_dx);
Omega = 0.5 * (du_dz - dv_dx);

Q = 0.5 * (Omega.^2 - (Sxx.^2 + 2*Sxy.^2 + Syy.^2));
Q_plot = Q.^2;

% === 4. Asegurar tamaños ===
[X, Z] = meshgrid(x, z);
if ~isequal(size(Q_plot), size(Z)), Q_plot = Q_plot'; end
if ~isequal(size(Q), size(Z)), Q = Q'; end

% === 5. Umbral y regiones ===
umbral = 1e-5;
Q_bin = Q_plot > umbral;

[etiquetas, num] = bwlabel(Q_bin, 8);
area_vortices = zeros(1, num);
centros_x = zeros(1, num);
centros_z = zeros(1, num);

fprintf('\nVÓRTICES DETECTADOS: %d\n', num);

for i = 1:num
    region_mask = etiquetas == i;
    area_vortices(i) = sum(region_mask(:)) * dx * mean(diff(z));
    
    % Obtener máximos dentro de esta región
    Q_max = imregionalmax(Q_plot) & region_mask;
    [rz, rx] = find(Q_max);
    
    if ~isempty(rx)
        centros_x(i) = mean(x(rx));
        centros_z(i) = mean(z(rz));
    else
        % Si no hay máximos, usamos el centroide aproximado
        [rz2, rx2] = find(region_mask);
        centros_x(i) = mean(x(rx2));
        centros_z(i) = mean(z(rz2));
    end
    
    fprintf('  Vórtice %d: Centro (X=%.0f m, Z=%.0f m), Área ≈ %.0f m²\n', ...
        i, centros_x(i), centros_z(i), area_vortices(i));
end

% === 6. Visualización ===
figure('Color', 'w');
qmin = 0;
qmax = 3e-4;

contourf(X, Z, Q_plot, 50, 'LineStyle', 'none');
colormap(turbo);
cb = colorbar;
ylabel(cb, 'Criterio Q² (1/s⁴)');
caxis([qmin qmax]);
hold on;

% Contornos Q original
contour(X, Z, Q, [0 0], 'k', 'LineWidth', 1);
contour(X, Z, Q, [0.5e-4, 1e-4, 2e-4], 'LineColor', [0.5 0 0], 'LineWidth', 0.8);
contour(X, Z, Q, [-2e-4, -1e-4, -0.5e-4], 'LineColor', [0 0 0.5], 'LineWidth', 0.8);

% Marcar centros de vórtices
plot(centros_x, centros_z, 'wo', 'MarkerFaceColor', 'k', 'MarkerSize', 6);

% Etiquetas y leyenda
xlabel('Distancia X (m)');
ylabel('Altura real Z (m)');
title(sprintf('Criterio Q² y Vórtices en Titán - Tiempo %d', time_index));
grid on;
legend({'Q = 0', 'Q positivo', 'Q negativo', 'Centro de vórtice'}, ...
    'Location', 'northeastoutside');

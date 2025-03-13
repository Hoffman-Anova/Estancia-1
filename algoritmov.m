% === 1. Leer datos desde el archivo NetCDF ===
filename = 'A3.nc';  

% Parámetros para seleccionar el tiempo y nivel vertical
time_index = 24;
z_level = 1;

% Leer la componente U (x-wind)
u = squeeze(ncread(filename, 'U', [1, 1, z_level, time_index], [Inf, 1, 1, 1]));

% Leer la componente W (z-wind)
w = squeeze(ncread(filename, 'W', [1, 1, z_level, time_index], [Inf, 1, 1, 1]));

% === 2. Ajustar tamaños usando la media ===
% Promediar U para igualar tamaño con W (por ser 'stagged')
u = 0.5 * (u(1:end-1) + u(2:end));  % Promediar U para que mida igual que W

% === 3. Definir la malla espacial ===
dx = 10; % Espaciado en x
dz = 10; % Espaciado en z (ajustar si es necesario)
x = (0:(size(u, 1)-1)) * dx;

% === 4. Calcular gradientes ===
du_dx = diff(u) / dx;
du_dx = [du_dx; du_dx(end)]; % Mantener tamaño

dw_dz = diff(w) / dz;
dw_dz = [dw_dz; dw_dz(end)]; % Mantener tamaño

% === 5. Estimar la componente V usando la ecuación de continuidad ===
divergence = du_dx + dw_dz;
v = -cumsum(divergence) * dx;  % Integración simple con respecto a x

% === 6. Calcular gradientes de V ===
dv_dx = diff(v) / dx;
dv_dx = [dv_dx; dv_dx(end)]; % Mantener tamaño

dv_dz = zeros(size(dv_dx));  % Suponiendo constante en Z

% === 7. Calcular los términos del tensor de gradiente de velocidad ===
Sxx = du_dx;
Syy = dv_dz;
Sxy = 0.5 * (du_dx + dv_dx);

Omega = 0.5 * (du_dx - dv_dx);  % Vorticidad

% === 8. Calcular el criterio Q ===
Q = 0.5 * (Omega.^2 - (Sxx.^2 + 2*Sxy.^2 + Syy.^2));

% === 9. Visualización ===
figure;
plot(x, Q);
title('Criterio Q - Detección de Vórtices');
xlabel('X');
ylabel('Q');
grid on;

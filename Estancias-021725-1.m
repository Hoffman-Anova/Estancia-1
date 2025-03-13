% Limpieza del espacio de trabajo
clear all; close all; clc

% Nombre del archivo a leer
archivo = 'A3.nc';

% Variables a extraer del archivo
Var1 = 'U';    % Velocidad horizontal
Var2 = 'T';    % Temperatura
Var3 = 'W';    % Velocidad vertical
Var4 = 'HFX';  % Flujo de calor

% Lectura de datos del archivo NetCDF
u = squeeze(ncread(archivo, Var1));    % Velocidad horizontal
t = squeeze(ncread(archivo, Var2));    % Temperatura
w = squeeze(ncread(archivo, Var3));    % Velocidad vertical
HFX = squeeze(ncread(archivo, Var4));  % Flujo de calor

% Lectura de variables de altura
PHB = squeeze(ncread(archivo, 'PHB')); % Base state geopotential
PH = squeeze(ncread(archivo, 'PH'));   % Perturbation geopotential

% Cálculo de altura
gravedad = 9.81;                       % Aceleración gravitacional estándar
gz = PHB + PH;                         % Geopotencial total
gz2d = mean(gz, 3, "omitmissing");     % Promedio en la dimensión temporal
ejez = mean(gz2d)./gravedad;           % Conversión a metros

% Selección del step temporal específico (23)
tiempo_especifico = 23;

% Extracción de datos para el tiempo específico
temp_slice = t(:,:,tiempo_especifico);
u_slice = u(:,:,tiempo_especifico);
w_slice = w(:,:,tiempo_especifico);

% Creación de la figura con dos subplots
figure('Position', [100 100 800 800]);

% Subplot 1: Mapa de calor de temperatura
subplot(2,1,1)
contourf(temp_slice, 20);  % 20 niveles de contorno
colorbar;
title('Distribución de Temperatura - Step 23', 'FontSize', 14);
xlabel('Posición Horizontal (puntos de malla)', 'FontSize', 12);
ylabel('Altura (puntos de malla)', 'FontSize', 12);
colormap(jet);  % Esquema de colores
grid on;

% Subplot 2: Perfil vertical
subplot(2,1,2)
posicion_x = round(size(temp_slice,1)/2);  % Posición central
plot(temp_slice(posicion_x,:), ejez(1:end-1), 'r-', 'LineWidth', 2);
title('Perfil Vertical de Temperatura', 'FontSize', 14);
xlabel('Temperatura (K)', 'FontSize', 12);
ylabel('Altura (m)', 'FontSize', 12);
ylim([0 max(ejez)]);
grid on;

% Ajustar espaciado entre subplots
set(gcf, 'Color', 'white');  % Fondo blanco
sgtitle('Análisis de Temperatura - A3.nc', 'FontSize', 16);

% Cálculos estadísticos básicos
temp_mean = mean(temp_slice, 'all');
temp_std = std(temp_slice, [], 'all');
temp_max = max(temp_slice, [], 'all');
temp_min = min(temp_slice, [], 'all');

% Mostrar estadísticas en la consola
fprintf('\nEstadísticas de temperatura para step 23:\n');
fprintf('Temperatura media: %.2f K\n', temp_mean);
fprintf('Desviación estándar: %.2f K\n', temp_std);
fprintf('Temperatura máxima: %.2f K\n', temp_max);
fprintf('Temperatura mínima: %.2f K\n', temp_min);
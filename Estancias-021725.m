% PASO 1: Limpieza del espacio de trabajo
clear all; close all; clc

% PASO 2: Definición del archivo y variables
archivo = 'A3.nc';
tiempo_especifico = 23;  % Step temporal a analizar

% PASO 3: Lectura de variables
% Cuando usamos squeeze, la dimensión y (que es 1) se elimina
% Por lo tanto, nuestras variables tendrán dimensiones diferentes
T = squeeze(ncread(archivo, 'T'));     % Ahora será 2000 x 99 x 43
P = squeeze(ncread(archivo, 'P'));     % 2000 x 99 x 43
PB = squeeze(ncread(archivo, 'PB'));   % 2000 x 99 x 43

% Variables para altura
PH = squeeze(ncread(archivo, 'PHB'));  % 2000 x 100 x 43
PHB = squeeze(ncread(archivo, 'PH'));  % 2000 x 100 x 43

gravedad = 1.352;
gz = PHB + PH;  % Geopotencial total
gz2d = gz(:,:,tiempo_especifico); % Seleccionar tiempo específico
ejez = gz2d./gravedad;  % Conversión a metros

% Extraer temperatura para el tiempo específico
temp_slice = T(:,:,tiempo_especifico);

%Ajustar niveles verticales
% Como T tiene 99 niveles y PH/PHB tienen 100, calculamos altura en los niveles de T
ejez_t = (ejez(:,1:end-1) + ejez(:,2:end))/2;

% PASO 7: Creación de la figura
posicion_x = 2000/2;  % Punto medio aproximado (de 2000 puntos)

% PASO 8: Graficar el perfil
plot(temp_slice(posicion_x,:), mean(ejez_t));

% PASO 9: Personalización de la gráfica

xlabel('Temperatura (K)');
ylabel('Altura (m)');
ylim([0 max(max(ejez_t))]);
grid on;
set(gcf, 'Color', 'white');


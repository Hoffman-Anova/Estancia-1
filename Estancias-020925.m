% Read data
filename = 'A3.nc';
ph = ncread(filename,'PH');   
phb = ncread(filename,'PHB'); 
temp = ncread(filename,'T');  

% Calculate actual temperature
theta0 = 100;  
actual_temp = temp + theta0;  
actual_temp = actual_temp - 273.15;  

% Calculate height
geopotential = (ph + phb);    
height = geopotential / 1.352;

% Create x-axis
x = linspace(0, 2000*10, 2000);

% Get slice for specific time step
time_step = 1;
temp_slice = squeeze(actual_temp(:,1,:,time_step));
height_slice = squeeze(height(:,1,1:end-1,time_step));



% Create the plot
figure
pcolor(x, height_slice', temp_slice')  % Quitamos la transpuesta aquí porque ya la hicimos antes
shading interp
colorbar
xlabel('Distance (m)')
ylabel('Height (m)')
title('Actual Temperature (°C)')
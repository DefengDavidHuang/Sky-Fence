clear all

plane_height = 1.0e3; %1.0e4;
orbit_height = 5.5e5; %5.5e7; %5.5e6; % %3.0e5;
r_Earth = 6.371e6;

T_satellite = 2*pi * sqrt((orbit_height + r_Earth)^3/((5.972*6.6723) * 1.0e13));

plane_speed = 250; %-125% -250; %-600; %-50; %600; %250; %0; % %-600; %-50; %  -250; % %600; % %600; % % %600; %200;   0; %50; % %250m/second = 900km/hour

%{
time_stepzoom = 0.00000002;
total_timezoom = 0.0006; % % second0.6; %  
start_timezoom = -0.0003; % %-0.3; %
end_timezoom = start_timezoom + total_timezoom;

time_step = time_stepzoom;
total_time = total_timezoom;
start_time = start_timezoom;
end_time = end_timezoom;
%}

time_step = 0.0005; %05;
total_time = 0.6; %0.1; %0.6; %0.6; %0.0006; % % second  
start_time = -0.3; %-0.05; %0.3; %-0.3; %-0.0003; % %
end_time = start_time + total_time;

%Nx = 4001; %0; %00; %1000;
%Ny = 4001; %0; %00; %1000;
%N_F = Nx * Ny + 1;
N_F = 401; %001;
%%%%%%%%%

%antenna size 0.4m x 0.4m
%antenna_halfwidth = 0.22; %0.2;
%antenna_halflength = 0.22; % 0.2;

%p1 = 0.1;
%p2 = 0.1;

wave_select = 1; %2; %3; %
wave_set = [1.0e-2, 2.0e-2, 0.5e-2];
band = [30, 15, 60];
wavelength =  wave_set(wave_select); %1.0e-2; %30GHZ 0.5e-2; %60GHz %2.0e-2; %15GHZ %

width_y = 20; %10; %4; %2; %0; %0; %10; %20; %20; %20.0;
length_x = 20; %10; %4; % 2; %0; %10; %0; %20; %20.0;    %x: parallel to ground and perpendicular to the satellite link

Num_data_point = round(total_time/time_step,0); 
%
z = plane_height;

y = -100:0.1:100;

ytemp = [-width_y/2:width_y/N_F:width_y/2]';
         
 
U_N = 1 / (wavelength * z) * abs(sum(exp(-1i * 2 * pi / wavelength / z * (kron(y, ytemp))))) * width_y /N_F * length_x; %x = 0
          
         
lambda = wavelength;
b = length_x;
a = width_y;
%z = dPlane;
         
beta= 0; %b*X/(lambda*z); 
gamma= a*y/(lambda*z); %a*Y/(lambda*z); 
U_F = 1 / (wavelength * z) * a * b * abs(sinc(gamma));

%figure
  
plot(y, U_N,'-xr', y, U_F,'-ok');
legend('Numerical', 'Frauhofer');
 
%title(strcat(strcat(num2str(band(wave_select)),strcat(strcat(" GHz    ", num2str(round(theta_initial_degree,2))), ' degree')), strcat("   ", strcat(num2str(plane_speed),' m/sec')))); 
 
 %save('FH40degree.mat','x_time','y_FH_F','y_FH');
  %alpha_initial = alpha_initial - alpha_between_satellite;
%end

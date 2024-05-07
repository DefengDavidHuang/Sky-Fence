clear all

%also add Doppler measurements

plane_height = 1.0e4; %1.0e3; %1.0e4;
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

time_step = 0.00051; %0.000051; % %0.0000051; %05;
total_time = 0.6; %0.3; %5e-3; % %0.1; %0.6; %0.6; %0.0006; % % second  
start_time = -0.3; %-0.15; %-0.05; % %-0.05; %0.3; %-0.3; %-0.0003; % %
end_time = start_time + total_time;

Nx = 401; %0; %00; %1000;
Ny = 401; %0; %00; %1000;
N_F = Nx * Ny + 1;
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

width_z = 10; %20; %10; %4; %2; %0; %0; %10; %20; %20; %20.0;
length_x = 10; %10; %4; % 2; %0; %10; %0; %20; %20.0;    %x: parallel to ground and perpendicular to the satellite link

Num_data_point = length([start_time : time_step : end_time]); %round(total_time/time_step,0); 
%

%alpha_station = d_station/r_Earth; % angle (relative to center of Earth) between two neighbour ground stations

%distanceX = zeros(Num_data_point+1,1); % + 10; % the movement of flying object in X direction

%alpha_initial = 0.089999999999999999999999999999999999; %elevation angle: 40 degrees

%alpha_initial = 0; %elevation angle: 90 degrees

alpha_initial = 0.03674999999999; %elevation angle: 65 degrees

%[1/60:-1/120:-1/60]*pi; %[1/30:-1/120:-1/30]*pi; % [1/60:-1/2400:-1/60]*pi; %  %[1/300:-1/600:-1/300]*pi
theta_ini = pi/2 - atan(((orbit_height + r_Earth)*sin(alpha_initial))/(((orbit_height+r_Earth)*cos(alpha_initial)) - r_Earth));
theta_initial_degree = theta_ini/pi * 180

A = plane_height/sin(theta_ini) * sqrt((orbit_height + r_Earth)^2 - (r_Earth * cos(theta_ini))^2 ) - r_Earth * sin(theta_ini);

temp = (orbit_height +  r_Earth);
         
dthetadt =  (temp * r_Earth * cos(alpha_initial) - temp * temp)/(temp * temp - 2 * temp * r_Earth * cos(alpha_initial) + r_Earth * r_Earth)  * (2 * pi / T_satellite);

%Note: the speed of the airplane is in the x direction (for ground static
%coordinate)
vsquare = plane_speed * plane_speed + (plane_height/sin(theta_ini) * dthetadt)^2; 
doppler_slope = vsquare/wavelength/(plane_height/sin(theta_ini))

%doppler = - plane_speed * cos(theta_ini)/wavelength
%doppler_offset = plane_speed * plane_speed * total_time/2/wavelength/plane_height*sin(theta_ini);

%var_theta = zeros(Num_groundstation,length(alpha_initial));

%intensity_all_output1 = zeros(Num_satellite,Num_groundstation,Num_data_point+1);
 
output_intensity = zeros(Num_data_point,1);
output_FK = output_intensity;
output_FN = output_intensity;
output_FN_F = output_intensity;
output_FH_F = output_intensity;
output_FH = output_intensity;
% distanceY = ([-total_time/2 : time_step : total_time/2]) * plane_speed; % - 1000;

%for m_S = 1:1:Num_satellite
     
    % alpha = -1.0 * [start_time : time_step : end_time] * 2 * pi/T_satellite + alpha_initial;   %the alpha angle relative to the first ground station. Note the minus is for, at the beginning the alpha is positive, which corrsponds to negative time
     %distanceY1st = -1.0 * ([start_time : time_step : end_time]) * plane_speed - d_satellite * (m_S -1); % relative to the first ground station, the minus sign due to the distance at the beginning is set to be positive
     
     distanceX1st = -1.0 * ([start_time : time_step : end_time]) * plane_speed; % relative to the satellite orbit plane, the minus sign due to the distance at the beginning is set to be positive
     
     %m_S
     
     %theta_plane = pi/2 - atan(((orbit_height + r_Earth + plane_height)*sin(alpha_initial))/(((orbit_height + r_Earth + plane_height)*cos(alpha_initial)) - r_Earth - plane_height))
     %Y_satellite_position = plane_height/tan(theta_plane)
     
    % for m_Y = 1:1:Num_groundstation
     
      %   m_Y   
         
     distanceX = distanceX1st; % - Gstation_location(m_Y); %distance relative to the line perpendicular to the ground station of interest
     
     alpha = -1.0 * [start_time : time_step : end_time] * 2 * pi/T_satellite + alpha_initial; % - Gstation_location(m_Y)/r_Earth;
    
     theta = pi/2 - atan(((orbit_height+r_Earth).*sin(alpha))./(((orbit_height+r_Earth).*cos(alpha)) - r_Earth)); %the elevation angle of satellite relative to the ground station of interest
     
     %theta_planetemp = pi/2 - atan(((orbit_height+r_Earth + plane_height).*sin(alpha))./(((orbit_height + r_Earth + plane_height).*cos(alpha)) - r_Earth)); %the elevation angle of satellite relative to the ground station
     %Y_satellite_positiontemp = plane_height./tan(theta_planetemp)
     
     %theta(Num_data_point+1)/pi*180
     
     test1 = strcat(num2str(band(wave_select)),'GHz');
     test2 = strcat(num2str(round(theta_initial_degree,2)), 'degree');
     test3 = strcat(test1,test2);
     test4 = strcat(num2str(plane_speed),'mps');
     test5 = strcat(test3,test4);
     total_time1 = round(total_time*10)
     test51 = strcat(test5,num2str(total_time1));
     savefilename = strcat(test51,'Cross.mat');
     
     %testtitle = strcat(,strcat(strcat("GHz", num2str(round(theta_initial_degree,2))), 'degree')), strcat(strcat(num2str(plane_speed),'.mat'))
     
       for m = 1:1:Num_data_point
         %m
         width_ztemp = width_z * sin(theta(m));
         %distancez = plane_height/sin(theta(m)) + cos(theta(m)) * (distanceY(m)  - plane_height / tan(theta(m))); %+ plane_height /tan(theta_ini));
         
         dPlane = plane_height/sin(theta_ini) * cos(theta(m) - theta_ini); %projected onto the optical axis, distance between flying object and groud receiver 
         dLEO = sqrt((orbit_height + r_Earth)^2 - (r_Earth * cos(theta(m)))^2 ) - r_Earth * sin(theta(m));     % distance between satellite and ground station
         
         dPlane2 = dPlane * dPlane;
         dLEO2 = (dLEO - dPlane) * (dLEO - dPlane);
                
        % funFresnel = @(r1,r2) A*exp(1i * 2 * pi / wavelength * (r1 + r2))./(r1 .* r2) ;
         
         tempX = distanceX(m);
         tempZ = plane_height/sin(theta_ini) * sin(theta(m) - theta_ini); %sin(theta(m)) * (distanceX(m)  - plane_height / tan(theta(m)));   %deviation of flying object center from the optical axis
         %fun1 = @(x1,z1) funFresnel(sqrt(x1.*x1 + (z1 + tempY).*(z1 + tempY) + dPlane2), sqrt(x1 .* x1 + (z1 + tempY) .* (z1 + tempY) + dLEO2));
         
         %q1 = quad2d(fun1,-length_x/2,length_x/2,-width_ytemp/2,width_ytemp/2); %, 'RelTol', 0.1); %'AbsTol',1); %1e-1);
         
         %q1 =integral2(fun1,-length_x/2,length_x/2,-width_ytemp/2,width_ytemp/2);
         
         pathloss = (orbit_height*orbit_height)/(dLEO * dLEO); %orbit_height*orbit_height/(sqrt((orbit_height + r_Earth)^2 - (r_Earth * cos(theta(m)))^2 ) - r_Earth * sin(theta(m)))^2;
         %output_intensity(m) = (1 + 1i*q1/wavelength) * conj(1 + 1i*q1/wavelength) * temp;
         
         
         %using Fresnel Kirchhoff formula
         xtemp = linspace(-length_x/2,length_x/2, Nx);
         ytemp = linspace(-width_ztemp/2,width_ztemp/2, Ny);
                
         A1 = (xtemp + tempX) .* (xtemp + tempX);
         B1 = (ytemp + tempZ) .* (ytemp + tempZ);
         C1 = bsxfun(@plus,A1,B1');
         
         R_Station = sqrt(C1 + dPlane2); %distance with ground station
         R_LEO = sqrt(C1 + dLEO2);   %distance with satellite
         
         temp11 = sum(exp(1i * 2 * pi / wavelength * (R_Station + R_LEO))./(R_Station .* R_LEO),'all');       
         tttemp = 1i*temp11/wavelength * length_x * width_ztemp /Nx /Ny * dLEO * exp(-1i * 2 * pi * dLEO /wavelength); %(dLEO - dPlane) * dPlane % * exp(-1i * 2 * pi * (dLEO - dPlane) /wavelength) * exp(-1i * 2 * pi * dPlane /wavelength); %  %      
         tempFK = (1 + tttemp) * conj(1 + tttemp) * pathloss;
         output_FK(m) = tempFK;
         
         
         %useing Fresnel approximation
         %C_F = bsxfun(@plus,A1,B1);
              
         temp11_F = sum(exp(1i * pi / wavelength * C1 /dPlane),'all'); 
         tttemp_F = 1i*temp11_F/wavelength/dPlane  * length_x * width_ztemp /Nx /Ny; % * dPlane; %* exp(1i * 2 * pi /wavelength * dPlane)
          tempFN = (1 + tttemp_F) * conj(1 + tttemp_F) * pathloss;
          output_FN(m)= tempFN;
          
          %using Fraunhofer approximation
           ytemp_F = linspace(-width_ztemp/2,width_ztemp/2, N_F);
         % tempFH1 = 1i / (wavelength * dPlane) * sum(exp(-1i * 2 * pi / wavelength / dPlane * (- tempY) * ytemp_F),'all') * exp(1i * pi / (wavelength * dPlane) * tempY * tempY) * width_ytemp /N_F * length_x;
          tempFH1 = 1i * exp(1i * 2 * pi / wavelength * dPlane) / (wavelength * dPlane) * exp(1i * pi / (wavelength * dPlane) * tempZ * tempZ) * sum(exp(-1i * 2 * pi / (wavelength * dPlane) * (- tempZ) * ytemp_F),'all') * width_ztemp /N_F * length_x;
          output_FH(m) = (1 + tempFH1) * conj(1 + tempFH1) * pathloss;
         
                  
         %Formula: Fresnel approximation
         tt1 = sqrt(2/wavelength/dPlane);
         F1 = fresnelc(tt1 * (length_x / 2 - tempX)) + fresnelc(tt1 * (length_x / 2 + tempX)) + 1i *  fresnels(tt1 * (length_x / 2 - tempX)) + 1i *  fresnels(tt1 * (length_x / 2 + tempX));
         F2 = fresnelc(tt1 * (width_ztemp / 2 - tempZ)) + fresnelc(tt1 * (width_ztemp / 2 + tempZ)) + 1i * fresnels(tt1 * (width_ztemp / 2 - tempZ)) + 1i * fresnels(tt1 * (width_ztemp / 2 + tempZ));
         
         tempFN_F = (1 + 1i / 2 * F1 * F2) * conj(1 + 1i / 2 * F1 * F2) * pathloss;
       
         output_FN_F(m) = tempFN_F;
         
         %using Fraunhofer approximation
         lambda = wavelength;
         b = length_x;
         a = width_ztemp;
         z = dPlane;
         
         beta= 0; %b*X/(lambda*z); 
         gamma= a*(-tempZ)/(lambda*z); %a*Y/(lambda*z); 
         temp1 = a*a*b*b/(lambda*z*lambda*z);
         temp2 =  2*a*b/(lambda *z);
         output_FH_F(m) = 1 + temp1*(sinc(beta)*sinc(gamma)).^2  - temp2 * (sinc(beta)) * (sinc(gamma)) * sin((pi*((tempZ * tempZ + 2*z*z)))/(lambda *z));
         tempFH = output_FH_F(m) * pathloss;
         output_FH_F(m) = tempFH; %considering free space path loss
         
        % temp = orbit_height*orbit_height/(sqrt((orbit_height + r_Earth)^2 - (r_Earth * cos(theta(m)))^2 ) - r_Earth * sin(theta(m)))^2; %deal with distance difference.
         %output_intensity(m) = mean(output, 'all') * temp;
     end
     
 %    intensity_all_output1(m_S, m_Y,:) = output_intensity;
          
    % output_intensity1 = output_intensity/mean(output_intensity);
    % output_dB = 10*log10(output_intensity1);
    % output_detrend = detrend(output_dB);  %remove the linear component of the signal intesnity
          
     % var_theta(m_Y,i_alpha_ini) = std(output_detrend)
  %figure
  
  x_time = [start_time : time_step : end_time];
 % y_intensity = 10*log10(output_intensity);
  y_FK = 10*log10(output_FK);
 % y_NUM_temp = 10*log10(output_NUM_temp);
  y_FH_F = 10*log10(output_FH_F);
  
  y_FH = 10*log10(output_FH);
  
  
  y_FN_F = 10*log10(output_FN_F);
  y_FN = 10*log10(output_FN);
  
 %y_intensity = output_intensity;
 
% spectrogram(y_intensity,128,20,128,1/time_step);
% plot(x_time, y_FK, 'b');
% legend('Fresnel Kirchhoff');

%{
lot(x_time,y_FH,'c', x_time,y_FN_F,'r',x_time,y_FN,'x' );
legend('Fraunhofer', 'Fresnel', 'Fresnel Numerical');
%}

%y_intensity = output_FN;
 %time_step = x_time(2) - x_time(1);
 % samples = 1;
 %y_intensity_update = y_intensity(1:samples:length(y_intensity));
 %spectrogram(y_intensity_update,128,127,512,1/(time_step*samples));
 
 figure
plot(x_time, y_FK, '-k',x_time,y_FH_F,'--', x_time, y_FH,'-.o',x_time,y_FN_F,'-.r', x_time,y_FN,'-xy');
legend('Fresnel Kirchhoff','Fraunhofer', 'Frauhofer Simulation','Fresnel', 'Fresnel Simulation');
 
title(strcat(strcat(num2str(band(wave_select)),strcat(strcat(" GHz    ", num2str(round(theta_initial_degree,2))), ' degree')), strcat("   ", strcat(num2str(plane_speed),' m/sec')))); 
 
 save(savefilename,'x_time','y_FK','y_FH_F','y_FH','y_FN_F','y_FN','doppler_slope');
  %alpha_initial = alpha_initial - alpha_between_satellite;
%end

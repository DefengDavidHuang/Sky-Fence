clear all

SNR = 20; %10; %0; %20; %-2000; %-200; %-100; %-40; %20; %-20; %20; %in dB
data_length = 1000; %100; % %100; %1e3; %1e4; %100; %1e3; 

%modulation scheme
%c = [-5 -5i 5 5i -3 -3-3i -3i 3-3i 3 3+3i 3i -3+3i -1 -1i 1 1i]; % 16-QAM constellation based on the V.29 standard for telephone-line modems.
c = [-3-3i -3-1i -3+1i -3+3i -1-3i -1-1i -1+1i -1+3i 1-3i 1-1i 1+1i 1+3i 3-3i 3-1i 3+1i 3+3i];  %standard 16QAM
M = length(c);

plane_height = 1.0e4;  %1.0e3; %
orbit_height = 5.5e5; %5.5e7; %5.5e6; % %3.0e5;
r_Earth = 6.371e6;

T_satellite = 2*pi * sqrt((orbit_height + r_Earth)^3/((5.9722*6.6743) * 1.0e13));

plane_speed = 250; %-250; % %-125% -250; %-600; %-50; %600; %250; %0; % %-600; %-50; %  -250; % %600; % %600; % % %600; %200;   0; %50; % %250m/second = 900km/hour

time_step = 0.0051; %0.0021; %0.00021; %0.00051; %0.000051; %05;
total_time = 3; %0.3; %3; % 0.3; %3; %0.3; %5e-3; % %0.1; %0.6; %0.6; %0.0006; % % second  
start_time = -1.5; %-0.15; %-1.5; %-0.15; %-1.5; %-0.15; %-0.05; % %-0.05; %0.3; %-0.3; %-0.0003; % %
end_time = start_time + total_time;

Nx = 401; %0; %00; %1000;
Ny = 401; %0; %00; %1000;
N_F = Nx * Ny + 1;
%%%%%%%%%

%antenna size 0.4m x 0.4m
%antenna_halfwidth = 0.22; %0.2;
%antenna_halflength = 0.22; % 0.2;

wave_select = 1; %2; %3; %
wave_set = [1.0e-2, 2.0e-2, 0.5e-2];
band = [30, 15, 60];
wavelength =  wave_set(wave_select); %1.0e-2; %30GHZ 0.5e-2; %60GHz %2.0e-2; %15GHZ %

width_z = 40; %40; %20; %20; %20; %10; %4; %2; %0; %0; %10; %20; %20; %20.0;
length_x = 20; %20; %10; %10; %10; %4; % 2; %0; %10; %0; %20; %20.0;    %x: parallel to ground and perpendicular to the satellite link

Num_data_point = length([start_time : time_step : end_time]); %round(total_time/time_step,0); 

%alpha_initial = 0.089999999999999999999999999999999999; %elevation angle: 40 degrees
%alpha_initial = 0; %elevation angle: 90 degrees
alpha_initial = 0.03674999999999; %elevation angle: 65 degrees

%[1/60:-1/120:-1/60]*pi; %[1/30:-1/120:-1/30]*pi; % [1/60:-1/2400:-1/60]*pi; %  %[1/300:-1/600:-1/300]*pi
theta_ini = pi/2 - atan(((orbit_height + r_Earth)*sin(alpha_initial))/(((orbit_height+r_Earth)*cos(alpha_initial)) - r_Earth));
theta_initial_degree = theta_ini/pi * 180

%Note: the speed of the airplane is in the x direction (for ground static coordinate)

% Plane parameters in ground static coordinate
v_xG = 0.0;
v_yG =  - plane_speed;
v_zG = 0.0;

%x_G = 0;
%y_G = plane_height / tan(theta_ini);
%z_G = plane_height;

% Plane parameters in satellite_ground link fixed coordinate
%x_F = x_G;
%temp = [cos(theta_ini) sin(theta_ini); -sin(theta_ini) cos(theta_ini)] * [y_G z_G]';
%y_F = temp(1);
%z_F = temp(2);
%d_Sat = sqrt((orbit_height + r_Earth)^2 - (r_Earth * cos(theta_ini))^2 ) - r_Earth * sin(theta_ini);

output_intensity = zeros(Num_data_point,1);
output_FK = output_intensity;
%output_FN = output_intensity;
%output_FN_F = output_intensity;
%output_FH_F = output_intensity;
%output_FH = output_intensity;
output_FKMPSK = output_intensity;
% distanceY = ([-total_time/2 : time_step : total_time/2]) * plane_speed; % - 1000;

%plane location along with time
distanceY1st = -1.0 * ([start_time : time_step : end_time]) * plane_speed + plane_height/tan(theta_ini); % relative to the reference ground station (y=0), the minus sign due to the distance at the beginning is set to be positive          
distanceY = distanceY1st; % - Gstation_location(m_Y); %distance relative to the line perpendicular to the ground station of interest

%satellite elevation angle
alpha = -1.0 * [start_time : time_step : end_time] * 2 * pi/T_satellite + alpha_initial; % - Gstation_location(m_Y)/r_Earth;
theta = pi/2 - atan(((orbit_height+r_Earth)*sin(alpha))./(((orbit_height+r_Earth)*cos(alpha)) - r_Earth)); %the elevation angle of satellite relative to the ground station of interest

%satellite elevation angle change rate
temp = (orbit_height +  r_Earth);     
Vdthetadt =  (temp * r_Earth * cos(alpha) - temp * temp)./(temp * temp - 2 * temp * r_Earth * cos(alpha) + r_Earth * r_Earth)  * (- 2 * pi / T_satellite);
    
%coordinate conversion
Vx_G = zeros(size(theta));
Vy_G =  distanceY;                              % plane_height ./ tan(theta);
Vz_G = plane_height * ones(size(theta));

     % Plane parameters in satellite_ground link fixed coordinate
Vx_F = Vx_G;   
Vy_F = cos(theta) .* Vy_G + sin(theta) .* Vz_G; % - y_F;
Vz_F = -sin(theta) .* Vy_G + cos(theta) .* Vz_G; % - z_F;

    %speed conversion
v_x = v_xG;   
Vv_y = cos(theta) * v_yG + sin(theta) * v_zG + ( -sin(theta) .* Vy_G + cos(theta) .* Vz_G) .* Vdthetadt;
Vv_z = -sin(theta) * v_yG + cos(theta) * v_zG + ( -cos(theta) .* Vy_G - sin(theta) .* Vz_G) .* Vdthetadt;

r1 = sqrt(plane_height * plane_height + distanceY .* distanceY);
%r1 = sqrt(plane_height * plane_height + plane_height*plane_height/tan(theta_ini)/tan(theta_ini));

Vy_S1 = sqrt((orbit_height + r_Earth)^2 - (r_Earth * cos(theta)).^2 ) - r_Earth * sin(theta);

%r2 = sqrt((orbit_height + r_Earth)^2 - (r_Earth * cos(theta_ini)).^2 ) - r_Earth * sin(theta_ini) - r1;
r2 =  Vy_S1 - r1;
Vy_S = Vy_S1; % - d_Sat;
    
Vv_yS = - cos(theta + alpha) * 2 * pi * (r_Earth + orbit_height) / T_satellite;

%{
Vdoppler1 = 1/wavelength * (x_F * v_x + y_F * Vv_y + z_F * Vv_z) / r1;
Vdoppler2 = 1/wavelength * ((x_F * v_x + (y_F - d_Sat) * (Vv_y - Vv_yS) + z_F * Vv_z) / r2 - Vv_yS);
Vdoppler3 = 1/wavelength * (Vx_F * v_x + Vy_F .* Vv_y + Vz_F .* Vv_z) / r1;
Vdoppler4 = 1/wavelength * (Vx_F * v_x + (Vy_F - Vy_S) .* (Vv_y - Vv_yS) + Vz_F .* Vv_z) / r2;
%}

%vector mode
%Vdoppler1 = 1/wavelength * (x_F * v_x + y_F * Vv_y + z_F * Vv_z) ./ r1;
%Vdoppler2 = 1/wavelength * ((x_F * v_x + (y_F - d_Sat) * (Vv_y - Vv_yS) + z_F * Vv_z) ./ r2 - Vv_yS);
Vdoppler3 = 1/wavelength * (Vx_F * v_x + Vy_F .* Vv_y + Vz_F .* Vv_z) ./ r1;
Vdoppler4 = 1/wavelength * ((Vx_F * v_x + (Vy_F - Vy_S) .* (Vv_y - Vv_yS) + Vz_F .* Vv_z) ./ r2 - Vv_yS);

%Vdoppler = Vdoppler1 + Vdoppler2 + Vdoppler3 + Vdoppler4;

Vdoppler = Vdoppler3 + Vdoppler4;
Vwidth_z = width_z * sin(theta);
Vbound = (1 ./ r1 + 1 ./ r2)/wavelength .* (abs(length_x * v_x) + abs(Vwidth_z .* Vv_z))/2;

test1 = strcat(num2str(SNR),'dB');
test2 = strcat(num2str(round(theta_initial_degree,2)), 'degree');
test3 = strcat(test1,test2);
test4 = strcat(num2str(plane_speed),'mps');
test5 = strcat(test3,test4);
%total_time1 = round(total_time*10);
test51 = strcat(test5,num2str(data_length));
savefilename = strcat(test51,'SameOrbitAmplitudeQAM.mat');
     
for m = 1:1:Num_data_point
         %m
         width_ztemp = Vwidth_z(m);
         %distancez = plane_height/sin(theta(m)) + cos(theta(m)) * (distanceY(m)  - plane_height / tan(theta(m))); %+ plane_height /tan(theta_ini));
         
         dPlane = plane_height/sin(theta(m)) + (distanceY(m)  - plane_height / tan(theta(m))) * cos(theta(m)); %projected onto the optical axis, distance between flying object and groud receiver 
         dLEO = sqrt((orbit_height + r_Earth)^2 - (r_Earth * cos(theta(m)))^2 ) - r_Earth * sin(theta(m));     % distance between satellite and ground station
         
         dPlane2 = dPlane * dPlane;
         dLEO2 = (dLEO - dPlane) * (dLEO - dPlane);
                
        % funFresnel = @(r1,r2) A*exp(1i * 2 * pi / wavelength * (r1 + r2))./(r1 .* r2) ;
         
         tempY = sin(theta(m)) * (distanceY(m)  - plane_height / tan(theta(m)));   %deviation of flying object center from the optical axis
         %fun1 = @(x1,z1) funFresnel(sqrt(x1.*x1 + (z1 + tempY).*(z1 + tempY) + dPlane2), sqrt(x1 .* x1 + (z1 + tempY) .* (z1 + tempY) + dLEO2));
         
         %q1 = quad2d(fun1,-length_x/2,length_x/2,-width_ytemp/2,width_ytemp/2); %, 'RelTol', 0.1); %'AbsTol',1); %1e-1);
         
         %q1 =integral2(fun1,-length_x/2,length_x/2,-width_ytemp/2,width_ytemp/2);
        
         pathloss_am = orbit_height/dLEO;
         pathloss = pathloss_am * pathloss_am; %orbit_height*orbit_height/(sqrt((orbit_height + r_Earth)^2 - (r_Earth * cos(theta(m)))^2 ) - r_Earth * sin(theta(m)))^2;
         %output_intensity(m) = (1 + 1i*q1/wavelength) * conj(1 + 1i*q1/wavelength) * temp;
              
         %using Fresnel Kirchhoff formula
         xtemp = linspace(-length_x/2,length_x/2, Nx);
         ytemp = linspace(-width_ztemp/2,width_ztemp/2, Ny);
                
         A1 = xtemp .* xtemp;
         B1 = (ytemp - tempY) .* (ytemp - tempY);
         C1 = bsxfun(@plus,A1,B1');
         
         R_Station = sqrt(C1 + dPlane2); %distance with ground station
         R_LEO = sqrt(C1 + dLEO2);   %distance with satellite
         
         temp11 = sum(exp(1i * 2 * pi / wavelength * (R_Station + R_LEO))./(R_Station .* R_LEO),'all');       
         tttemp = 1i*temp11/wavelength * length_x * width_ztemp /Nx /Ny * dLEO * exp(-1i * 2 * pi * dLEO /wavelength); %(dLEO - dPlane) * dPlane % * exp(-1i * 2 * pi * (dLEO - dPlane) /wavelength) * exp(-1i * 2 * pi * dPlane /wavelength); %  %      
         RXgain = (1 + tttemp) * sqrt(pathloss);
         tempFK = RXgain * conj(RXgain); % * pathloss;
         %tempFK_am = (1 + tttemp) * pathloss_am;
         output_FK(m) = sqrt(tempFK);
         
         %here about digital communications???
         data = randi([0 M-1],data_length,1);
         modData = genqammod(data,c);
         rxSig = awgn(modData*RXgain,SNR,'measured','db');
         output_FKMPSK(m) = mean(abs(rxSig));
         % = tempFK;
         % obtain the amplitudes and averaging out
     end
     
 %    intensity_all_output1(m_S, m_Y,:) = output_intensity;
          
    % output_intensity1 = output_intensity/mean(output_intensity);
    % output_dB = 10*log10(output_intensity1);
    % output_detrend = detrend(output_dB);  %remove the linear component of the signal intesnity
          
     % var_theta(m_Y,i_alpha_ini) = std(output_detrend)
  figure
  
  x_time = [start_time : time_step : end_time];
  %{
 % y_intensity = 10*log10(output_intensity);
  y_FK = 10*log10(output_FK);
 % y_NUM_temp = 10*log10(output_NUM_temp);
  y_FH_F = 10*log10(output_FH_F);
  
  y_FH = 10*log10(output_FH);
  
  
  y_FN_F = 10*log10(output_FN_F);
  y_FN = 10*log10(output_FN);
  %}
  
 %y_intensity = output_intensity;
 
% spectrogram(y_intensity,128,20,128,1/time_step);
% plot(x_time, y_FK, 'b');
% legend('Fresnel Kirchhoff');

%{
lot(x_time,y_FH,'c', x_time,y_FN_F,'r',x_time,y_FN,'x' );
legend('Fraunhofer', 'Fresnel', 'Fresnel Numerical');
%}
%{
plot(x_time, y_FK, '-k',x_time,y_FH_F,'--', x_time, y_FH,'-.o',x_time,y_FN_F,'-.r', x_time,y_FN,'-xy');
legend('Fresnel Kirchhoff','Fraunhofer', 'Frauhofer Simulation','Fresnel', 'Fresnel Simulation');
%}

plot(x_time, abs(Vdoppler),'-k', x_time, abs(Vdoppler + Vbound), '--k', x_time, abs(Vdoppler - Vbound),'-.k');
legend('Median','Upper Bound','Lower Bound');
%title(strcat(strcat(num2str(band(wave_select)),strcat(strcat(" GHz    ", num2str(round(theta_initial_degree,2))), ' degree')), strcat("   ", strcat(num2str(plane_speed),' m/sec')))); 
save(savefilename,'width_z','length_x','x_time','output_FK','output_FKMPSK','Vdoppler','Vbound');
%save(savefilename,'width_z','length_x','x_time','output_FK','output_FH_F','output_FH','output_FN','Vdoppler','Vbound');
%save(savefilename,'x_time','y_FK','y_FH_F','y_FH','y_FN_F','y_FN','Vdoppler'); %, 'doppler_slope','doppler_offset');
 %save('FK90degree2010update.mat','x_time','y_FK','y_FH_F','y_FH','y_FN_F','y_FN');
  %alpha_initial = alpha_initial - alpha_between_satellite;
%end

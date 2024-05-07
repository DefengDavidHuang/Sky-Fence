clear all
format long
%This code generates the signals in H0 (no target) and H1(with target)

plane_height = 500.1; %101.1; %500; %1.0e4;  %1.0e3; %
plane_speed = 100.1; %300.1; %100.1; %10; %30; %plane speed in m/s

width_z = 6.4; %1; %2.1; %6.1; %length of plane
length_x = 8.1; %3.1; %8.1; %width of plane  %x: parallel to ground and perpendicular to the satellite link

% offset from x=0
x_track_initial = 20.1; %0.0; %100.1; %200.1; %0.0; %150.1; %111.1; %5; % 101.1; %in the x direction, the offset of the center of the airplane to the plane scanned by the satellite-reciever link 

%setup the elevation angle

%alpha_initial = 0.089999999999999999999999999999999999; %elevation angle: 40 degrees
%alpha_initial = 0; %elevation angle: 90 degrees
alpha_initial = 0.03674999999999; %elevation angle: 65 degrees

time_step = 5.0e-8; %5.01e-7; %0.000501; %0.000501; %0.000501; %0.0021; %0.00021; %0.00051; %0.000051; %05;
total_time = 0.4; %0.3; %0.3; %2; %5; %2; %0.6; %0.6; %1; %8; %10.0; %5.0; %1.0; %0.1; %0.25; %0.5; %2; %5; %1; %3; %0.3; %3; % 0.3; %3; %0.3; %5e-3; % %0.1; %0.6; %0.6; %0.0006; % % second  
start_time = -total_time/2; %-1; -2.5; %-0.5; %-1.5; %-0.15; %-1.5; %-0.15; %-1.5; %-0.15; %-0.05; % %-0.05; %0.3; %-0.3; %-0.0003; % %
end_time = start_time + total_time;

orbit_height = 5.5e5; %5.5e7; %5.5e6; % %3.0e5;
r_Earth = 6.371e6;

T_satellite = 2*pi * sqrt((orbit_height + r_Earth)^3/((5.9722*6.6743) * 1.0e13));

%number of samples used for diffraction modelling of target
Nx = 401; %0; %00; %1000;
Nz = 401; %0; %00; %1000;

%antenna size 0.4m x 0.4m
%antenna_halfwidth = 0.22; %0.2;
%antenna_halflength = 0.22; % 0.2;

wave_select = 1; %2; %3; %
wave_set = [1.0e-2, 2.0e-2, 0.5e-2];
band = [30, 15, 60];
wavelength =  wave_set(wave_select); %1.0e-2; %30GHZ 0.5e-2; %60GHz %2.0e-2; %15GHZ %

Num_data_point = length([start_time : time_step : end_time]); %round(total_time/time_step,0); 

theta_ini = pi/2 - atan(((orbit_height + r_Earth)*sin(alpha_initial))/(((orbit_height+r_Earth)*cos(alpha_initial)) - r_Earth));
theta_initial_degree = theta_ini/pi * 180

%Note: the speed of the airplane is in the yG direction (i.e., parallel with the ground in the ground static coordinate)

%vectors for output
output_intensity = zeros(Num_data_point,1);
output_FK = output_intensity;
output_N0 = output_intensity;
outputgain_N1 = output_intensity;
outputgain_N0 = output_intensity;
outputgain_target = output_intensity;

%plane location along with time
distanceY = -1.0 * ([start_time : time_step : end_time]) * plane_speed + plane_height/tan(theta_ini); % relative to the reference ground station (y=0), the minus sign due to the distance at the beginning is set to be positive          
%distanceY = distanceY1st; % - Gstation_location(m_Y); %distance relative to the line perpendicular to the ground station of interest

%satellite elevation angle 
alpha = -1.0 * [start_time : time_step : end_time] * 2 * pi/T_satellite + alpha_initial; % - Gstation_location(m_Y)/r_Earth;
theta = pi/2 - atan(((orbit_height+r_Earth)*sin(alpha))./(((orbit_height+r_Earth)*cos(alpha)) - r_Earth)); %the elevation angle of satellite relative to the ground station of interest for each sampling time

%satellite elevation angle change rate
temp = (orbit_height +  r_Earth);     
Vdthetadt =  (temp * r_Earth * cos(alpha) - temp * temp)./(temp * temp - 2 * temp * r_Earth * cos(alpha) + r_Earth * r_Earth)  * (- 2 * pi / T_satellite);


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

dthetadt =  (temp * r_Earth * cos(alpha_initial) - temp * temp)/(temp * temp - 2 * temp * r_Earth * cos(alpha_initial) + r_Earth * r_Earth)  * (- 2 * pi / T_satellite);

    
%coordinate conversion
%coordinate attached with the ground
Vx_G = zeros(size(theta)) + x_track_initial;
Vy_G =  distanceY;                              % plane_height ./ tan(theta);
Vz_G = plane_height * ones(size(theta));

% Plane parameters in coordinate attached to satellite_ground link
Vx_F = Vx_G;   
Vy_F = cos(theta) .* Vy_G + sin(theta) .* Vz_G; % - y_F;
Vz_F = -sin(theta) .* Vy_G + cos(theta) .* Vz_G; % - z_F;

% Speed of target in ground static coordinate
v_xG = 0.0;
v_yG =  - plane_speed;
v_zG = 0.0;

%speed in coordinate attached to satellite-ground station link
v_x = v_xG;   
Vv_y = cos(theta) * v_yG + sin(theta) * v_zG + ( -sin(theta) .* Vy_G + cos(theta) .* Vz_G) .* Vdthetadt;
Vv_z = -sin(theta) * v_yG + cos(theta) * v_zG + ( -cos(theta) .* Vy_G - sin(theta) .* Vz_G) .* Vdthetadt;

% distance related to the center of the target
r1 = sqrt(Vx_G .* Vx_G + Vy_G .* Vy_G + Vz_G .* Vz_G);  %sqrt(plane_height * plane_height + distanceY .* distanceY);

Vy_S = sqrt((orbit_height + r_Earth)^2 - (r_Earth * cos(theta)).^2 ) - r_Earth * sin(theta);

r2 = sqrt(Vx_F .* Vx_F + (Vy_S - Vy_F) .* (Vy_S - Vy_F) + Vz_F .* Vz_F);

%Vy_S = Vy_S1; % - d_Sat;
    
Vv_yS = - cos(theta + alpha) * 2 * pi * (r_Earth + orbit_height) / T_satellite;

Vdoppler3 = 1/wavelength * (Vx_F * v_x + Vy_F .* Vv_y + Vz_F .* Vv_z) ./ r1;
Vdoppler4 = 1/wavelength * ((Vx_F * v_x + (Vy_F - Vy_S) .* (Vv_y - Vv_yS) + Vz_F .* Vv_z) ./ r2 - Vv_yS);
Vdoppler = Vdoppler3 + Vdoppler4;

Vwidth_z = width_z * sin(theta);
Vbound = (1 ./ r1 + 1 ./ r2)/wavelength .* (abs(length_x * v_x) + abs(Vwidth_z .* Vv_z))/2;

%Note: the speed of the airplane is in the x direction (for ground static
%coordinate)
vsquare = plane_speed * plane_speed + (plane_height/sin(theta_ini) * dthetadt)^2; 
doppler_slope = vsquare/wavelength/(plane_height/sin(theta_ini))

test1 = strcat(num2str(x_track_initial),'offsetinX');
test2 = strcat(num2str(round(theta_initial_degree,2)), 'degree');
test3 = strcat(test1,test2);
test4 = strcat(num2str(plane_speed),'mps');
test5 = strcat(test3,test4);
savefilename = strcat(test5,'signal02secondlower20M.mat');

x_time = [start_time : time_step : end_time];
   
for m = 1:1:Num_data_point
         %m
         width_ztemp = Vwidth_z(m);
         %distancez = plane_height/sin(theta(m)) + cos(theta(m)) * (distanceY(m)  - plane_height / tan(theta(m))); %+ plane_height /tan(theta_ini));
         
                %dPlane = plane_height/sin(theta(m)) + (distanceY(m)  - plane_height / tan(theta(m))) * cos(theta(m)) %projected onto the optical axis, distance between flying object and groud receiver 
         dPlane = Vy_F(m);
         
         %dLEO = sqrt((orbit_height + r_Earth)^2 - (r_Earth * cos(theta(m)))^2 ) - r_Earth * sin(theta(m));     % distance between satellite and ground station
         dLEO = Vy_S(m);       
         dPlane2 = dPlane * dPlane;
         dLEO2 = (dLEO - dPlane) * (dLEO - dPlane);
                
        % funFresnel = @(r1,r2) A*exp(1i * 2 * pi / wavelength * (r1 + r2))./(r1 .* r2) ;
         
         tempZ = sin(theta(m)) * (distanceY(m)  - plane_height / tan(theta(m)));   %deviation of flying object center from the optical axis
         %fun1 = @(x1,z1) funFresnel(sqrt(x1.*x1 + (z1 + tempY).*(z1 + tempY) + dPlane2), sqrt(x1 .* x1 + (z1 + tempY) .* (z1 + tempY) + dLEO2));
         %q1 = quad2d(fun1,-length_x/2,length_x/2,-width_ytemp/2,width_ytemp/2); %, 'RelTol', 0.1); %'AbsTol',1); %1e-1);
                  %q1 =integral2(fun1,-length_x/2,length_x/2,-width_ytemp/2,width_ytemp/2);
        
         pathloss_am = orbit_height/dLEO; % normalized by LEO overhead
         pathloss = pathloss_am * pathloss_am; %orbit_height*orbit_height/(sqrt((orbit_height + r_Earth)^2 - (r_Earth * cos(theta(m)))^2 ) - r_Earth * sin(theta(m)))^2;
         %output_intensity(m) = (1 + 1i*q1/wavelength) * conj(1 + 1i*q1/wavelength) * temp;
              
         %using Fresnel Kirchhoff formula
         xtemp = linspace(-length_x/2,length_x/2, Nx);
         ztemp = linspace(-width_ztemp/2,width_ztemp/2, Nz);
                
         A1 = (xtemp + x_track_initial) .* (xtemp + x_track_initial); %square of the x-coordinate
         B1 = (ztemp - tempZ) .* (ztemp - tempZ); %square of the z-coordinate
         C1 = bsxfun(@plus,A1,B1');
         
         R_Station = sqrt(C1 + dPlane2); %distance with ground station (i.e. r1)
         R_LEO = sqrt(C1 + dLEO2);   %distance with satellite (i.e. r2)
         
         temp11 = sum(exp(1i * 2 * pi / wavelength * (R_Station + R_LEO))./(R_Station .* R_LEO),'all');       
         tttemp = 1i*temp11/wavelength * length_x * width_ztemp /Nx /Nz * dLEO * exp(-1i * 2 * pi * dLEO /wavelength); %(dLEO - dPlane) * dPlane % * exp(-1i * 2 * pi * (dLEO - dPlane) /wavelength) * exp(-1i * 2 * pi * dPlane /wavelength); %  %      
         RXgain = (1 + tttemp)/dLEO * orbit_height * exp(1i * 2 * pi * dLEO /wavelength); %/dLEO * orbit_height; %normalized by LEO overhead
         outputgain_N1(m) = RXgain;  %with flying object, total signal (complex-valued)
         
         outputgain_target(m) = tttemp/dLEO * orbit_height;   %gain due to target only
         outputgain_N0(m) = exp(1i * 2 * pi * dLEO /wavelength)/dLEO * orbit_height; %no flying object, normalized by LEO overhead sqrt(pathloss);
         
         tempFK = RXgain * conj(RXgain); % * pathloss;
         %tempFK_am = (1 + tttemp) * pathloss_am;
         output_FK(m) = sqrt(tempFK);  %with flying oject, amplitude of total signal
         output_N0(m) = sqrt(outputgain_N0(m) * conj(outputgain_N0(m))); %without flying object
         if mod(m,100000) == 0
             m
              y_FK = output_FK;
             save(savefilename,'x_time','y_FK','output_N0','Vdoppler','Vbound','outputgain_N1','outputgain_target','outputgain_N0', 'doppler_slope');
         end
     end
     
  %figure
  

  
  y_FK = output_FK;
 % y_intensity = 10*log10(output_intensity);
 %y_FK = output_FKMPSK; %
 % y_FK = 20*log10(output_FKMPSK);
  %y_FK = 20*log10(output_FK);
  
  %save(savefilename,'x_time','y_FK','output_N0','outputgain_N1','outputgain_target','outputgain_N0');
  save(savefilename,'x_time','y_FK','output_N0','Vdoppler','Vbound','outputgain_N1','outputgain_target','outputgain_N0', 'doppler_slope');
  %plot(x_time, abs(outputgain_target), 'b');
  %figure
  
  %plot(x_time, y_FK, 'b',x_time,output_N0,'.k');
  %output_target = sqrt(outputgain_target.*conj(outputgain_target));
  
  y_FK_detrend = detrend(y_FK,'linear');
  output_N0detrend = detrend(output_N0,'linear');
  %figure
  plot(x_time, y_FK, 'b',x_time,output_N0,'.k');
legend('Fresnel Kirchhoff','channel gain in power');
xlabel('time (second)');
ylabel('received signal (linear)');

figure
plot(x_time,abs(Vdoppler + Vbound),'-.k', x_time, abs(Vdoppler - Vbound), '--k');
legend('bound 1', 'bound 2');
%plot(x_time,Vdoppler, 'b');
%ylabel('received signal (dB)');
 % y_NUM_temp = 10*log10(output_NUM_temp);
  %y_FH_F = 10*log10(output_FH_F);
  
  %y_FH = 10*log10(output_FH);
  
  
  %y_FN_F = 10*log10(output_FN_F);
  %y_FN = 10*log10(output_FN);
  
  
 %y_intensity = output_intensity;
 
% spectrogram(y_intensity,128,20,128,1/time_step);


%{
lot(x_time,y_FH,'c', x_time,y_FN_F,'r',x_time,y_FN,'x' );
legend('Fraunhofer', 'Fresnel', 'Fresnel Numerical');
%}
%{
plot(x_time, y_FK, '-k',x_time,y_FH_F,'--', x_time, y_FH,'-.o',x_time,y_FN_F,'-.r', x_time,y_FN,'-xy');
legend('Fresnel Kirchhoff','Fraunhofer', 'Frauhofer Simulation','Fresnel', 'Fresnel Simulation');
%}

%plot(x_time, abs(Vdoppler),'-k', x_time, abs(Vdoppler + Vbound), '--k', x_time, abs(Vdoppler - Vbound),'-.k');
%legend('Median','Upper Bound','Lower Bound');
%title(strcat(strcat(num2str(band(wave_select)),strcat(strcat(" GHz    ", num2str(round(theta_initial_degree,2))), ' degree')), strcat("   ", strcat(num2str(plane_speed),' m/sec')))); 
%save(savefilename,'width_z','length_x','x_time','output_FK','output_FKMPSK','Vdoppler','Vbound');
%save(savefilename,'width_z','length_x','x_time','output_FK','output_FH_F','output_FH','output_FN','Vdoppler','Vbound');
%save(savefilename,'x_time','y_FK','y_FH_F','y_FH','y_FN_F','y_FN','Vdoppler'); %, 'doppler_slope','doppler_offset');
 %save('FK90degree2010update.mat','x_time','y_FK','y_FH_F','y_FH','y_FN_F','y_FN');
  %alpha_initial = alpha_initial - alpha_between_satellite;
%end

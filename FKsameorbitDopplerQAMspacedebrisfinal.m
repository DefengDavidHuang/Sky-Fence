clear all

SNR = 10; %20; %-200; %-100; %-40; %20; %-20; %20; %in dB

% Assume 100MSymbols/second rate, then a data_length of 1e3 can measure
% Doppler frequency up to 50kHz (1 Hz needs 2 samples)
data_length = 1e3; %1e4; %1e4; %100; %1e3; 

%modulation scheme
%c = [-5 -5i 5 5i -3 -3-3i -3i 3-3i 3 3+3i 3i -3+3i -1 -1i 1 1i]; % 16-QAM constellation based on the V.29 standard for telephone-line modems.
c = [-3-3i -3-1i -3+1i -3+3i -1-3i -1-1i -1+1i -1+3i 1-3i 1-1i 1+1i 1+3i 3-3i 3-1i 3+1i 3+3i];  %standard 16QAM
M = length(c);

Debris_height = 5.4e5; %5e5; %5.4e5; %5.4e5; %4e5; %5.4999e5; % % %4.0e5; %4.0e5;  %1.0e3; %
%plane_height = Debris_height;
orbit_height = 5.5e5; %5.5e7; %5.5e6; % %3.0e5;
r_Earth = 6.371e6;

flagdirection = 1; %-1; % 1 for same direction; -1 for opposite direction

T_satellite = 2*pi * sqrt((orbit_height + r_Earth)^3/((5.9722*6.6743) * 1.0e13)); % G = 6.6743e{-11} m^3kg^{-1}s^{-2}, M_{Earth}=5.9722e24 kg
T_Debris = 2*pi * sqrt((Debris_height + r_Earth)^3/((5.9722*6.6743) * 1.0e13));

time_step = 0.00141; %0.00041; %0.00021; %0.000054001; %0.00041; % % %2e-5; %0.00021; %2e-5; %0.021; %0.51e-3; %0.021; %0.99e-4; %0.99e-4 %0.021; %0.0021; %0.00021; %0.00021; %0.00051; %0.000051; %05;
total_time = 2; %0.015; %10; %0.02; %50; %2; %20; %3; %0.3; %3; % 0.3; %3; %0.3; %5e-3; % %0.1; %0.6; %0.6; %0.0006; % % second  
start_time = -1; %-0.0075; %-5; %-0.01; %-25; %3; %-10; %-1.5; %-0.15; %-1.5; %-0.15; %-1.5; %-0.15; %-0.05; % %-0.05; %0.3; %-0.3; %-0.0003; % %
end_time = start_time + total_time;
x_time = [start_time : time_step : end_time];
  
Nx = 401; %0; %00; %1000;
Ny = 401; %0; %00; %1000;
N_F = Nx * Ny + 1;
%%%%%%%%%

wave_select = 1; %2; %3; %
wave_set = [1.0e-2, 2.0e-2, 0.5e-2];
band = [30, 15, 60];
wavelength =  wave_set(wave_select); %1.0e-2; %30GHZ 0.5e-2; %60GHz %2.0e-2; %15GHZ %

width_z = 40; %20; %50; %20; %50; %0.5; %50; %20; %5; %10; %2; %20; %40; %50; %20; %40; %0.4; %40; %20; %20; %20; %10; %4; %2; %0; %0; %10; %20; %20; %20.0;
length_x = 20; %10; %25; %10; %25; %0.2; %25; %10; %2; %5; %1; %10; %20; %25; %10; %20; %0.2; %20; %10; %10; %10; %4; % 2; %0; %10; %0; %20; %20.0;    %x: parallel to ground and perpendicular to the satellite link

Num_data_point = length([start_time : time_step : end_time]); %round(total_time/time_step,0); 

alpha_initial = 0.089999999999999999999999999999999999; %elevation angle: 40 degrees
%alpha_initial = 0; %elevation angle: 90 degrees
%alpha_initial = 0.03674999999999; %elevation angle: 65 degrees

%[1/60:-1/120:-1/60]*pi; %[1/30:-1/120:-1/30]*pi; % [1/60:-1/2400:-1/60]*pi; %  %[1/300:-1/600:-1/300]*pi
theta_ini = pi/2 - atan(((orbit_height + r_Earth)*sin(alpha_initial))/(((orbit_height+r_Earth)*cos(alpha_initial)) - r_Earth)); %the same for satellite and debris
theta_initial_degree = theta_ini/pi * 180

%note the following is for theta_ini less than pi/2 only
%same initial elevation angle for the debris and the LEO
temp = Debris_height +  r_Earth;
alpha_debris_ini = acos((r_Earth * cos(theta_ini) * cos(theta_ini) + sqrt(sin(theta_ini)^2 * (temp^2 - r_Earth^2 * cos(theta_ini)^2 )))/temp);

test1 = strcat(num2str(band(wave_select)),'GHz');
test2 = strcat(num2str(round(theta_initial_degree,2)), 'degree');
test3 = strcat(test1,test2);
test4 = strcat(num2str(Debris_height),'mSNR');
test5 = strcat(test3,test4);
%total_time1 = SNR; round(total_time*10);
test52 = strcat(test5,num2str(SNR));
test53 = strcat(test52,'dB');
test51 = strcat(test53,num2str(flagdirection));
savefilename = strcat(test51,'directionQAMDebris.mat');

%satellite and Debris elevation angles
alpha = -1.0 * [start_time : time_step : end_time] * 2 * pi/T_satellite + alpha_initial; % - Gstation_location(m_Y)/r_Earth;
temp = orbit_height +  r_Earth;
theta = pi/2 - atan(temp*sin(alpha)./(temp*cos(alpha) - r_Earth)); %the elevation angle of satellite relative to the ground station of interest
VdLEO = sqrt(temp^2 - (r_Earth * cos(theta)).^2) - r_Earth * sin(theta);  %distance between LEO and ground station
 
alpha_Debris = -1.0 * flagdirection * [start_time : time_step : end_time] * 2 * pi/T_Debris + alpha_debris_ini;
temp = Debris_height +  r_Earth;
theta_Debris = pi/2 - atan(temp*sin(alpha_Debris)./(temp*cos(alpha_Debris) - r_Earth)); %elevaton angle of Debris
VdDebris = sqrt(temp^2 - (r_Earth * cos(theta_Debris)).^2) - r_Earth * sin(theta_Debris); %distance between Debris and ground station

Vwidth_z = width_z * cos(theta - theta_Debris); %the width of the Debris projected onto the plane that is perpendicular to the LEO-station link

%Doppler shift and Doppler spread

%coordinatse of Debris

Vy_F = VdDebris .* cos(theta - theta_Debris);
Vz_F = - VdDebris .* sin(theta - theta_Debris);
Vx_F = 0.0 * Vy_F; % Debris coordinate in x axis

%Distance of Debris to ground station
%temp = Debris_height +  r_Earth;
%Vd_Debris = sqrt(temp^2 - (r_Earth * cos(theta_Debris)).^2) - r_Earth * sin(theta_Debris);

%coordinate of Satellite 
%temp = orbit_height +  r_Earth;
Vy_S = VdLEO; %sqrt(temp^2 - (r_Earth * cos(theta)).^2) - r_Earth * sin(theta);

r1 = VdDebris;
r2 = sqrt(Vy_S.^2 + VdDebris.^2 - 2 * Vy_S .* VdDebris .* cos(theta - theta_Debris));

%Velocity of satellite
Vv_yS = - cos(theta + alpha) * 2 * pi * (r_Earth + orbit_height) / T_satellite;

%velocity of Debris
temp = orbit_height +  r_Earth;     
Vdthetadt =  (temp * r_Earth * cos(alpha) - temp * temp)./(temp * temp - 2 * temp * r_Earth * cos(alpha) + r_Earth * r_Earth)  * (- 2 * pi / T_satellite);
temp = Debris_height +  r_Earth;
VdthetaDebrisdt = (temp * r_Earth * cos(alpha_Debris) - temp * temp)./(temp * temp - 2 * temp * r_Earth * cos(alpha_Debris) + r_Earth * r_Earth)  * (- 2 * pi /T_Debris) * flagdirection;

%change rate of distance between debris and station
VdDebrisdt = - cos(theta_Debris + alpha_Debris) * 2 * pi * (r_Earth + Debris_height) / T_Debris * flagdirection;

Vv_y = cos(theta - theta_Debris) .* VdDebrisdt - VdDebris .* sin(theta - theta_Debris) .* (Vdthetadt - VdthetaDebrisdt);
Vv_z = -sin(theta - theta_Debris) .* VdDebrisdt - VdDebris .* cos(theta - theta_Debris) .* (Vdthetadt - VdthetaDebrisdt);
Vv_x = 0.0 * Vv_y;

%Doppler

Vdoppler3 = 1/wavelength * (Vx_F .* Vv_x + Vy_F .* Vv_y + Vz_F .* Vv_z) ./ r1;
Vdoppler4 = 1/wavelength * ((Vx_F .* Vv_x + (Vy_F - Vy_S) .* (Vv_y - Vv_yS) + Vz_F .* Vv_z) ./ r2 - Vv_yS);
         
Vdoppler = Vdoppler3 + Vdoppler4;

Vbound = (1./r1 + 1./r2)/wavelength .* (abs(length_x * Vv_x) + abs(Vwidth_z .* Vv_z))/2;

Deviation = VdDebris .* sin(theta - theta_Debris);
Deviation_projection = VdDebris .* cos(theta - theta_Debris);
%plot(x_time, Deviation, '-k');
for m = 1:1:Num_data_point
         %m
         width_ztemp = Vwidth_z(m);
         %distancez = plane_height/sin(theta(m)) + cos(theta(m)) * (distanceY(m)  - plane_height / tan(theta(m))); %+ plane_height /tan(theta_ini));
         
         dPlane = Deviation_projection(m); %VdDebris(m) * cos(theta(m)-theta_Debris(m)); %projected onto the optical axis, distance between Debris and groud receiver 
         dLEO = VdLEO(m);     % distance between satellite and ground station
                 
         dPlane2 = dPlane * dPlane;
         dLEO2 = (dLEO - dPlane) * (dLEO - dPlane);
                
         tempZ =  Deviation(m); %VdDebris(m) * sin(theta(m)-theta_Debris(m));   %deviation of Debris center from the optical axis
               
         pathloss_am = orbit_height/dLEO;
         pathloss = pathloss_am * pathloss_am; %orbit_height*orbit_height/(sqrt((orbit_height + r_Earth)^2 - (r_Earth * cos(theta(m)))^2 ) - r_Earth * sin(theta(m)))^2;
         %output_intensity(m) = (1 + 1i*q1/wavelength) * conj(1 + 1i*q1/wavelength) * temp;
              
         %using Fresnel Kirchhoff formula
         xtemp = linspace(-length_x/2,length_x/2, Nx);
         ztemp = linspace(-width_ztemp/2,width_ztemp/2, Ny);
                
         A1 = xtemp .* xtemp;
         B1 = (ztemp - tempZ) .* (ztemp - tempZ);
         C1 = bsxfun(@plus,A1,B1');
         
         
         R_Station = sqrt(C1 + dPlane2); %distance with ground station
         R_LEO = sqrt(C1 + dLEO2);   %distance with satellite
                 
         %cos_alpha1 = Deviation_projection(m)./R_Station;
         %cos_alpha2 = (dLEO - Deviation_projection(m))./R_LEO;
         
         temp11 = sum(exp(1i * 2 * pi / wavelength * (R_Station + R_LEO))./(R_Station .* R_LEO),'all'); 
         %temp11 = sum(exp(1i * 2 * pi / wavelength * (R_Station + R_LEO))./(R_Station .* R_LEO).* (cos_alpha1 + cos_alpha2)/2,'all');        % 
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

%plot(x_time,(theta - theta_Debris)*180/pi,'-k');
%figure
plot(x_time,output_FK,'-k',x_time,output_FKMPSK,'-.k');
legend('without noise','with noise');
figure
plot(x_time, abs(Vdoppler),'-k', x_time, abs(Vdoppler + Vbound), '--k', x_time, abs(Vdoppler - Vbound),'-.k');
legend('Median','Upper Bound','Lower Bound');

%title(strcat(strcat(num2str(band(wave_select)),strcat(strcat(" GHz    ", num2str(round(theta_initial_degree,2))), ' degree')), strcat("   ", strcat(num2str(plane_speed),' m/sec')))); 
save(savefilename,'width_z','length_x','x_time','output_FK','output_FKMPSK','Vdoppler','Vbound');
%save(savefilename,'width_z','length_x','x_time','output_FK','output_FKMPSK','Vdoppler','Vbound');

%save(savefilename,'width_z','length_x','x_time','output_FK','output_FH_F','output_FH','output_FN','Vdoppler','Vbound');
%save(savefilename,'x_time','y_FK','y_FH_F','y_FH','y_FN_F','y_FN','Vdoppler'); %, 'doppler_slope','doppler_offset');
 %save('FK90degree2010update.mat','x_time','y_FK','y_FH_F','y_FH','y_FN_F','y_FN');
  %alpha_initial = alpha_initial - alpha_between_satellite;
%end

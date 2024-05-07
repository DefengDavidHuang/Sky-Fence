clear all
format long
%This code generates the signals in H0 (no target) and H1(with target)

%SNR = 10; %20; %10; %0; %20; %-2000; %-200; %-100; %-40; %20; %-20; %20; %in dB
%data_length = 10000; %100; % %100; %1e3; %1e4; %100; %1e3; 

%modulation scheme
%c = [-5 -5i 5 5i -3 -3-3i -3i 3-3i 3 3+3i 3i -3+3i -1 -1i 1 1i]; % 16-QAM constellation based on the V.29 standard for telephone-line modems.
%c = [-3-3i -3-1i -3+1i -3+3i -1-3i -1-1i -1+1i -1+3i 1-3i 1-1i 1+1i 1+3i 3-3i 3-1i 3+1i 3+3i]/sqrt(9.75);  %standard 16QAM
%M = length(c);

%plane_height = 500.00001; %101.1; %500; %1.0e4;  %1.0e3; %
%{
plane_height = 150.00001;

%plane_speed = 100.00001; %10; %30; %108km/hour %0; %-250; %0; %150; %250; %-250; % %-125% -250; %-600; %-50; %600; %250; %0; % %-600; %-50; %  -250; % %600; % %600; % % %600; %200;   0; %50; % %250m/second = 900km/hour
plane_speed = 13.0000001;

%LENGTH = 6.15; % object length in metres
%WIDTH = 8.62;  % object width in metres
LENGTH = 0.3; %
WIDTH = 0.35; %
width_z = LENGTH; %6.15; %length %8.1; %40; %40; %20; %20; %20; %10; %4; %2; %0; %0; %10; %20; %20; %20.0;
length_x = WIDTH; %8.62; %width 1.01% 6.1; %20; %20; %10; %10; %10; %4; % 2; %0; %10; %0; %20; %20.0;    %x: parallel to ground and perpendicular to the satellite link

%x_initial_set = 0:10.01:450; %0.1:10:11; %
x_initial_set = 0:5.01:200.1; %
%}

plane_height = 500; %100; %500; 100; %500.1; %101.1; %500; %1.0e4;  %1.0e3; %
plane_speed = -20; %20; %100.1; %300.1; %100.1; %10; %30; %plane speed in m/s

width_z = 0.64; %0.13; % %0.33; %*2; %0.17; %0.33; %6.4; %1; %2.1; %6.1; %length of plane
length_x =0.81; %0.16; %  %0.45; %*2; % 0.25; %0.45; % 8.1; %3.1; %8.1; %width of plane  %x: parallel to ground and perpendicular to the satellite link

%x_initial_set = 0:2:100; %0:0.5:10; %0:1:250; %0:5.01:2000.1;
x_initial_set = 0:2:100; %0:2:20; %100
% offset from x=0
%x_track_initial = 20.1; %0.0; %100.1; %200.1; %0.0; %150.1; %111.1; %5; % 101.1; %in the x direction, the offset of the center of the airplane to the plane scanned by the satellite-reciever link 
orbit_height = 5.5e5; %5.5e7; %5.5e6; % %3.0e5;
r_Earth = 6.371e6;

T_satellite = 2*pi * sqrt((orbit_height + r_Earth)^3/((5.9722*6.6743) * 1.0e13));

time_step = 0.0004001; %0.00041; %0.000041; %0.0021; %0.00021; %0.00051; %0.000051; %05;
total_time = 0.4; %5; %0.3; %2; %5; %2; %0.6; %0.6; %1; %8; %10.0; %5.0; %1.0; %0.1; %0.25; %0.5; %2; %5; %1; %3; %0.3; %3; % 0.3; %3; %0.3; %5e-3; % %0.1; %0.6; %0.6; %0.0006; % % second  
start_time = -total_time/2; %-1; -2.5; %-0.5; %-1.5; %-0.15; %-1.5; %-0.15; %-1.5; %-0.15; %-0.05; % %-0.05; %0.3; %-0.3; %-0.0003; % %
end_time = start_time + total_time;

%number of samples used for diffraction modelling of target
Nx = 401; %0; %00; %1000;
Nz = 401; %0; %00; %1000;
%N_F = Nx * Ny + 1;
%%%%%%%%%

%antenna size 0.4m x 0.4m
%antenna_halfwidth = 0.22; %0.2;
%antenna_halflength = 0.22; % 0.2;

wave_select = 1; %2; %3; %
wave_set = [1.0e-2, 2.0e-2, 0.5e-2];
band = [30, 15, 60];
wavelength =  wave_set(wave_select); %1.0e-2; %30GHZ 0.5e-2; %60GHz %2.0e-2; %15GHZ %


% offset from x=0

Num_data_point = length([start_time : time_step : end_time]); %round(total_time/time_step,0); 

%setup the elevation angle

%alpha_initial = 0.089999999999999999999999999999999999; %elevation angle: 40 degrees
%alpha_initial = 0; %elevation angle: 90 degrees
alpha_initial = 0.03674999999999; %elevation angle: 65 degrees

%[1/60:-1/120:-1/60]*pi; %[1/30:-1/120:-1/30]*pi; % [1/60:-1/2400:-1/60]*pi; %  %[1/300:-1/600:-1/300]*pi
theta_ini = pi/2 - atan(((orbit_height + r_Earth)*sin(alpha_initial))/(((orbit_height+r_Earth)*cos(alpha_initial)) - r_Earth));
theta_initial_degree = theta_ini/pi * 180

test1 = strcat(num2str(total_time),'seconds');
test2 = strcat(num2str(round(theta_initial_degree,2)), 'degree');
test3 = strcat(test1,test2);
test4 = strcat(num2str(plane_speed),'mps');
test5 = strcat(test3,test4);
%total_time1 = round(total_time*10);
%test51 = strcat(test5,num2str(data_length));
savefilename = strcat(test5,'signal02seconddroneset250mintotalDJI100msmallupdate064times081.mat');
%Note: the speed of the airplane is in the yG direction (i.e., parallel with the ground in the ground static coordinate)



 %plane location along with time
    distanceY1st = -1.0 * ([start_time : time_step : end_time]) * plane_speed + plane_height/tan(theta_ini); % relative to the reference ground station (y=0), the minus sign due to the distance at the beginning is set to be positive          
    distanceY = distanceY1st; % - Gstation_location(m_Y); %distance relative to the line perpendicular to the ground station of interest
%x_track_initial = 201.1; %10.1; %20.1; %30.1; %50.1; %100.1; %10300.1; %301.1; %201.1; %301.9; %290; %280; %51; %301.9; %240; %280; %301.9; %51; %101; %181; %200; %240; %270; %290.1; %301.9; %55.9; %290; %280; %240; %301.9; %240; %120; %25.9; %105.5; %105.5; %0; %159.5; %


%satellite elevation angle 
alpha = -1.0 * [start_time : time_step : end_time] * 2 * pi/T_satellite + alpha_initial; % - Gstation_location(m_Y)/r_Earth;
theta = pi/2 - atan(((orbit_height+r_Earth)*sin(alpha))./(((orbit_height+r_Earth)*cos(alpha)) - r_Earth)); %the elevation angle of satellite relative to the ground station of interest

%satellite elevation angle change rate
temp = (orbit_height +  r_Earth);     
Vdthetadt =  (temp * r_Earth * cos(alpha) - temp * temp)./(temp * temp - 2 * temp * r_Earth * cos(alpha) + r_Earth * r_Earth)  * (- 2 * pi / T_satellite);

%coordinate attached with the ground
%Vx_G = zeros(size(theta)) + x_track_initial;
Vy_G =  distanceY;                              % plane_height ./ tan(theta);
Vz_G = plane_height * ones(size(theta));

% Plane parameters in coordinate attached to satellite_ground link
%Vx_F = Vx_G;   
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

Vy_S1 = sqrt((orbit_height + r_Earth)^2 - (r_Earth * cos(theta)).^2 ) - r_Earth * sin(theta);
Vy_S = Vy_S1; % - d_Sat;

Vv_yS = - cos(theta + alpha) * 2 * pi * (r_Earth + orbit_height) / T_satellite;
Vwidth_z = width_z * sin(theta);

output_signalset = zeros(Num_data_point,length(x_initial_set));
output_N0_set = output_signalset;
Vdoppler_set = output_signalset;
Vbound_set = output_signalset;


for index = 1:1:length(x_initial_set)
    index
    x_track_initial = x_initial_set(index);

    output_intensity = zeros(Num_data_point,1);
    output_FK = output_intensity;
    output_N0 = output_intensity;
    outputgain_N1 = output_intensity;
    outputgain_N0 = output_intensity;
    outputgain_target = output_intensity;
       
    %coordinate attached with the ground
    Vx_G = zeros(size(theta)) + x_track_initial;
    %Plane parameters in coordinate attached to satellite_ground link
    Vx_F = Vx_G; 
   
    % distance related to the center of the target
    r1 = sqrt(Vx_G .* Vx_G + Vy_G .* Vy_G + Vz_G .* Vz_G);  %sqrt(plane_height * plane_height + distanceY .* distanceY);
    %r1 = sqrt(plane_height * plane_height + plane_height*plane_height/tan(theta_ini)/tan(theta_ini));

    %r2 = sqrt((orbit_height + r_Earth)^2 - (r_Earth * cos(theta_ini)).^2 ) - r_Earth * sin(theta_ini) - r1;
    %r2 =  Vy_S1 - r1;
     r2 = sqrt(Vx_F .* Vx_F + (Vy_S1 - Vy_F) .* (Vy_S1 - Vy_F) + Vz_F .* Vz_F);
  
    Vdoppler3 = 1/wavelength * (Vx_F * v_x + Vy_F .* Vv_y + Vz_F .* Vv_z) ./ r1;
    Vdoppler4 = 1/wavelength * ((Vx_F * v_x + (Vy_F - Vy_S) .* (Vv_y - Vv_yS) + Vz_F .* Vv_z) ./ r2 - Vv_yS);

    %Vdoppler = Vdoppler1 + Vdoppler2 + Vdoppler3 + Vdoppler4;

     Vdoppler = Vdoppler3 + Vdoppler4;
     Vbound = (1 ./ r1 + 1 ./ r2)/wavelength .* (abs(length_x * v_x) + abs(Vwidth_z .* Vv_z))/2;
    
     for m = 1:1:Num_data_point
         %m
         width_ztemp = Vwidth_z(m);
         %distancez = plane_height/sin(theta(m)) + cos(theta(m)) * (distanceY(m)  - plane_height / tan(theta(m))); %+ plane_height /tan(theta_ini));
         
         %Vy_F(m) %r1(m)
         %dPlane = plane_height/sin(theta(m)) + (distanceY(m)  - plane_height / tan(theta(m))) * cos(theta(m)) %projected onto the optical axis, distance between flying object and groud receiver 
         dPlane = Vy_F(m);
         
         %dLEO = sqrt((orbit_height + r_Earth)^2 - (r_Earth * cos(theta(m)))^2 ) - r_Earth * sin(theta(m));     % distance between satellite and ground station
         dLEO = Vy_S(m);
         %Vy_S(m) - Vy_F(m)
         %r2(m)
         %dLEO - dPlane
         
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
         outputgain_N1(m) = RXgain;
         
         outputgain_target(m) = tttemp/dLEO * orbit_height;
         outputgain_N0(m) = exp(1i * 2 * pi * dLEO /wavelength)/dLEO * orbit_height; %normalized by LEO overhead sqrt(pathloss);
         
         tempFK = RXgain * conj(RXgain); % * pathloss;
         %tempFK_am = (1 + tttemp) * pathloss_am;
         output_FK(m) = sqrt(tempFK);
         output_N0(m) = sqrt(outputgain_N0(m) * conj(outputgain_N0(m))); %(1 + 0.0) * sqrt(pathloss);
         
         %here about digital communications???
         
         %{
         data = randi([0 M-1],data_length,1);
         modData = genqammod(data,c);
         rxSig = awgn(modData*RXgain,SNR,'measured','db');
         output_FKMPSK(m) = mean(abs(rxSig));
         %}
         % = tempFK;
         % obtain the amplitudes and averaging out
     end
     output_signalset(:,index) = output_FK;
     output_N0_set(:,index) = output_N0;
     Vdoppler_set(:,index) = Vdoppler;
     Vbound_set(:,index) = Vbound;
end
     
 %    intensity_all_output1(m_S, m_Y,:) = output_intensity;
          
    % output_intensity1 = output_intensity/mean(output_intensity);
    % output_dB = 10*log10(output_intensity1);
    % output_detrend = detrend(output_dB);  %remove the linear component of the signal intesnity
          
     % var_theta(m_Y,i_alpha_ini) = std(output_detrend)
  figure
  
  x_time = [start_time : time_step : end_time];
  
  y_FK = output_FK;
 % y_intensity = 10*log10(output_intensity);
 %y_FK = output_FKMPSK; %
 % y_FK = 20*log10(output_FKMPSK);
  %y_FK = 20*log10(output_FK);
  save(savefilename,'x_initial_set','x_time','output_signalset','output_N0_set','Vdoppler_set','Vbound_set'); %,'outputgain_N1','outputgain_target','outputgain_N0');
  
  mesh(x_time, x_initial_set, output_signalset');
  xlabel('Time (second)');
  ylabel('Minimum distances to links (metres)');
  zlabel('Signal');
 
  figure
  plot(x_time, y_FK, 'b',x_time,output_N0,'.k');
  % y_FK_detrend = detrend(y_FK,'linear');
  % plot(x_time, y_FK_detrend, 'b'); %,x_time,output_N0,'.k');
legend('Fresnel Kirchhoff','channel gain in power');
xlabel('time (second)');
ylabel('received signal (linear)');
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

figure
plot(x_time, abs(Vdoppler),'-k', x_time, abs(Vdoppler + Vbound), '--k', x_time, abs(Vdoppler - Vbound),'-.k');
legend('Median','Upper Bound','Lower Bound');
%title(strcat(strcat(num2str(band(wave_select)),strcat(strcat(" GHz    ", num2str(round(theta_initial_degree,2))), ' degree')), strcat("   ", strcat(num2str(plane_speed),' m/sec')))); 
%save(savefilename,'width_z','length_x','x_time','output_FK','output_FKMPSK','Vdoppler','Vbound');
%save(savefilename,'width_z','length_x','x_time','output_FK','output_FH_F','output_FH','output_FN','Vdoppler','Vbound');
%save(savefilename,'x_time','y_FK','y_FH_F','y_FH','y_FN_F','y_FN','Vdoppler'); %, 'doppler_slope','doppler_offset');
 %save('FK90degree2010update.mat','x_time','y_FK','y_FH_F','y_FH','y_FN_F','y_FN');
  %alpha_initial = alpha_initial - alpha_between_satellite;
%end
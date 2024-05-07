clear all

plane_height = 1.0e4 + 0.051;
orbit_height = 5.5e5; %3.0e5;
r_Earth = 6.371e6;

T_satellite = 2*pi * sqrt((orbit_height + r_Earth)^3/((5.972*6.6723) * 1.0e13));

plane_speed = 900; %600; %200; %50; %200; %100; %200; %200; %500; %200; %125; %250;  %250m/second = 900km/hour

time_step = 0.01; %0.002; %0.0000002; % second per step
total_time = 280; %40; %0; %0.02; %0.2; %1; %0.3; % 0.5; %1; %0.5; %0.5 second
start_time = -140; %-20;
end_time = start_time + total_time;
Satellite_in_total = 360; %0;
alpha_between_satellite = 2 * pi /Satellite_in_total;

%d_satellite = (r_Earth + plane_height)*sin(alpha_between_satellite/2)*2 % at plane height, the distance between two links of neighboring satelites

Num_satellite = 5; %3;
Num_groundstationtemp = 2000; %1600; %800; %200; %1500; %500; %ls5; %10;

d_station = 50; %100; %40;% 100; %distance between neighbour ground stations in metres

Gstation_location = -[-Num_groundstationtemp/2:Num_groundstationtemp/2]*d_station; %location of Ground Stations
Num_groundstation = length(Gstation_location);

%%%%%%%%
GSSAMPLE = 1 %64 %16
%%%%%%%%

alpha_station = d_station/r_Earth; % angle (relative to center of Earth) between two neighbour ground stations

%alpha_initial = alpha_between_satellite; %0.0; %Reference (virtual) line. The line perpendicular to the ground at the first ground station. The flying object and satellites are all referenced to this line.

%[1/60:-1/120:-1/60]*pi; %[1/30:-1/120:-1/30]*pi; % [1/60:-1/2400:-1/60]*pi; %  %[1/300:-1/600:-1/300]*pi
%theta_ini = pi/2 - atan(((orbit_height + r_Earth).*sin(alpha_initial))./(((orbit_height+r_Earth).*cos(alpha_initial)) - r_Earth));
%theta_initial_degree = theta_initial/pi * 180

load tempfile600; %tempfile_800G_50m_360S_5;
temp1 = size(intensity_all_output1);
end_time = start_time + total_time;

output_dB = 10*log10(intensity_all_output1);

sampletime = 0.2;
tempinteger = ceil(sampletime/time_step);
breakpoint = tempinteger:tempinteger:length([start_time : time_step : end_time]);

Num_samples = floor((length([start_time : time_step : end_time]) - 1)/(tempinteger-1));

Y_max = zeros(temp1(1),Num_samples);
I_max = Y_max;

for sn = 1:temp1(1)
    %{
    figure
    hold on
    %}
    temp = reshape(output_dB(sn,:,:),temp1(2),temp1(3)); % the nth row is for the received signal of the nth ground station
    output_detrend = detrend(temp','linear',breakpoint);   % nth column is for the detrended signal of the nth ground stationsiz 
    output_temp = output_detrend(2:(1+Num_samples*(tempinteger-1)),:);
    output_temp1 = reshape(output_temp, tempinteger-1, Num_samples*temp1(2));
    output_std_temp = reshape(std(output_temp1), Num_samples, temp1(2));
    output_std = output_std_temp(:,1:GSSAMPLE:temp1(2));
    [Y_max(sn,:), I_max(sn,:)] = max(output_std'); 
    
end

Y_sat = zeros(Num_samples,Num_satellite);
Z_sat = zeros(Num_samples,Num_satellite);

Num_satellite = temp1(1);

x_plane = zeros(Num_samples,2);
     %distanceY1st = -1.0 * ([start_time : time_step : end_time]) * plane_speed; % relative to the first ground station, the minus sign due to the distance at the beginning is set to be positive          
for i_Y = 1:1:Num_samples  %Num_groundstation
   
    A = zeros(2*Num_satellite,Num_satellite + 2);
    g = zeros(2*Num_satellite,1);
    Vn = zeros(2*Num_satellite);
    
   alpha_initial = alpha_between_satellite;
   for m_S = 1:1:Num_satellite
         m_Y = (I_max(m_S,i_Y) - 1) * GSSAMPLE + 1; %which ground station has the maximum
        %distanceY = distanceY1st + Gstation_location(m_Y); %distance relative to the line perpendicular to the ground station of interest
        alpha1 = -1.0 * start_time * 2 * pi/T_satellite - 1.0 * ((tempinteger - 1) * (i_Y - 1) + (tempinteger - 1)/2) * time_step * 2 * pi/T_satellite + alpha_initial - Gstation_location(m_Y)/r_Earth;
        theta1 = pi/2 - atan(((orbit_height+r_Earth).*sin(alpha1))./(((orbit_height+r_Earth).*cos(alpha1)) - r_Earth)); %the elevation angle of satellite relative to the ground station of interest
        d_SG = sqrt((orbit_height + r_Earth)^2 - (r_Earth * cos(theta1))^2 ) - r_Earth * sin(theta1);
        Y_sat = Gstation_location(m_Y)/d_SG -  cos(-theta1);
       % Gstation_location(m_Y) - Y_sat(i_Y,m_S);
       % Z_sat(i_Y,m_S) = d_SG * sin(theta1);
        
        A((m_S-1)*2+1,1) = 1;
        A((m_S-1)*2+2,2) = 1;
        A((m_S-1)*2+1,m_S+2) =  -cos(theta1); % - d_SG * cos(theta1); %the y direction, the '-' due to to the right is positive, and to the left is negative 
        A((m_S-1)*2+2,m_S+2) = -sin(theta1); %d_SG * sin(theta1); % the z direction
        
        g((m_S-1)*2+1) =  Gstation_location(m_Y);
        g((m_S-1)*2+2) = 0;
        
        Vn((m_S-1)*2+1,(m_S-1)*2+1) = Y_max(m_S, i_Y) * Y_max(m_S, i_Y); %the inverse of Vn actually
        Vn((m_S-1)*2+2,(m_S-1)*2+2) = Y_max(m_S, i_Y) * Y_max(m_S, i_Y);
         
        alpha_initial = alpha_initial - alpha_between_satellite;
   end
   
   btemp = inv(A' * Vn * A) * A' * Vn * g; %pinv(A) * g; % % % %pinv(A) * g; % % %pinv(A) * g; %
   
   x_plane(i_Y,1) = btemp(1);
   x_plane(i_Y,2) = btemp(2);
end
%(figure
 %  hold on  
 %plot(x_plane(:,1),x_plane(:,2))
%}
 % drawing the output?

 TIME1 = start_time + [0:Num_samples-1] * (tempinteger - 1) * time_step + (tempinteger - 1)/2 * time_step;
 
 plane_locationx = - TIME1 * plane_speed;
 
 plane_locationy = ones(Num_samples,1)* plane_height;

figure
hold on
%for  sn = 1:temp1(1)
       plot3(TIME1,x_plane(:,1), x_plane(:,2));  %sn*ones(Num_samples,1),Y_max(sn,:))
       plot3(TIME1, plane_locationx, plane_locationy);
      %plot3(1:Num_samples,sn*ones(Num_samples,1),I_max(sn,:))
%end
xlabel('time');
ylabel('x axis');
zlabel('height');
D_STA = GSSAMPLE * d_station;
%title(strcat('Sample Duration ', num2str(sampletime)));
title(strcat('distance between ground stations in metres = ', num2str(D_STA)));
 view(3)
         
%{
for sn = 1:temp1(1)
    %figure 
    %hold on
    
for i = 650:800 %1:temp1(2) %1:50 % %100:110 %50 %1:50  %
   plot3([start_time : time_step : end_time],i*ones(temp1(3),1),reshape(intensity_all_output1(sn,i,:),temp1(3),1))
  
end
%view(3)
end
%image([1:1:temp1(2)],[1:1:temp1(3)],temp)
%}
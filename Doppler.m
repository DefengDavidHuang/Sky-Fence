clear all

plane_height = 1.0e4;
orbit_height = 5.5e5; %3.0e5;
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

time_step = 0.00002;
total_time = 6; %0.6; %0.0006; % % second  
start_time = -3; %-0.3; %-0.0003; % %
end_time = start_time + total_time;

%%%%%%%%%

%antenna size 0.4m x 0.4m
antenna_halfwidth = 0.22; %0.2;
antenna_halflength = 0.22; % 0.2;

p1 = 0.1;
p2 = 0.1;

wave_select = 1; %2; %3; %
wave_set = [1.0e-2, 2.0e-2, 0.5e-2];
band = [30, 15, 60];
wavelength =  wave_set(wave_select); %1.0e-2; %30GHZ 0.5e-2; %60GHz %2.0e-2; %15GHZ %

width_y = 10; %20; %20; %20.0;
length_x = 10; %20; %20.0;    %x: parallel to ground and perpendicular to the satellite link

Num_data_point = round(total_time/time_step,0); 
%

%alpha_station = d_station/r_Earth; % angle (relative to center of Earth) between two neighbour ground stations

distanceX = zeros(Num_data_point+1,1); % + 10; % the movement of flying object in X direction

%alpha_initial = 0.089999999999999999999999999999999999; %elevation angle: 40 degrees

alpha_initial = 0; %elevation angle: 90 degrees

%alpha_initial = 0.03674999999999; %elevation angle: 65 degrees

%[1/60:-1/120:-1/60]*pi; %[1/30:-1/120:-1/30]*pi; % [1/60:-1/2400:-1/60]*pi; %  %[1/300:-1/600:-1/300]*pi
theta_ini = pi/2 - atan(((orbit_height + r_Earth)*sin(alpha_initial))/(((orbit_height+r_Earth)*cos(alpha_initial)) - r_Earth));
theta_initial_degree = theta_ini/pi * 180;

doppler = - plane_speed * cos(theta_ini)/wavelength
doppler_offset = plane_speed * plane_speed * total_time/2/wavelength/plane_height*sin(theta_ini)

%var_theta = zeros(Num_groundstation,length(alpha_initial));

%intensity_all_output1 = zeros(Num_satellite,Num_groundstation,Num_data_point+1);
 
output_intensity = zeros(Num_data_point+1,1);
% distanceY = ([-total_time/2 : time_step : total_time/2]) * plane_speed; % - 1000;

%for m_S = 1:1:Num_satellite
     
    % alpha = -1.0 * [start_time : time_step : end_time] * 2 * pi/T_satellite + alpha_initial;   %the alpha angle relative to the first ground station. Note the minus is for, at the beginning the alpha is positive, which corrsponds to negative time
     %distanceY1st = -1.0 * ([start_time : time_step : end_time]) * plane_speed - d_satellite * (m_S -1); % relative to the first ground station, the minus sign due to the distance at the beginning is set to be positive
     
     distanceY1st = -1.0 * ([start_time : time_step : end_time]) * plane_speed + plane_height/tan(theta_ini); % relative to the reference ground station (y=0), the minus sign due to the distance at the beginning is set to be positive
     
     %m_S
     
     %theta_plane = pi/2 - atan(((orbit_height + r_Earth + plane_height)*sin(alpha_initial))/(((orbit_height + r_Earth + plane_height)*cos(alpha_initial)) - r_Earth - plane_height))
     %Y_satellite_position = plane_height/tan(theta_plane)
     
    % for m_Y = 1:1:Num_groundstation
     
      %   m_Y   
         
     distanceY = distanceY1st; % - Gstation_location(m_Y); %distance relative to the line perpendicular to the ground station of interest
     
     alpha = -1.0 * [start_time : time_step : end_time] * 2 * pi/T_satellite + alpha_initial; % - Gstation_location(m_Y)/r_Earth;
    
     theta = pi/2 - atan(((orbit_height+r_Earth).*sin(alpha))./(((orbit_height+r_Earth).*cos(alpha)) - r_Earth)); %the elevation angle of satellite relative to the ground station of interest
     
     %theta_planetemp = pi/2 - atan(((orbit_height+r_Earth + plane_height).*sin(alpha))./(((orbit_height + r_Earth + plane_height).*cos(alpha)) - r_Earth)); %the elevation angle of satellite relative to the ground station
     %Y_satellite_positiontemp = plane_height./tan(theta_planetemp)
     
     %theta(Num_data_point+1)/pi*180
     
     
       for m = 1:1:Num_data_point+1
         width_ytemp = width_y * sin(theta(m));
         distancez = plane_height/sin(theta(m)) + cos(theta(m)) * (distanceY(m)  - plane_height / tan(theta(m))); %+ plane_height /tan(theta_ini));
         in = [wavelength, width_ytemp, length_x, distancez];
         
         tempY = sin(theta(m)) * (distanceY(m)  - plane_height / tan(theta(m)));    % distanceY(m) * sin(theta(m)) - plane_height * cos(theta(m)) + plane_height * sin(theta(m))/tan(theta_ini);
         Y = ((tempY - antenna_halflength):p1:(tempY + antenna_halflength));  %Note: find proper Y and X, and do the averaging directly, by updating the upper limit and the lower limit 
         X = ((distanceX(m) - antenna_halfwidth):p2:(distanceX(m) + antenna_halfwidth));
         
         output = snapshot_function(in, Y, X);
         temp = orbit_height*orbit_height/(sqrt((orbit_height + r_Earth)^2 - (r_Earth * cos(theta(m)))^2 ) - r_Earth * sin(theta(m)))^2; %deal with distance difference.
         output_intensity(m) = mean(output, 'all') * temp;
     end
     
 %    intensity_all_output1(m_S, m_Y,:) = output_intensity;
     
     
    % output_intensity1 = output_intensity/mean(output_intensity);
    % output_dB = 10*log10(output_intensity1);
    % output_detrend = detrend(output_dB);  %remove the linear component of the signal intesnity
     
     
     % var_theta(m_Y,i_alpha_ini) = std(output_detrend)
  
 
  figure
  
  x_time = [start_time : time_step : end_time];
 % y_intensity = 10*log10(output_intensity);
  
 y_intensity = output_intensity;
 
 spectrogram(y_intensity,128,20,128,1/time_step);
 

 % plot(x_time, y_intensity);
  
  %plot([start_time : time_step : end_time], output_detrend); %dB); %detrend);
  
  %title(['Satellite #', num2str(m_S), 'Ground Station #', num2str(m_Y)]);
  %xlabel('time (second)');
  %ylabel('normalized intensity (dB)');
  %ylabel('intensity');
  title(strcat(strcat(num2str(band(wave_select)),strcat(strcat(" GHz    ", num2str(round(theta_initial_degree,2))), ' degree')), strcat("   ", strcat(num2str(plane_speed),' m/sec')))); 
  
  %axes('position',[.65 .175 .25 .25])
%box on % put box around new pair of axes
%indexOfInterest = (t < 11*pi/8) & (t > 9*pi/8); % range of t near perturbation
%plot([start_timezoom:timestep:end_timezoom,signal(indexOfInterest)) % plot on new axes
%axis tight
  
  %alpha = alpha - alpha_station; % note that the left side is negative angle, so the minus (-) sign
  
     %end
 
  save('onesatellite.mat','output_intensity','start_time','total_time','time_step','distanceY1st','plane_height');
  %alpha_initial = alpha_initial - alpha_between_satellite;
%end

%Data Processing

%figure
%plot(theta_initial_degree, var_theta);
%xlabel('elevation angle (degree)');
%ylabel('normalized standard deviation (dB)');

%axis([38 142 1.5 9]);
%plot([time_step:time_step:total_time], output_detrend); %10*log10(output_detrend)); 

%I1 = I/max(max(I))*64;
%image(Y, X, output)
%colorbar
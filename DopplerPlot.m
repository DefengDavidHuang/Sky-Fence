%clear all
%doppler_offset = plane_speed * plane_speed * total_time/2/wavelength/plane_height*sin(theta_ini)

%var_theta = zeros(Num_groundstation,length(alpha_initial));

%intensity_all_output1 = zeros(Num_satellite,Num_groundstation,Num_data_point+1);
 
%load FK65degree2010update
%load 30GHz65degree250mpsSameOrbit
%load 30GHz90degree125mpsSameOrbit
%load 30GHz90degree250mpsSameOrbit
%load 30GHz90degree250mpsSameOrbit
%load 30GHz65degree125mpsCross 
%load 30GHz90degree125mpsCross 

%load 30GHz90degree125mpsCross
%load 30GHz90degree250mpsCross    %FK2010cross90degree0mps
%load 30GHz65degree250mpsCross
%load 30GHz65degree250mps09Cross
%load 30GHz90degree250mps9Cross

%load 30GHz65degree250mps6Cross
load 
%load 30GHz65degree250mps30SameOrbit
%load 30GHz90degree250mps30SameOrbit
%load 30GHz40degree250mps30SameOrbit

%load 30GHz90degree125mps9SameOrbit
%load 30GHz90degree250mps9SameOrbit
%load 30GHz65degree250mps9SameOrbit
%load 30GHz65degree250mps30SameOrbit

%load FK90degree2010update
%load FK65degree2010cross65degree
%load FK65degree2010cross90degree
  
% y_intensity = 10.^(y_FN_F/10);
 
  y_intensity = 10.^(y_FK/10);
 time_step = x_time(2) - x_time(1);
 
 samples = 3;
 y_intensity_update = y_intensity(1:samples:length(y_intensity));
 
 [s,w,t] = spectrogram(y_intensity_update,128,127,256,1/(time_step*samples),'yaxis');
  spectrogram(y_intensity_update,128,127,256,1/(time_step*samples),'yaxis','centered');
  h = gca;
  h.XTickLabel = string(h.XTick + min(x_time));  %note in seconds, no x 1000 is required
  

 doppler_slope
 
 hold on
 t1 = t + min(x_time);
 plot(t, doppler_slope*t1,'-k');
 hold off
 
%spectrogram(y_intensity,128,110,1024,1/time_step);
 %spectrogram(y_intensity_update,1024,1000,1024,1/(time_step*samples));
 

 % plot(x_time, y_intensity);
  
  %plot([start_time : time_step : end_time], output_detrend); %dB); %detrend);
  
  %title(['Satellite #', num2str(m_S), 'Ground Station #', num2str(m_Y)]);
  %xlabel('time (second)');
  %ylabel('normalized intensity (dB)');
  %ylabel('intensity');
  %title(strcat(strcat(num2str(band(wave_select)),strcat(strcat(" GHz    ", num2str(round(theta_initial_degree,2))), ' degree')), strcat("   ", strcat(num2str(plane_speed),' m/sec')))); 
  
  %axes('position',[.65 .175 .25 .25])
%box on % put box around new pair of axes
%indexOfInterest = (t < 11*pi/8) & (t > 9*pi/8); % range of t near perturbation
%plot([start_timezoom:timestep:end_timezoom,signal(indexOfInterest)) % plot on new axes
%axis tight
  
  %alpha = alpha - alpha_station; % note that the left side is negative angle, so the minus (-) sign
  
     %end
 
 % save('onesatellite.mat','output_intensity','start_time','total_time','time_step','distanceY1st','plane_height');
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
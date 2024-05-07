%clear all


%No NOISE

%load 30GHz40degree-250mps30SameOrbitAmplitudeMPSK
%load 30GHz65degree-250mps30SameOrbitAmplitudeMPSK
%load 30GHz90degree-250mps30SameOrbitAmplitudeMPSK
%load 30GHz90degree250mps30SameOrbitAmplitudeMPSK
%load 30GHz65degree250mps30SameOrbitAmplitudeMPSK
%load 30GHz40degree250mps30SameOrbitAmplitudeMPSK

%with Noise and QAM
%100samples

%load('0offsetinX90degree100.1mpssignal02secondlower.mat'); % 30GHz65degree250mps30SameOrbitAmplitudeQAMn100
%load('20.1offsetinX90degree100.1mpssignal02secondlower.mat');
load('20.1offsetinX65degree100.1mpssignal02secondlower.mat');
%load('20.1offsetinX40degree100.1mpssignal02secondlower.mat');

%load 0dB65degree250mps1000SameOrbitAmplitudeQAM
%load 10dB65degree250mps100SameOrbitAmplitudeQAM
%load 10dB65degree250mps1000SameOrbitAmplitudeQAM
%load 20dB65degree250mps100SameOrbitAmplitudeQAM
%load 20dB65degree250mps1000SameOrbitAmplitudeQAM
 
% y_intensity = output_FK; %10.^(y_FK/10);
  y_intensity = y_FK; %output_FKMPSK;
 time_step = x_time(2) - x_time(1);
 
 figure
 plot(x_time, y_intensity, '-k');
 xlabel('Time (secs)');
 ylabel('Amplitude of received signal');
 
 samples = 4; %3;
 y_intensity_update = y_intensity(1:samples:length(y_intensity));
 
 figure
 %
 lengthFFT = 128; %64; %128; %64; %128; %64; %256; %128; %64; %256; %2048; %512;
 [s,w,t] = spectrogram(y_intensity_update,lengthFFT,lengthFFT-1,lengthFFT*2,1/(time_step*samples),'yaxis');
 spectrogram(y_intensity_update,lengthFFT,lengthFFT-1,lengthFFT*2,1/(time_step*samples),'yaxis');
  %spectrogram(y_intensity_update,128,127,256,1/(time_step*samples),'yaxis','centered');
  h = gca;
  h.XTickLabel = string(h.XTick + min(x_time));  %note in seconds, no x 1000 is required
 
% figure 
% plot(x_time, abs(Vdoppler),'-k', x_time, abs(Vdoppler + Vbound), '--k', x_time, abs(Vdoppler - Vbound),'-.k'); %plot(x_time, abs(Vdoppler),'-k');
%figure 
st = lengthFFT/2; %floor((length(Vdoppler) - length(t))/2);
 hold on
  t1 = t; % + min(x_time);
 range = [st:1:st+length(t1)-1]*samples;
  
  plot(t1, abs(Vdoppler(range)),'-k', t1, abs(Vdoppler(range) + Vbound(range)), '--k', t1,abs(Vdoppler(range) - Vbound(range)) ,'-.k'); %plot(x_time, abs(Vdoppler),'-k');
  %plot(t1, abs(Vdoppler(range))/1000,'-k', t1, abs(Vdoppler(range) + Vbound(range))/1000, '--k', t1,abs(Vdoppler(range) - Vbound(range))/1000 ,'-.k');
  legend('Doppler shift','bound 1','bound 2');
   
 hold off
 

 %}
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
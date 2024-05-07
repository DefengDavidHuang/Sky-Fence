clear all

%LoverA = [0.01, 0.02, 0.05, 0.1];

alphadegree = 3.5;
alpha = alphadegree/180*pi;

%h = 500; %100;
h = 100;

%S = [500 800 1000 1500]; %[100 400 800 1500];
S = [100 200 400 800]; 

LoverA = alpha * h ./ S;

K = 1:1:100;

p_D = ones(length(LoverA),length(K));

for n = 1:1:length(LoverA)
   p_D(n,:) = 1 - (1 - LoverA(n)).^K;
end
    
plot(K, p_D(1,:), '-k', K, p_D(2,:), '-.k', K, p_D(3,:),'--k', K, p_D(4,:),':k');

%legend("S = 500 meters", "S = 800 meters", "S = 1000 meters", "S = 1500 meters");
legend("S = 100 meters", "S = 200 meters", "S = 400 meters", "S = 800 meters");
xlabel('Number of links');
ylabel('Probability of Dection');
%plot(time, kai1, '-', time, kai2, 'x')






     %theta_planetemp = pi/2 - atan(((orbit_height+r_Earth + plane_height).*sin(alpha))./(((orbit_height + r_Earth + plane_height).*cos(alpha)) - r_Earth)); %the elevation angle of satellite relative to the ground station
     %Y_satellite_positiontemp = plane_height./tan(theta_planetemp)
     
     %theta(Num_data_point+1)/pi*180
   

     

     
     %{
     output_intensity1 = output_intensity/mean(output_intensity);
     output_dB = 10*log10(output_intensity1);
     output_detrend = detrend(output_dB);  %remove the linear component of the signal intesnity
     %}
     
     % var_theta(m_Y,i_alpha_ini) = std(output_detrend)
  
 %{
  figure
  plot([start_time : time_step : end_time], output_detrend); %dB); %detrend);
  title(['Satellite #', num2str(m_S), 'Ground Station #', num2str(m_Y)]);
  xlabel('time (second)');
  ylabel('normalized intensity (dB)');
 %}
  
  %alpha = alpha - alpha_station; % note that the left side is negative angle, so the minus (-) sign
  
 % save('tempfile1.mat','intensity_all_output1','Num_satellite', 'Num_groundstation', 'd_station','start_time','total_time','time_step','distanceY1st','plane_height');
 % alpha_initial = alpha_initial - alpha_between_satellite;
 
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
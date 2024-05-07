%load tempfile_500station_200m_360satellite;
%load tempfile_500station_50m_3600satellite;
%load tempfile_500station_100m_360satellite_3;

%load tempfile200station_200m_360satellite_2;

clear all
load tempfile600


%load tempfile_800G_50m_360S_5;
temp1 = size(intensity_all_output1);
%total_time = 280; %0; %0.02; %0.2; %1; %0.3; % 0.5; %1; %0.5; %0.5 second
%start_time = -140;
end_time = start_time + total_time;
%temp = reshape(intensity_all_output1(1,:,:),temp1(2),temp1(3));
%temp = reshape(intensity_all_output1(2,:,:),temp1(2),temp1(3));

time1 = 50000;
test = 1:time1;
%sn = 3;  %sn=3 920:990
%sn = 2; %965:1035
%sn = 1;  %1010:1080

figure                      
hold on
for sn = 1:temp1(1)
for i = 1:temp1(2) %900:1100 %965:1035 %500:1500 %920:990  %800:1200 %1:temp1(2) %1:50 % %100:110 %50 %1:50  %
   %plot3(test,i*ones(time1,1),reshape(intensity_all_output1(sn,i,test),time1,1))
  % plot3([start_time : time_step : end_time][1:1:temp1(3)]*time_step,i*ones(temp1(3),1),reshape(intensity_all_output1(sn,i,:),temp1(3),1))
  plot3([start_time : time_step : end_time],i*ones(temp1(3),1),reshape(intensity_all_output1(sn,i,:),temp1(3),1))
end
end
view(3)
%image([1:1:temp1(2)],[1:1:temp1(3)],temp)
%}
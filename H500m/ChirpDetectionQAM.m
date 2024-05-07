clear all


%load('301.9offsetinX65degree250mpssignal.mat');
%load('280offsetinX65degree250mpssignal.mat');
%load('240offsetinX65degree250mpssignal.mat');
%load('101offsetinX65degree250mpssignal.mat');
%load('51offsetinX65degree250mpssignal.mat');
%INPUT from the file: outputgain_N1, outputgain_N0
load('201.1offsetinX65degree100.1mpssignal02second.mat');

Num_data_point = length(x_time);

%{
figure
plot(x_time, y_FK, 'b',x_time,output_N0,'.k');
legend('Fresnel Kirchhoff');
xlabel('time (second)');
ylabel('received signal (linear)');

figure
y_FK_detrend = detrend(y_FK,'linear');
output_N0_detrend = detrend(output_N0,'linear');
plot(x_time, y_FK_detrend, 'b',x_time,output_N0_detrend,'.k');
legend('N1','N0');
title('signal only');
%}

%modulation scheme
%c = [-5 -5i 5 5i -3 -3-3i -3i 3-3i 3 3+3i 3i -3+3i -1 -1i 1 1i]; % 16-QAM constellation based on the V.29 standard for telephone-line modems.


c = [-3-3i -3-1i -3+1i -3+3i -1-3i -1-1i -1+1i -1+3i 1-3i 1-1i 1+1i 1+3i 3-3i 3-1i 3+1i 3+3i]; % /sqrt(9.75);  %standard 16QAM
M = length(c);
temp = sqrt(c*c'/M);
c = c/temp; %16QAM, power normalize to 1


%c = [-1-1i -1+1i 1-1i 1+1i]/sqrt(2);
%M = length(c);

output_FKMPSK = y_FK;
output_MPSKN0 = y_FK;
SNR = 10; %20; %10; %3; %10; %10; %in dB
data_length = 10000; %100; %5000; %10000; %5000; %1000; %100; %10000; %1000; %5000; %1000; %10000; %1000; %10000; %10000; %5000; %1000; % 5000; % 50; %10; %100; %5000; %1000; %5000; %1000; %1000; %5000; %10000; %5000; %100; %5000; %1; %10000; %1000; %10000;
%Vec_T = 0.0025:0.0002:0.0050;
%Vec_T = 0.0028:0.0001:0.0055; %5000
%Vec_T = 0.020:0.0002:0.035;  %100
%Vec_T = 0.15:0.008:0.35;
%Vec_T = 0.002:0.0002:0.004
%Vec_T = 0.002:0.0004:0.006; %75; %0.01:0.0003:0.02; %0.015; %0.02; %0.015; %0.004:0.0001:0.005; %0.004
%Vec_T = 0.0060:0.0005:0.0090; %1000 symbols
%Vec_T = 0.0020:0.0001:0.003; %7; %10000 symbols 16QAM
%Vec_T = 0.0016:0.0002:0.004; %10000 symbols 16QAM
%Vec_T = 0.0055:0.0004:0.0090; %1000 symbols 16QAM
%Vec_T = 0.0028:0.0002:0.0050; %0.0020:0.0004:0.0065; %10,000 symbols
%Vec_T = 0.06:0.002:0.11;
%Vec_T = 0.035:0.001:0.055;
%Vec_T = 0.06:0.002:0.09;

%Vec_T = 0.0045; %0.003:0.0002:0.005; %0.004:0.0002:0.0080;
Vec_T = 0.011:0.0001:0.019;
%Vec_T = 0.0018:0.0001:0.005; %0.0030:0.0002:0.0094; %0.0030:0.0001:0.0045; %0.0068:0.0002:0.010; % % 0.010:0.001:0.018;
%Vec_T = 0.0045:0.0001:0.009;

%Vec_T = 0.0032:0.0001:0.0050;
%Vec_T = 0.0016:0.0001:0.0030; 
%Vec_T = 0.005:0.0002:0.010; 
%Vec_T = 0.0020:0.0001:0.0050;
%Vec_T = 0.0055:0.0002:0.0080;
%Vec_T = 0.0025:0.0001:0.0040;
%Vec_T = 0.03:0.001:0.045;
%Vec_T  = 0.0030:0.0001:0.0050; %0.0040:0.0004:0.0080;
%Vec_T = 0.0025:0.0001:0.003; %0.0020:0.0004:0.005; %7; %10000 symbols QPSK

Num_N1 = zeros(length(Vec_T),1);
Num_N0 = Num_N1;

Max_trials = 1e6;
Total_trials = Max_trials * ones(length(Vec_T),1);

Max_Num_errors = 100; %5; %10; %0;

var_noise = 10^(-SNR/10);
%var_total = (var(c) + var_noise)/data_length;
for j2 = 1:1:length(Vec_T)
 T = Vec_T(j2)
 %Num_N1
 %Num_N0
 for j1 = 1:1:Max_trials
  for m = 1:1:Num_data_point
     data = randi([0 M-1],data_length,1);
     modData = genqammod(data,c);
    
     noise = sqrt(var_noise/2) * (randn(data_length,1) + 1j * randn(data_length,1));
     RXgain = outputgain_N1(m);
     rxSig = modData*RXgain + noise;
     
     %rxSig = abs(modData) * abs(outputgain_N0(m)) + abs(modData) * real(outputgain_target(m)) + real(noise .* exp(-1i * phase(modData)));
     %rxSig = awgn(modData*RXgain,SNR,'measured','db');
     
     output_FKMPSK(m) = mean(abs(rxSig));
     %temp1 = mean(abs(rxSig))
     %tempN1 = real(RXgain) + mean(real(noise .* exp(-angle(modData))))
     
     RXgainN0 =  outputgain_N0(m);
     rxSigN0 = modData * RXgainN0 + noise;
     %rxSigN0 = abs(modData) * abs(RXgainN0) + real(noise .* exp(-1i * phase(modData)));
     
     %rxSigN0 = awgn(modData*RXgainN0,SNR,'measured','db');
     output_MPSKN0(m) = mean(abs(rxSigN0));
     %temp0 = mean(abs(rxSigN0))
     %tempN0 = RXgainN0 + mean(real(noise .* exp(-angle(modData))))
  end
  output_FKMPSK_detrend = detrend(output_FKMPSK,'linear');
  
  temp = mean(output_FKMPSK);
  TN1 = std(output_FKMPSK_detrend)/temp;
  if TN1 > T
      Num_N1(j2) = Num_N1(j2) + 1;
  end
  
  output_MPSKN0_detrend = detrend(output_MPSKN0,'linear');
  temp = mean(output_MPSKN0);
  TN0 = std(output_MPSKN0_detrend)/temp;
  if TN0 > T
      Num_N0(j2) = Num_N0(j2) + 1;
  end
  if Num_N0(j2) == Max_Num_errors
      Total_trials(j2) = j1
      save('201ChirpRoC10000D10dBQAMtemp','Num_N0','Num_N1','Total_trials','Vec_T');
      break;
  end
 end
 if j1 == Max_trials
       save('201ChirpRoC10000D10dBQAMtemp','Num_N0','Num_N1','Total_trials','Vec_T');
     break;
 end
end


%{
figure
plot(x_time, output_FKMPSK, '-b', x_time, output_MPSKN0, '-k');
legend('N1','N0');
title('signal with modulation and noise');
%}

%{
figure
output_FKMPSK_detrend = detrend(output_FKMPSK,'linear');
output_MPSKN0_detrend = detrend(output_MPSKN0,'linear');
plot(x_time,output_FKMPSK_detrend,'-b',x_time, output_MPSKN0_detrend,'-k');
legend('N1','N0');
title('after LEO compensation');
%}



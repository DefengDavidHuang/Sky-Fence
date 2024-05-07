clear all

%load('290offsetinX65degree250mpssignal.mat');
%load('301.9offsetinX65degree250mpssignal.mat');
%load('290.1offsetinX65degree250mpssignal.mat');
%load('280offsetinX65degree250mpssignal.mat');
%load('240offsetinX65degree250mpssignal.mat');
%load('200offsetinX65degree250mpssignal.mat');
%load('181offsetinX65degree250mpssignal.mat');
%load('101offsetinX65degree250mpssignal.mat');
%load('51offsetinX65degree250mpssignal.mat');
%INPUT from the file: outputgain_N1, outputgain_N0

%load('301.9offsetinX65degree250mpssignal5second.mat');
%load('301.9offsetinX65degree250mpssignal2second.mat');
%load('301.9offsetinX65degree250mpssignal05second.mat');
%load('301.9offsetinX65degree250mpssignal025second.mat');
%load('301.9offsetinX65degree250mpssignal01second.mat');
%load('301.9offsetinX65degree150mpssignal8second.mat')
%load('301.9offsetinX65degree0mpssignal06second.mat')
load('301.9offsetinX65degree-250mpssignal03second.mat')

Num_data_point = length(x_time);

%{
figure
plot(x_time, abs(outputgain_N1), 'b',x_time,abs(outputgain_N0),'.k');
legend('Fresnel Kirchhoff');
xlabel('time (second)');
ylabel('received signal (linear) modified');

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

%{
c = [-3-3i -3-1i -3+1i -3+3i -1-3i -1-1i -1+1i -1+3i 1-3i 1-1i 1+1i 1+3i 3-3i 3-1i 3+1i 3+3i]; % /sqrt(9.75);  %standard 16QAM
M = length(c);
temp = sqrt(c*c'/M);
c = c/temp; %16QAM, power normalize to 1
%}

c = [-1-1i -1+1i 1-1i 1+1i]/sqrt(2);
M = length(c);

output_FKMPSK = y_FK;
output_MPSKN0 = y_FK;
SNR = 10; %20; %10; %10; %3; %10; %10; %in dB
data_length = 10000; %1000; %5000; %10000; %5000; %1000; %10000; %5000; %1000; %10000; %5000; %150; %10; %50; %10; %100; %100; %10; %100; %1000; %5000; %1000; %10000; %5000; %1000; %10000; %1000; %5000; %1000; %1000; % 5000; %100; %10000; %1000; %10000;

%Vec_T = 0.012:0.0001:0.016; %10000symbols
%Vec_T = 0.0010:0.0002:0.0090;   %1000 symbols, 20dB
%Vec_T = 0.0030:0.0001:0.0050;   %5000 symbols
Vec_T = 0.0022:0.00002:0.0045;
%Vec_T = 0.015:0.0002:0.025; %0.06:0.002:0.09;
%Vec_T = 0.025:0.001:0.035;
%Vec_T = 0.0080:0.0002:0.0110; %1000 symbols

%Vec_T = 0.0020:0.0001:0.0035;
%Vec_T = 0.02:0.001:0.035; %0.002: 0.0001:0.005; %0.003: 0.0001:0.005; %0.006:0.0002:0.01; %0.0086;

%Vec_T = 0.0020:0.0001:0.015; %7; %10000 symbols 16QAM
%Vec_T = 0.01:0.001:0.016; %7; %5000 symbols 16QAM
%Vec_T = [0.0055:0.0004:0.007, 0.0071:0.0001:0.0090]; %1000 symbols 16QAM
Num_N1 = zeros(length(Vec_T),1);
Num_N0 = Num_N1;

Max_trials = 1e6;
Total_trials = Max_trials * ones(length(Vec_T),1);

Max_Num_errors = 100; %100; %5; %10; %0;

%The variance of the amplitude
c_amplitude = abs(c);

var_amplitude = var(c_amplitude,1) % (M-1)/M;
A = mean(c_amplitude,'all');
%var(c_amplitude,1)
var_noise = 0.5* 10^(-SNR/10); %0.7737 * 10^(-SNR/10); %0.5* 10^(-SNR/10); % %%0.7737 * 10^(-SNR/10); % 10^(-SNR/10); %0.5 * 10^(-SNR/10); % 0.7737 * 10^(-SNR/10)   %variance reduced due to digital modulation
var_total = (var_amplitude * (outputgain_N0 .* conj(outputgain_N0)) + var_noise)/data_length;
%{
var_noiseN0 = 10^(-SNR/10);
var_totalN0 = (var_amplitude * (outputgain_N0 .* conj(outputgain_N0)) + var_noiseN0)/data_length
%}
%sqrt(var_total)

for j2 = 1:1:length(Vec_T)
 T = Vec_T(j2)
 %Num_N1
 %Num_N0
 for j1 = 1:1:Max_trials
 
 %{
 for m = 1:1:Num_data_point
     %data = randi([0 M-1],data_length,1);
     %modData = genqammod(data,c);
     RXgain = outputgain_N1(m);
     rxSig = awgn(modData*RXgain,SNR,'measured','db');
     output_FKMPSK(m) = mean(abs(rxSig));
     
     RXgainN0 =  outputgain_N0(m);
     rxSigN0 = awgn(modData*RXgainN0,SNR,'measured','db');
     output_MPSKN0(m) = mean(abs(rxSigN0));
  end
  %}
     AWGNnoise = sqrt(var_total) .* randn(Num_data_point,1);
    % output_FKMPSK = abs(outputgain_N1) + AWGNnoise;
    output_FKMPSK = A * abs(outputgain_N0) + A * real(outputgain_target) + AWGNnoise;
     output_FKMPSK_detrend = detrend(output_FKMPSK,'linear');
  
  temp = mean(output_FKMPSK);
  TN1 = std(output_FKMPSK_detrend)/temp;
  if TN1 > T
      Num_N1(j2) = Num_N1(j2) + 1;
  end
  
  output_MPSKN0 = A * abs(outputgain_N0) + AWGNnoise; % sqrt(var_total) * randn(Num_data_point,1);
  
  %{
  AWGNnoiseN0 = sqrt(var_totalN0) .* randn(Num_data_point,1);
  output_MPSKN0 = A * abs(outputgain_N0) + AWGNnoiseN0;
  %}
  
  output_MPSKN0_detrend = detrend(output_MPSKN0,'linear');
  temp = mean(output_MPSKN0);
  TN0 = std(output_MPSKN0_detrend)/temp;
  if TN0 > T
      Num_N0(j2) = Num_N0(j2) + 1;
  end
  if Num_N0(j2) == Max_Num_errors
      Total_trials(j2) = j1
      save('301ChirpRoCnum10000D10dBQAMtemp03second-250mps','Num_N0','Num_N1','Total_trials','Vec_T');
      break;
  end
 end
 if j1 == Max_trials
     save('301ChirpRoCnum10000D10dBQAMtemp03second-250mps','Num_N0','Num_N1','Total_trials','Vec_T');
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



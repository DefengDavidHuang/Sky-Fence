clear all

%load('205.5offsetinX65degree250mpssignal.mat');
%load('25.9offsetinX65degree250mpssignal.mat');
%load('105.5offsetinX65degree250mpssignal.mat');

%load('120offsetinX65degree250mpssignal.mat');
%load('301.9offsetinX65degree250mpssignal.mat');
%load('240offsetinX65degree250mpssignal.mat');

load('20.1offsetinX65degree100.1mpssignal02secondlower.mat');
%INPUT from the file: outputgain_N1, outputgain_N0

plot(x_time', abs(outputgain_N1),'-k');

Num_data_point = length(x_time);
directpowr = norm(output_N0)^2/length(output_N0);
%directpowrDB = 10 * log10(directpowr)
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
c = [-3-3i -3-1i -3+1i -3+3i -1-3i -1-1i -1+1i -1+3i 1-3i 1-1i 1+1i 1+3i 3-3i 3-1i 3+1i 3+3i]; % /sqrt(9.75);  %standard 16QAM
%c = [-5 -5i 5 5i -3 -3-3i -3i 3-3i 3 3+3i 3i -3+3i -1 -1i 1 1i]; % 16-QAM constellation based on the V.29 standard for telephone-line modems.
%c = [-3-3i -3-1i -3+1i -3+3i -1-3i -1-1i -1+1i -1+3i 1-3i 1-1i 1+1i 1+3i 3-3i 3-1i 3+1i 3+3i]/sqrt(9.75);  %standard 16QAM
%c = [-1-i -1-1i -1+1i -1+1i -1-i -1-1i -1+1i -1+1i 1-1i 1-1i 1+1i 1+1i 1-1i 1-1i 1+1i 1+1i]/sqrt(2);

%c = [1 1i -1 -1i];
%M = length(c);
M = length(c);
temp = sqrt(c*c'/M);
c = c/temp;

output_FKMPSK = zeros(Num_data_point,1); %y_FK;
output_MPSKN0 = output_FKMPSK; %y_FK;
DNR_dB = 10; %8; %10; %3; %10; %10; %in dB
DNR = 10^(DNR_dB/10);

data_length = 200; %100; %1e3; %1000; %10000; %1000; %10000;

Vec_T = 0.1100:0.002:0.170; %570; %0.09:0.001:0.12; %0.021:0.0001:0.022; %0.005:0.001:0.017; %0.0035:0.0002:0.0075; %0.01:0.0003:0.02; %0.015; %0.02; %0.015; %0.004:0.0001:0.005; %0.004
Num_N1 = zeros(length(Vec_T),1);
Num_N0 = Num_N1;

Max_trials = 1e6; %1e5; %1e6; %1e5; %1e6;
Total_trials = Max_trials * ones(length(Vec_T),1);

Max_Num_errors = 100; %0; %100; %5; %10; %0;

for j2 = 1:1:length(Vec_T)
 T = Vec_T(j2)
 %Num_N1
 %Num_N0
 for j1 = 1:1:Max_trials
  for m = 1:1:Num_data_point
     data = randi([0 M-1],data_length,1);
     modData = genqammod(data,c);
     %RXgain = ;
      rxSig = sqrt(2*DNR/directpowr)*modData*outputgain_N1(m) +  randn(data_length,1) + 1j * randn(data_length,1); %awgn(modData*outputgain_N1(m),SNR); %,'measured','db');
    % rxSig = awgn(modData*RXgain,SNR,'measured','db');
     output_FKMPSK(m) = mean(abs(rxSig));
     
    % RXgainN0 =  ;
    rxSigN0 = sqrt(2*DNR/directpowr)*modData*outputgain_N0(m) +  randn(data_length,1) + 1j * randn(data_length,1);
    % rxSigN0 = awgn(modData*outputgain_N0(m),SNR); %,'measured','db');
     output_MPSKN0(m) = mean(abs(rxSigN0));
  end
  output_FKMPSK_detrend = detrend(output_FKMPSK,'linear');
  
  %temp = mean(output_FKMPSK);
  TN1 = std(output_FKMPSK_detrend); %/temp;
  if TN1 > T
      Num_N1(j2) = Num_N1(j2) + 1;
  end
  
  output_MPSKN0_detrend = detrend(output_MPSKN0,'linear');
  %temp = mean(output_MPSKN0);
  TN0 = std(output_MPSKN0_detrend); %/temp;
  if TN0 > T
      Num_N0(j2) = Num_N0(j2) + 1;
  end
  if Num_N0(j2) == Max_Num_errors
      Total_trials(j2) = j1
      Num_N0
      Num_N1
      save('ChirpRoCupdate10dB16QAMNA200','Num_N0','Num_N1','Total_trials','Vec_T');
      break;
  end
 end
 if j1 == Max_trials
     Num_N0
     Num_N1
     break;
 end
end

save('MPSKChirpRoCtest1update10dB16QAMNA200','Num_N0','Num_N1','Total_trials','Vec_T');

semilogx(Num_N0./Total_trials, Num_N1./Total_trials, 'k.-');
xlabel('Prob of false alarm');
ylabel('Probability of detection (P_{d})');
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


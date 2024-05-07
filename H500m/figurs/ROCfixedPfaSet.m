% Signal template parameters
%load('301.1offsetinX65degree100.1mpssignal02second.mat');

%load('20.1offsetinX90degree100.1mpssignal02second.mat');
%load('200.1offsetinX90degree100.1mpssignal02second.mat');

%load('200.1offsetinX40degree100.1mpssignal02second.mat');
%load('200.1offsetinX40degree100.1mpssignal02secondlower.mat');
%load('200.1offsetinX40degree100.1mpssignal02secondlower.mat');

%load('0.3seconds65degree100mpssignal02second.mat');
%load('0.3seconds65degree13mpssignal02seconddrone.mat');
%load('0.5seconds65degree13mpssignal02seconddrone.mat');

%load('0.5seconds65degree100.1mpssignal02seconddroneset.mat');
%load('0.5seconds65degree100.1mpssignal02seconddroneset250mintotal.mat');
%load('0.5seconds65degree20mpssignal02seconddroneset250mintotalDJI.mat');
%load('0.5seconds90degree20mpssignal02seconddroneset250mintotalDJI.mat');
%load('0.5seconds65degree20mpssignal02seconddroneset250mintotalDJImini.mat');
%load('0.5seconds65degree20mpssignal02seconddroneset250mintotalDJI500m.mat');
%load('0.5seconds65degree20mpssignal02seconddroneset250mintotalDJI500mlarge.mat'); %0.66m x 0.9m
%load('0.4seconds65degree20mpssignal02seconddroneset250mintotalDJI500mlargeupdate1times.mat'); %0.33m x 0.45m, 0.4second
%load('0.4seconds65degree20mpssignal02seconddroneset250mintotalDJI500mlargeupdate064times081.mat'); %0.64m x 0.81m;
%load('0.4seconds65degree20mpssignal02seconddroneset250mintotalDJI100msmallupdate013times016.mat'); %100m height, 0.13mx0.16m
%load('0.4seconds65degree-20mpssignal02seconddroneset250mintotalDJI100msmallupdate013times016.mat'); %100m height, 0.13mx0.16m, -20m/s
load('0.4seconds65degree-20mpssignal02seconddroneset250mintotalDJI100msmallupdate064times081.mat'); %500m height, 0.64mx0.81m,-20m/s

output_N0 = output_N0_set(:,1);
directpowr = norm(output_N0)^2/length(output_N0);
directpowrDB = 10 * log10(directpowr)

signal_O = detrend(output_signalset); %,'linear');
%signal_template = abs(outputgain_N0); %... % Define your signal template here
%signalPower = norm(signal_O)^2 %/length(signal_template)  %; % Signal power

data_length = 1e5; %1e5; %2e4; %1e4; %1e4; %1e4; %1e4; %1e5; %1e6; %1e5; %1e4; %100000; %1e4; %10000; %10000;

%c = [-1-1i -1+1i 1-1i 1+1i];
c = [-3-3i -3-1i -3+1i -3+3i -1-3i -1-1i -1+1i -1+3i 1-3i 1-1i 1+1i 1+3i 3-3i 3-1i 3+1i 3+3i]; % /sqrt(9.75);  %standard 16QAM
M = length(c);
temp = sqrt(c*c'/M);
c = c/temp; %16QAM, power normalize to 1

%c = [-1-1i -1+1i 1-1i 1+1i]/sqrt(2);
%M = length(c);

%output_FKMPSK = y_FK;
%output_MPSKN0 = y_FK;

%The variance of the amplitude
c_amplitude = abs(c);

var_amplitude = var(c_amplitude,1) % (M-1)/M;
A = mean(c_amplitude,'all');
%var(c_amplitude,1)
%SNR_dB = -15:5.1:15; %15; %20; %in dB, this SNR is the SNR at the receiver, not relevant to the signal produced by the flying object
DNR_dB = 10;
SNR_dB = DNR_dB - directpowrDB; %9.8; %-5:5.1:15;
SNR = 10.^(SNR_dB./10); % SNR in linar form
var_noise = 0.5 ./ SNR ; %0.7737 * 10^(-SNR/10); %0.5* 10^(-SNR/10); % %%0.7737 * 10^(-SNR/10); % 10^(-SNR/10); %0.5 * 10^(-SNR/10); % 0.7737 * 10^(-SNR/10)   %variance reduced due to digital modulation


 %Ot = detrend(y_FK,'linear');
  PFAtemp = -9:1:-6; %-6; %:0.3:-1;
  PFA = 10.^PFAtemp;
 index1 = size(signal_O);
 PD_all = zeros(index1(2), length(PFA));
 PD_matched = PD_all;
 
  N_D = length(signal_O(:,1));
  
 for i = 1:index1(2)
     Noisepowerduetosignal = mean( var_amplitude * (output_N0 .* conj(output_N0)));
    var_total = (Noisepowerduetosignal + var_noise)/data_length;
   lamda = signal_O(:,i)' * signal_O(:,i) ./var_total;
  
% Probability of false alarm (PFA)

%PFA = 0:0.01:1;

% Initialize arrays to store PD and PFA for each SNR value

%PFA_all = zeros(length(PFA), length(SNR));



% Compute PD and PFA for each SNR value
    for i1 = 1: length(PFA)
    
     PD_all(i,i1) = 1 - ncx2cdf(chi2inv(1 - PFA(i1),N_D), N_D, lamda);
    %temp = 
     PD_matched(i,i1) = qfunc(qfuncinv(PFA(i1)) - sqrt(lamda));
    %PD_all(:,i) = qfunc(qfuncinv(PFA) - sqrt(signalPower / var_total(i)));
    %PFA_all(:,i) = PFA;
    end
end

% Plot ROC Curves for each SNR value
%{
figure;
plot(SNR_dB + directpowrDB, PD_all, 'k.-',SNR_dB + directpowrDB, PD_matched, 'r-');
xlabel('DNR (dB)');
ylabel('Probability of detection (P_{d})');
legend(cellstr(num2str(PFA', 'Prob of false alarm = %0.7f')));
title('Prob of Detection for HueristicDetector (Theoretical)');
%}
%x_initial_set = 0:10.01:450;

%semilogx(x_initial_set, PD_all, 'k.-');
%semilogx(x_initial_set, PD_matched, 'k.-');
semilogx(x_initial_set, PD_all(:,1), 'kx-',x_initial_set, PD_matched(:,1), 'ks:');

%semilogx(x_initial_set, PD_all(:,1), 'k.-', x_initial_set, PD_all(:,2), 'r-',x_initial_set, PD_all(:,3), 'kx-',x_initial_set, PD_all(:,4), 'm.-');
%semilogx(x_initial_set, PD_matched(:,1), 'k.-', x_initial_set, PD_matched(:,2), 'r-',x_initial_set, PD_matched(:,3), 'kx-',x_initial_set, PD_matched(:,4), 'm.-');
xlabel('Minimum distance between target and communication link (metres)');
ylabel('Probability of detection (P_{D})');
%legend(cellstr(num2str((PFAtemp)', 'PFAtemp = %0.1f')));
legend('Heuristic Detector','CVD')
title('Hueristic Detector and CVD, P_{fa}=10^{-9}');

%grid on;
%title('ROC Curve for Matched Filter (Theoretical)');

%semilogx(PFA, PD_all, 'b.-');
%xlabel('Probability of false alarm (P_{fa})');
%ylabel('Probability of detection (P_{d})');



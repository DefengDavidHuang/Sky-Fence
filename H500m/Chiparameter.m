% Signal template parameters
%load('301.1offsetinX65degree100.1mpssignal02second.mat');

%load('20.1offsetinX90degree100.1mpssignal02second.mat');
%load('200.1offsetinX90degree100.1mpssignal02second.mat');

load('200.1offsetinX40degree100.1mpssignal02second.mat');
%load('200.1offsetinX40degree100.1mpssignal02secondlower.mat');

directpowr = norm(output_N0)^2/length(output_N0);
directpowrDB = 10 * log10(directpowr)
signal_O = detrend(y_FK,'linear');
%signal_template = abs(outputgain_N0); %... % Define your signal template here
%signalPower = norm(signal_O)^2 %/length(signal_template)  %; % Signal power

data_length = 1e5; %1e4; %100000; %1e4; %10000; %10000;

c = [-1-1i -1+1i 1-1i 1+1i];
%c = [-3-3i -3-1i -3+1i -3+3i -1-3i -1-1i -1+1i -1+3i 1-3i 1-1i 1+1i 1+3i 3-3i 3-1i 3+1i 3+3i]; % /sqrt(9.75);  %standard 16QAM
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
SNR_dB = -6:1:10; %20; %in dB, this SNR is the SNR at the receiver, not relevant to the signal produced by the flying object
SNR = 10.^(SNR_dB./10); % SNR in linar form
var_noise = 0.5 ./ SNR ; %0.7737 * 10^(-SNR/10); %0.5* 10^(-SNR/10); % %%0.7737 * 10^(-SNR/10); % 10^(-SNR/10); %0.5 * 10^(-SNR/10); % 0.7737 * 10^(-SNR/10)   %variance reduced due to digital modulation
Noisepowerduetosignal = mean( var_amplitude * (outputgain_N0 .* conj(outputgain_N0)));
var_total = (Noisepowerduetosignal + var_noise)/data_length;

 %Ot = detrend(y_FK,'linear');
 lamda = signal_O' * signal_O ./var_total(10);
 N_D = 1:10:10000; %[1,10,100,1000,10000]; %length(signal_O);
%Noisetotal 
% SNR values
 % Range of SNR values in dB
%SNR = 10.^(SNR_dB/10); % Convert dB to linear scale

% Probability of false alarm (PFA)
temp = -7:0.5:-6.5;
PFA = 10.^temp;
%PFA = 0:0.01:1;

% Initialize arrays to store PD and PFA for each SNR value
PD_all = zeros(length(PFA), length(N_D));
PFA_all = zeros(length(PFA), length(N_D));

% Compute PD and PFA for each SNR value
for i = 1:length(N_D)
    
    PD_all(:,i) = 1 - ncx2cdf(chi2inv(1 - PFA,N_D(i)), N_D(i), lamda);
    %PD_all(:,i) = qfunc(qfuncinv(PFA) - sqrt(signalPower / var_total(i)));
    %PFA_all(:,i) = PFA;
end

% Plot ROC Curves for each SNR value
figure;
plot(N_D, PD_all, 'k.-');
xlabel('Degree of freedom');
ylabel('Probability of detection (P_{d})');
%title('ROC Curve for Matched Filter (Theoretical)');

%semilogx(PFA, PD_all, 'b.-');
%xlabel('Probability of false alarm (P_{fa})');
%ylabel('Probability of detection (P_{d})');
title('Prob of Detection for HueristicDetector (Theoretical)');
%grid on;
legend(cellstr(num2str(PFA', 'Prob of false alarm = %0.7f')));

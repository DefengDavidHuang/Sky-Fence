% Signal template parameters
%load('301.1offsetinX65degree100.1mpssignal02second.mat');

%load('20.1offsetinX90degree100.1mpssignal02second.mat');
%load('200.1offsetinX90degree100.1mpssignal02second.mat');

%load('200.1offsetinX40degree100.1mpssignal02second.mat');
%load('200.1offsetinX40degree100.1mpssignal02secondlower.mat');
%load('200.1offsetinX40degree100.1mpssignal02secondlower.mat');
%load('20.1offsetinX65degree100.1mpssignal02secondlower.mat');
%load('20.1offsetinX65degree100.1mpssignal02secondlower.mat');

load('20.1offsetinX65degree100.1mpssignal02secondlower20Mupdate.mat');

%20.1offsetinX65degree100.1mpssignal02secondlower20M
%{
plane_height = 500.1; %101.1; %500; %1.0e4;  %1.0e3; %
plane_speed = 100.1; %300.1; %100.1; %10; %30; %plane speed in m/s
width_z = 6.4; %1; %2.1; %6.1; %length of plane
length_x = 8.1; %3.1; %8.1; %width of plane  %x: parallel to ground and perpendicular to the satellite link
x_track_initial = 20.1; %0.0; %100.1; %200.1; %0.0; %150.1; %111.1; %5; % 101.1; %in the x direction, the offset of the center of the airplane to the plane scanned by the satellite-reciever link 
%}

directpowr = norm(output_N0)^2/length(output_N0);
directpowrDB = 10 * log10(directpowr)
signal_O_temp = detrend(y_FK,'linear');
%signal_template = abs(outputgain_N0); %... % Define your signal template here
%signalPower = norm(signal_O)^2 %/length(signal_template)  %; % Signal power

%data_length = 100; %1e3; %1e5; %1e4; %100000; %1e4; %10000; %10000;
NA = [10 100 1000 10000 1e5]; %[1 10 20 50];  % 1000 10000 100000];
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
%A = mean(c_amplitude,'all');

%var(c_amplitude,1)
DNR_dB = -5; %8; 
SNR_dB = DNR_dB - directpowrDB; %10.8; %5; %11; %10; %0:1:5; %15; %20; %in dB, this SNR is the SNR at the receiver, not relevant to the signal produced by the flying object
SNR = 10.^(SNR_dB./10); % SNR in linar form
var_noise = 0.5 ./ SNR ; %0.7737 * 10^(-SNR/10); %0.5* 10^(-SNR/10); % %%0.7737 * 10^(-SNR/10); % 10^(-SNR/10); %0.5 * 10^(-SNR/10); % 0.7737 * 10^(-SNR/10)   %variance reduced due to digital modulation

%var_total = (Noisepowerduetosignal + var_noise)/data_length;
 %Ot = detrend(y_FK,'linear');
 
 

%Noisetotal 
% SNR values
 % Range of SNR values in dB
%SNR = 10.^(SNR_dB/10); % Convert dB to linear scale

% Probability of false alarm (PFA)
temp = -10:0.5:0; %-1;
PFA = 10.^temp;
%PFA = 0:0.01:1;

% Initialize arrays to store PD and PFA for each SNR value
PD_all = zeros(length(PFA), length(NA));
%PFA_all = zeros(length(PFA), length(SNR));
PD_matched = PD_all;


% Compute PD and PFA for each SNR value
for i = 1:length(NA)
    data_length = NA(i);
    %indextemp = 1:data_length/10:length(outputgain_N0);
    indextemp = 1:data_length:length(outputgain_N0)-1;
    signal_O = signal_O_temp(indextemp);
     lamda = signal_O' * signal_O; % ./var_total;
     N_D = length(signal_O);
    outputgain_N0_temp = outputgain_N0(indextemp);
    Noisepowerduetosignal = var_amplitude * mean(outputgain_N0_temp .* conj(outputgain_N0_temp));
    var_total = (Noisepowerduetosignal + var_noise)/data_length;
    PD_all(:,i) = 1 - ncx2cdf(chi2inv(1 - PFA,N_D), N_D, lamda/var_total);
    %data_length * lamda/var_total
    %temp = 
    PD_matched(:,i) = qfunc(qfuncinv(PFA) - sqrt(lamda/var_total));
    %PD_all(:,i) = qfunc(qfuncinv(PFA) - sqrt(signalPower / var_total(i)));
    %PFA_all(:,i) = PFA;
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

semilogx(PFA', PD_all(:,1), 'k.:',...
    PFA', PD_all(:,2), 'ko:',...
    PFA', PD_all(:,3), 'kx:',...
    PFA', PD_all(:,4), 'k*:', ...
    PFA', PD_all(:,5), 'ks:',...
    PFA', PD_matched,'r-');
xlabel('Prob of false alarm (P_{fa})');
ylabel('Prob of detection (P_{D})');
legend('$N_A$=10','$N_A$=100', '$N_A$=1,000', '$N_A$=10,000', '$N_A$=100,000', 'CVD', 'Interpreter','latex');
title('Receiver Operating Characteristic (DNR = -5dB)');

save('Theoretical10dB16QAMNAfinal.mat','PFA','PD_all')
%grid on;
%title('ROC Curve for Matched Filter (Theoretical)');

%semilogx(PFA, PD_all, 'b.-');
%xlabel('Probability of false alarm (P_{fa})');
%ylabel('Probability of detection (P_{d})');




%load('MPSKChirpRoCtest1.mat'); %,'Num_N0','Num_N1','Total_trials','Vec_T');
load('MPSKChirpRoCtest1update.mat');
PFA10dBsim = Num_N0./Total_trials;
PD10dBsim = Num_N1./Total_trials;

load('Theoretical10.8dB.mat');
PFA10dB = PFA;
PD_all10dB = PD_all;

load('Theoretical8dB.mat'); %,'PFA','PD_all')
PFA8dB = PFA;
PD_all8dB = PD_all;


load('MPSKChirpRoCtest1update8dB.mat');
PFA8dBsim = Num_N0./Total_trials; %PFA;
PD8dBsim = Num_N1./Total_trials;
%PD_all;

load('Theoretical8dB16QAM.mat'); %,'PFA','PD_all')
PFA8dBQAM = PFA;
PD_all8dBQAM = PD_all;

load('MPSKChirpRoCtest1update8dB16QAM.mat');
%MPSKChirpRoCtest1update8dB16QAM
PFA8dBsimQAM = Num_N0./Total_trials; %PFA;
PD8dBsimQAM = Num_N1./Total_trials;

load('Theoretical10dB16QAM.mat'); %,'PFA','PD_all')
PFA10dBQAM = PFA;
PD_all10dBQAM = PD_all;

load('MPSKChirpRoCtest1update10dB16QAM.mat');
PFA16dBsimQAM = Num_N0./Total_trials; %PFA;
PD16dBsimQAM = Num_N1./Total_trials;

load('Theoretical10dB16QAMNA1000.mat');
PFA10dBQAMNA200 = PFA;
PD_all10dBQAMNA200 = PD_all;

load('MPSKChirpRoCtest1update10dB16QAMNA200.mat');
PFA110dBsimQAMNA200 = Num_N0./Total_trials; %PFA;
PD110dBsimQAMNA200 = Num_N1./Total_trials;

semilogx(PFA10dBsim, PD10dBsim, 'k:', PFA10dB, PD_all10dB, 'k-', ...
    PFA10dBsim, PD10dBsim, 'ks', ...
    PFA8dB, PD_all8dB, 'kv', ...
    PFA10dBQAM, PD_all10dBQAM, 'ko', ...
    PFA8dBQAM, PD_all8dBQAM, 'kx', ...
    PFA10dBQAMNA200, PD_all10dBQAMNA200, 'k+', ...
     PFA10dBsim, PD10dBsim, 'ks:', PFA10dB, PD_all10dB, 'ks-',...
     PFA8dB, PD_all8dB, 'kv-', PFA8dBsim, PD8dBsim, 'kv:',...
      PFA10dBQAM, PD_all10dBQAM, 'ko-', PFA16dBsimQAM, PD16dBsimQAM, 'ko:',...
       PFA8dBQAM, PD_all8dBQAM, 'kx-',  PFA8dBsimQAM, PD8dBsimQAM, 'kx:',...
        PFA10dBQAMNA200, PD_all10dBQAMNA200, 'k+-', PFA110dBsimQAMNA200, PD110dBsimQAMNA200, 'k+:');

xlabel('Probability of false alarm');
ylabel('Probability of detection');
legend('Simulation','Numerical','DNR = 10dB, QPSK','DNR = 8dB, QPSK','DNR=10dB, 16QAM','DNR=8dB,16QAM', 'DNR=10dB, 16QAM, N_A=200');  %'theoretical, 10.5dB',
axis([1e-4, 1, 0, 1])
clear all

%f = 0:10:1e6;

B = 1e8;
temp = 0.5:0.01:10; %[0.5:0.1:10];
Ts = temp/B;
Gnef0 = 2 * (sinint(2 * pi * B *Ts) - (sin(pi * B * Ts)).^2./(pi * B * Ts))/pi;
plot(temp, Gnef0,'-k');
xlabel('BTs');
ylabel('PSD at f=0');



f = -10000:10.1:10000;
%f = -1e8:100000.1:1e8; %0:(1e7+1):10*B;
Ts = 0.5/B; 
G1ne2 = (sinint(2 * pi * (f + B)*Ts) - sinint(2 * pi * (f - B)*Ts) + (sin(pi * (f - B) * Ts)).^2./(pi * (f - B) * Ts) - (sin(pi * (f + B) * Ts)).^2./(pi * (f + B) * Ts))/pi;

Ts = 0.75/B;
%f = -1e4:100:1e4;

G3ne4 = (sinint(2 * pi * (f + B)*Ts) - sinint(2 * pi * (f - B)*Ts) + (sin(pi * (f - B) * Ts)).^2./(pi * (f - B) * Ts) - (sin(pi * (f + B) * Ts)).^2./(pi * (f + B) * Ts))/pi;

Ts = 1/B;
G1ne1 = (sinint(2 * pi * (f + B)*Ts) - sinint(2 * pi * (f - B)*Ts) + (sin(pi * (f - B) * Ts)).^2./(pi * (f - B) * Ts) - (sin(pi * (f + B) * Ts)).^2./(pi * (f + B) * Ts))/pi;

Ts = 2/B;
G2ne1 = (sinint(2 * pi * (f + B)*Ts) - sinint(2 * pi * (f - B)*Ts) + (sin(pi * (f - B) * Ts)).^2./(pi * (f - B) * Ts) - (sin(pi * (f + B) * Ts)).^2./(pi * (f + B) * Ts))/pi;

Ts = 20/B;
G20ne1 = (sinint(2 * pi * (f + B)*Ts) - sinint(2 * pi * (f - B)*Ts) + (sin(pi * (f - B) * Ts)).^2./(pi * (f - B) * Ts) - (sin(pi * (f + B) * Ts)).^2./(pi * (f + B) * Ts))/pi;

figure
plot(f, G1ne2, '-k', f, G3ne4, '-r', f, G1ne1, '-b', f, G2ne1, '-xk', f, G20ne1, '-.k');
legend('1/2','3/4','1', '2','20');
xlabel('frequency in Hz');
ylabel('PSD');
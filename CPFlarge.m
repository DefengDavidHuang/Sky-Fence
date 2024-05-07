clear all

std_noise = sqrt(0.01);

SNR_dB_set = 0:5:40

step = 0.005; %1;
xmax = 1;
x1 = -xmax:step:0; %xmax;
x2 = 0:step:xmax;

[X1, X2] = meshgrid(x1,x2);

Y1 = fresnelc(X2) - fresnelc(X1);
Y2 = fresnels(X2) - fresnels(X1);

total = 1000;

zlower = -0.5; %-1; %-0.5; %-0.27; %-2.5, -0.27; %-0.027
zupper = 0.5; %1; %0.5; %0.15; %2, 0.15; %0.015
leng = zupper - zlower

y1lower1 = fresnelc(zupper) - fresnelc(zlower);
y2upper1 = fresnels(zupper) - fresnels(zlower);

mse1set = zeros(length(SNR_dB_set),1);
mse2set = mse1set;
mseset = mse1set;
msebenchset = mse1set;

for i2 = 1:length(SNR_dB_set)
 SNR_dB = SNR_dB_set(i2);
 mse1 = 0;
 mse2 = 0;
 mse = 0;
 msebenchmark = 0;
% total1 = 0;
 for monto = 1:1:total
    %monto
   y1lower = awgn(y1lower1,SNR_dB,'measured');
   y2upper = awgn(y2upper1,SNR_dB,'measured');

   tempa1 = (Y1 - y1lower).*(Y1 - y1lower);
   tempa2 = (Y2 - y2upper).*(Y2 - y2upper);
   A = tempa1 + tempa2;
   minimum = min(min(A));
   [indx1,indx2] = find(A==minimum);

   output1 = mean(x2(indx1));
   output2 = mean(x1(indx2));
   lengthoutput = output1 - output2;
   
  % if length(indx1) == 1
     mse1 = mse1 + (output1 - zlower) * (output1 - zlower);
     mse2 = mse2 + (output2 - zupper) * (output2 - zupper);
     mse = mse + (lengthoutput - leng) * (lengthoutput - leng);
     
     msebenchmark = msebenchmark + (sqrt(y2upper * y2upper + y1lower * y1lower) - leng)^2;
   %  total1 = total1 + 1;
   %end
 end

mse1set(i2) = mse1/total; %/leng/leng;
mse2set(i2) = mse2/total; %/leng/leng;
mseset(i2)  = mse/total; %/leng/leng;
msebenchset(i2) = msebenchmark/total; %/leng/leng;
end

figure
plot(SNR_dB_set, mse1set,'-x',SNR_dB_set, mse2set,'-o',SNR_dB_set, mseset,'-b', SNR_dB_set, msebenchset,'--r'); 
legend('lower bound','upper bound','new method','benchmark')

xlabel('SNR in dB');
ylabel('MSE');
%figure
%surf(x1,x2,Y1);

%hold on
%figure
%surf(x1,x2,Y2);

%Y1R = Y1 + std_noise * randn(size(Y1));
%Y2R = Y2 + std_noise * randn(size(Y2));


%{
figure
plot(xome', kapp1k, '-k');
figure
plot(xome', kapp2k, '-xk');
%}


      
%{
      temp = sin(((kapp1k .* kapp1k)-(kapp2k .* kapp2k))*pi/2);
      
      kapp1kupdate = kapp1k - (sin(pi* (kapp2k .* kapp2k)/2) .* temp1  - (cos(pi* (kapp2k .* kapp2k)/2) .* temp2)) ./temp;
      kapp2kupdate = kapp2k - (sin(pi* (kapp1k .* kapp1k)/2) .* temp1 - (cos(pi* (kapp1k .* kapp1k)/2) .* temp2)) ./temp;
      
      kapp1k = kapp1kupdate;
      kapp2k = kapp2kupdate;
     % plot(xome', kapp1k, '-k');  
  
 %}
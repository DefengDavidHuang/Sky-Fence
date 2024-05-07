alpha1 = 1:0.1:10;
alpha = alpha1/180*pi;
theta1 = 65; %45; %65;
theta = theta1/180*pi;
hvec = [100 200 500 1000]; %100:100:1000;

L1 = zeros(length(alpha),length(hvec));
L2 = L1;
index = 0;
for h = hvec
    index = index + 1;
    L1(:, index) = h/2 * (1./tan(theta - alpha/2) - 1./tan(theta + alpha/2));
    L2(:,index) = h * alpha/2/sin(theta)/sin(theta);
end

%plot(alpha1, L1, 'kx-', alpha1, L2, 'r.:');
plot(alpha1, L1(:,1), 'k-', alpha1, L1(:,2), 'k:',alpha1, L1(:,3), 'k-.',alpha1, L1(:,4), 'k--');
xlabel('Receive antenna beamwidth (in degree)');
ylabel('Coverage range (in metres)');
legend('target height = 100m','target height = 200m','target height = 500m','target height = 1000m');
%legend(cellstr(num2str(hvec', 'Object height = %0.0f')));
%title('coverage for 65 degree elevation angle');
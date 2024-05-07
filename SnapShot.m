%------------------------------------------------------------------------
%-----Fraunhofer diffraction from a rectangular aperture-----------------
%------------------------------------------------------------------------

clc
close all
clear all

%------------------------------------------------------------------------

lambda = 1.0e-2; %5.0e-7; %1e-2; 
k = (2*pi)/lambda; % wavelength of light in vaccuum

a = 20.0; %1.0e-4; %10; 
b = 20.0; %1.0e-4; %10; % dimensions of diffracting rectangular aperture
                % a is along Y and b is along X
                
Io = 15; %300000; %8; %1.0e1; %1000.0; % relative intensity
z = 1.0e4 + 191.1; %1.0; %3.0e4; % distance of screen from aperture

Y = (-40:0.33:40);
%X = Y;
X = (-40:0.17:40); % coordinates of screen

beta=b*X/(lambda*z); 
gamma=a*Y/(lambda*z); % intermediate variable

 % diffracted intensity

I = zeros(length(Y),length(X));

temp1 = a*a*b*b/(lambda*z*lambda*z);

temp2 =  2*a*b/(lambda *z);

for i=1:length(Y)
    for j=1:length(X)
        % I(i,j) = Io * (sinc(beta(j))*sinc(gamma(i))).^2;
        I(i,j) = Io * (1 + temp1*(sinc(beta(j))*sinc(gamma(i))).^2  - temp2 * (sinc(beta(j)) * sinc(gamma(i)) * sin(pi*((X(j)*X(j) + Y(i)*Y(i) + 2*z*z)))/(lambda *z)));
    % I(i,j)=Io*(1 + ((sinc(beta(i)).^2).*(sinc(gamma(j))).^2)*a*a*b*b/(lambda*z*lambda*z) - sinc(beta(i)) * sinc(gamma(j)) * sin(pi*(X(i)*X(i) + Y(j)*Y(j) + 2*z*z)/(lambda *z)) * 2 *a * b/(lambda *z) );
    end
end

%------------------------------------------------------------------------
   
 %figure(1)
 %imshow(Y, X, I)
 %caxis('auto');
 
 %colormap summer
 caxis([min(min(I)) max(max(I))])
 
 %I1 = I/max(max(I))*64;
 image(Y, X, I)
 colorbar
 
 %title('Fraunhofer Diffraction','fontsize',14)
 %fh = figure(1);
 %set(fh, 'color', 'white'); 
 %plot(Y,I(501,:));
 

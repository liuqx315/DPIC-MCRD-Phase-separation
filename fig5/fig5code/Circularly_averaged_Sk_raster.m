function [SK] = Circularly_averaged_Sk_raster(Data,k)
% % %
% Copyright: Zhenpeng Ge, Quan-xing Liu, East China Normal University.
% Licence: MIT Licence
% 
% This function computes circularly averaged structure factor for 2D data in
% matrix form.
% Input:     
% Data: a 2D matrix must be a binary image consisting of 0 and 1.
% k:       a vector, wave number. Usually, k =0:1:floor(min(size(Data))/2)
% Output:  
% SK, a vector, circularly averaged structure factor.
% % % 

img = double(Data) ;
[M1,M2] = size(img) ;
% img(img>0) = 1 ;
% img(img==0) = -1 ;
S = abs(fftshift(fft2(img))).^2/(M1*M2) ;
[m,n] = meshgrid(-M2/2:M2/2-1,-M1/2:M1/2-1) ;
phi = (0.05:0.05:100)'/100*2*pi ;
SK = zeros(length(k),1) ;
for ii=1:length(k)
    SK(ii,1) = ii;
    SK(ii,2) = mean(interp2(m,n,S,k(ii)*cos(phi),k(ii)*sin(phi)));
    
end

SK=SK(3:end,:); % remove lower wavenumber

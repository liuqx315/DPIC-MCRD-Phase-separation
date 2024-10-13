% Cite this code as:
% Penny, G., Daniels, K.E., Thompson, S.E. (2013). Local pattern properties. 
% Local properties of patterned vegetation: quantifying endogenous and 
% exogenous effects. Proceedings of the Royal Society A: Mathematical, 
% Physical and Engineering Sciences. doi:10.xxxx/rspa.xxxx.xxxx
% ************************************************************************

function [imagecrop,L,D,params] = LocalPattProps(image,params)
% LocalPattProps.m
% This function calculates the local wavelength and direction of a 
% patterned image by recursively applying a Fourier algorithm to square
% sections of the image.
% The image is divided into a set of overlapping windows, and the
% wavelength and direction of the pattern in each window are calculated.
% The results from each window are combined into output arrays.
%
% INPUT:
% image: A binary patterned image (eg, vegetation / no vegetation)
% params: Parameters of the window and the min/max allowable frequencies
%   .w: width of the square window in px (~6-8 wavelengths)
%   .dL: length of increment for moving window in px
%   .minL: minimum allowable wavelength in px (typically set to 10-15 to ignore noise in image)
%   .maxL: maximum allowable wavelength in px (typ. w/4 < maxL <w/3, ie reject window if dominant
%       wavelength only shows three periods in the window)
%
% OUTPUT:
% L. (D.)
%   powmean: mean power of the pattern forming wavelengths from the power spectrum
%   powmax: max power for given radial band (angular band)
%   tile: resulting field of wavelengths in px (directions in rad)
%   uniqueness: uniqueness metric for the wavelengths (directions)
% L.within_range: 1 if the minL =< wavelength =< maxL, 0 otherwise


% PREPARE ADDITIONAL PARAMETERS
% M and N are the number of rows and columns in the window array
params.M = floor((size(image,1)-params.w)/params.dL)+1; %number of rows in window array 
params.N = floor((size(image,2)-params.w)/params.dL)+1; %number of columns in window array
params.O = params.M+params.w/params.dL-1; %M+w/dL-1; number of row steps in the image
params.P = params.N+params.w/params.dL-1; %N+w/dL-1; number of column steps in the image
params.n = params.w/params.dL; %w/dL; number of steps in a window

% Convert image to double precision
image = double(image);

% CALCULATE WAVELENGTH AND DIRECTION FOR EACH WINDOW
tstart = tic; %initialize close for the loop
for i = 1:params.M
    for j = 1:params.N
        
        % This function does the calculation
        [Lwin,Dwin] = LPPwindow(image,params,i,j); %calculate properties of the window

        % Store results of max, mean power and distance between furthest peaks
        L.powmax(i,j) = Lwin.powmax; %maximum radial power for window
        L.powmean(i,j) = Lwin.powmean; %mean radial power for window
        L.uniqueness(i,j) = Lwin.uniqueness; %uniqueness metric
        D.powmax(i,j) = Dwin.powmax; %maximum angular power for window
        D.powmean(i,j) = Dwin.powmean; %mean angular power for window
        D.uniqueness(i,j) = Dwin.uniqueness; %uniqueness metric
        L.within_range(i,j) = Lwin.within_range;
            
        % Store wavelength and direction of maximum power
        L.tile(i,j) = Lwin.tile; %wavelength
        D.tile(i,j) = Dwin.tile; %direction

    end
    fprintf([num2str(round(i/params.M*100*10)/10),'%% complete, ',num2str(toc(tstart)),' s\n'])
end

% Trim the image so that it aligns with the last row/column of windows
imagecrop = image(1:params.O*params.dL,1:params.P*params.dL);

end


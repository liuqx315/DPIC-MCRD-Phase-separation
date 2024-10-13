% Cite this code as:
% Penny, G., Daniels, K.E., Thompson, S.E. (2013). Local pattern properties. 
% Local properties of patterned vegetation: quantifying endogenous and 
% exogenous effects. Proceedings of the Royal Society A: Mathematical, 
% Physical and Engineering Sciences. doi:10.xxxx/rspa.xxxx.xxxx
% ************************************************************************
%
% This script is used to analyze the local pattern wavelength and
% orientation of a binarized pattern image. It does so by:
% Importing a binarized, patterned image and defining parameters
% for the analysis (#1). Calling LocalPattProps.m to calculate the 
% wavelength and orientation (direction) for each window specified by the 
% image and window size (#2). Plotting the mean power of pattern-forming 
% wavelengths for each window (#3). Allowing the user to select a minimum
% power threshold (MinPower) and minimum patch size (MinPatchSize to remove
% undesired windows from the results. The window results are then merged to
% create a final array of the local pattern wavelength, direction, and
% uniqueness (#5). The results are then plotted (#6).
%
% INSTRUCTIONS:
% Run the block of code from each step sequentially to perform the analysis
% and view the results. ExamplePattern1.tiff and ExamplePattern2.tiff can
% both be run with the same parameters.

%% 1. SET IMAGE AND CHOOSE WINDOW PARAMETERS
image = flipud(imread('ExamplePattern1.tiff')); %this should be a binary image of vegetation/novegetation. convert to double to use pcolor
%THESE VALUES SHOULD BE SET MANUALLY
params.w = 60; %set the width of the window (in pixels)
params.dL = 10; %set the step length (px)
params.minL = 5; %set the max wavelength (px)
params.maxL = 18; %set the min wavelength (typ w/3) (px)

%% 2. CALCULATE WAVELENGTH, DIRECTION, UNIQUENESS FOR EACH WINDOW
[imcrop,L,D,params2] = LocalPattProps(image,params);

%% 3. PLOT MEAN POWER (FROM POWER SPECTRUM) OF PATTERN-FORMING FREQUENCIES
% resize imagecrop if larger that 2000^2 pixels
numpx = size(imcrop,1) * size(imcrop,2);
if numpx > 2000^2; imcrop = imresize(imcrop,2000/sqrt(numpx),'nearest'); end

hf = figure; subplot(1,2,1)
pcolor(imcrop),shading flat, colorbar,title('Patterned image')
subplot(1,2,2)
pcolor(L.powmean),shading flat,colorbar,title('Corresponding mean power of pattern forming frequencies')

%% 4. DETERMINE MINIMUM POWER THRESHOLD AND MINIMUM PATCH SIZE THRESHOLD
MinPower = 2.5e10; %THIS VALUE SET MANUALLY
MinPatchSize = 20; %THIS VALUE SET MANUALLY

%% 5. MERGE THE OVERLAPPING WINDOWS
[L_merge,D_merge] = LPPmerge(L,D,params2,MinPower,MinPatchSize);

%% 6. PLOT THE RESULTS
figure,pcolor(L_merge.final),shading flat,colorbar,title('Local pattern wavelength')
figure,pcolor(D_merge.final),shading flat,colorbar,title('Local pattern direction')
figure,pcolor(imcrop),shading flat, colorbar,title('Patterned image')


% Cite this code as:
% Penny, G., Daniels, K.E., Thompson, S.E. (2013). Local pattern properties. 
% Local properties of patterned vegetation: quantifying endogenous and 
% exogenous effects. Proceedings of the Royal Society A: Mathematical, 
% Physical and Engineering Sciences. doi:10.xxxx/rspa.xxxx.xxxx
% ************************************************************************

function [Lwin,Dwin] = LPPwindow(image,params,i,j)
% LPPwindow.m
% This function takes a patterned image as input, crops it to a single 
% window (specified by i,j) and computes 2D FFT and power spectrum. It then 
% extracts the dominant wavelength and orientation from the power spectrum 
% as outlined in Penny et al (2013). The uniqueness metric is also computed 
% for both the wavelength and orientation. The mean power and max power of 
% the window are also returned.
%
% INPUT:
% image: binary image of the pattern.
% params: parameters for the calculations. See LocalPattProps.m for details
% i,j: The indices of the window to be analyzed
%
% OUTPUT:
% Lwin. (Dwin.)
%   powmean: mean power of the pattern forming wavelengths from the power spectrum
%   powmax: max power for given radial band (angular band)
%   tile: resulting field of wavelengths in px (directions in rad)
%   uniqueness: uniqueness metric for the wavelengths (directions)
% Lwin.within_range: 1 if the minL =< wavelength =< maxL, 0 otherwise

w = params.w;
dL = params.dL;
minL = params.minL;
maxL = params.maxL;

%% GET POWER SPECTRUM
% [~,fftpower] = FFTwindows(image,w,dL,i,j); 
im_hold = image(dL*(i-1)+1:dL*(i-1)+w,dL*(j-1)+1:dL*(j-1)+w); %retrieve the window
image_window = im_hold - mean2(im_hold); %subtract the mean from the window
imfft = fftshift(fft2(image_window)); %get FFT
fftpower = abs(imfft).^2; %calculate power

%% FIND RADIAL TOTAL POWER AND ANGULAR AVERAGE POWER
% [pow_avg_th,pow_tot_r,radii,theta] = FFTmaxpower(fftpower,minL,maxL);

%note: IMAGE MUST BE SQUARE
mindim = min(size(fftpower));
fmin_index = mindim/maxL; %minimum allowable frequency
fmax_index = mindim/minL; %maximum allowable frequency

%create grids of radii and angles corresponding to the power array
dimy = size(fftpower,1);
dimx = size(fftpower,2);
[X,Y] = meshgrid(-dimx/2:dimx/2-1,-dimy/2:dimy/2-1);
aa = atan(Y./X);
aa(X<0) = aa(X<0) + pi;
aa(aa<0) = aa(aa<0) + 2*pi; %grid of angles 0 to 2pi
rr = sqrt(X.^2 + Y.^2); %grid of radii

rtol = 1.5; %rtol is the tolerance for a specified radius (frequency) it is needed to smooth the results due to the binning of the power spectrum. the radial increment is 1.
dtol = pi/16; %dth is the tolerance on the angle. or 2*dtol is the width of the sliver. it is also the width of the angle increment
        

%*******************
%DETERMINE ANGULAR AVERAGE POWER
theta = 0+dtol:dtol:pi; % theta is center of sliver
for k = 1:length(theta) %iterate from zero to pi
   %find points that lie within sliver
   index = aa >= theta(k)-dtol & ...
       aa < theta(k)+dtol & rr >= fmin_index & rr <= fmax_index;
   pts_th{k} = fftpower(index);
   pow_avg_th(k) = mean(pts_th{k});
end

%***********************
%DETERMINE RADIAL TOTAL POWER
radii = ceil(fmin_index)-1:floor(fmax_index)+1; %NOTE: 1 extra point on each side
for r = 1:length(radii)
    index = rr < radii(r)+rtol & rr >= radii(r)-rtol;
    pts_r = fftpower(index);
    pow_tot_r(r) = sum(pts_r); %mean(pts_r);
    pow_a_r(r) = mean(pts_r);
end

%% CHECK FOR VALID WINDOW
valid_window = ... %check to ensure that each window has exactly 1 radial and angular power maximum
    length(find(max(pow_tot_r) == pow_tot_r)) == 1 & ...
    length(find(max(pow_avg_th) == pow_avg_th)) == 1;

if valid_window
      
%% CALCULATE WAVELENGTH, DIRECTION, UNIQUENESS, MEAN & MAX POWER
% [th_a,r_t] = FFTcalcs(fftpower,pow_avg_th,pow_tot_r,radii,theta,minL,maxL);

%MAX THETA AND R
th_a.th_indmax = find(max(pow_avg_th) == pow_avg_th);
th_a.th_max = theta(th_a.th_indmax); %theta max
r_t.r_indmax = find(max(pow_tot_r) == pow_tot_r);
r_t.r_max = radii(r_t.r_indmax); % r max

%FIND MAX R (weighted average)
poi_rat = 0.75;
[r_t.r_wavg,r_t.poi_ind,r_t.furthest_peak] = wavgpoi(pow_tot_r,radii,r_t.r_indmax,poi_rat);

%FIND WAVELENGTH OF MAX POWER
N = size(fftpower,1);
%r_t lengthscale
f_r_twa = 1/N * r_t.r_wavg;
r_t.L_wavg = 1/f_r_twa;

%FIND DIRECTION OF MAX POWER from MAX THETA (weighted average)
thL = length(theta);
shift = thL/2 - th_a.th_indmax;
theta_shift = theta; %determine how far to shift theta to center it on band of max power
if shift < 0; theta_shift(1:-shift) = theta_shift(1:-shift)+pi; %shift first part by + pi
else theta_shift(thL-shift+1:end)= theta_shift(thL-shift+1:end)-pi; end %shift last part by - pi

theta_shift = circshift(theta_shift,[0, shift]); %shift theta to be centered on angular band of max power
pow_avg_th_shift = circshift(pow_avg_th,[0, shift]); %shift power to correspond to shifted theta
[th_a.th_wavg,th_a.poi_ind,th_a.furthest_peak] = wavgpoi(pow_avg_th_shift,theta_shift,thL/2,poi_rat);
%shift indices back to where they were
th_a.poi_ind = th_a.poi_ind - shift; 
th_a.poi_ind(th_a.poi_ind < 1) = th_a.poi_ind(th_a.poi_ind < 1) + thL;
th_a.poi_ind(th_a.poi_ind > thL) = th_a.poi_ind(th_a.poi_ind > thL) - thL;
%shift weighted average if necessary
if th_a.th_wavg > pi; th_a.th_wavg = th_a.th_wavg - pi; end
if th_a.th_wavg < 0; th_a.th_wavg = th_a.th_wavg + pi; end
if th_a.th_wavg > pi
    N = N; end

%CALCULATE UNIQUENESS
r_t.uniqueness = 1 - abs(r_t.furthest_peak - r_t.r_max)/max(r_t.furthest_peak,r_t.r_max); %max possible distance is max(r_t.furthest_peak,r_t.r_wavg)
th_a.uniqueness = 1 - abs(th_a.furthest_peak - th_a.th_max)/pi; %max possible distance is pi

%check to ensure wavelength is not first or last
if r_t.r_indmax ~= 1 & r_t.r_indmax ~= length(radii) %LENGTHSCALE CANNOT BE FIRST OR LAST
    Lwin.within_range = 1; %max power falls within acceptable range
else
    Lwin.within_range = 0; %max power fall outside acceptable range
end
%% RETURN THE PARAMETERS
    Lwin.powmax = max(pow_tot_r); %maximum radial power for window
    Lwin.powmean = mean(pow_tot_r(2:end-1)); %mean radial power for window
    Lwin.uniqueness = r_t.uniqueness; %uniqueness metric
    Dwin.powmax = max(pow_avg_th); %maximum angular power for window
    Dwin.powmean = mean(pow_avg_th); %mean angular power for window
    Dwin.uniqueness = th_a.uniqueness; %uniqueness metric
    
    %Store wavelength and direction of maximum power
    Lwin.tile = r_t.L_wavg; %lengthscale of max power (from weighted average of peak in radial power)
    Dwin.tile = th_a.th_wavg; %direction of max power(from weighted average of peak in angular power)
    
else
    
    Lwin.powmax = NaN;
    Lwin.powmean = NaN;
    Lwin.uniqueness = NaN;
    Dwin.powmax = NaN;
    Dwin.powmean = NaN;
    Dwin.uniqueness = NaN;
    Lwin.within_range = NaN;
    
    %Store wavelength and direction of maximum power
    Lwin.tile = NaN;
    Dwin.tile = NaN;
        
end
end


%% WEIGHTED AVERAGE FROM POINTS OF INTEREST
function [wavg,poi_ind,furthest_peak] = wavgpoi(yvals,xvals,indmax,poi_rat)
%This functions finds a set of points in the peak of yvals based on the
%criteria specified below. It then finds the corresponding points in xvals
%and calculates the weighted average, weighting by the points in yvals.
%Note: xvals must be in increasing order and must correspond to yvals!
%
%yvals: vector of dependent values
%xvals: vector of independent values or indices
%indmax: index of max(yvals)
%poi_rat: ratio of points to include in Points of Interest (eg., 0.75)
%
%wavg: is the weighted average of x within the POI
%poi_ind: indices of the POI
%uniqueness_dist: largest distance between max(yvals) and another local
%                 maximum > 0.75*max(yvals). Difference in indices, not
%                 units of xvals.

%Identify range around indmax from local minimum to local minimum
[maxima,maxima_rat,minima] = findmaxmin(yvals); %find minima
minlo = indmax - min(indmax - minima(minima<indmax)); %find minimum just below max
if isempty(minlo); minlo = 1; end;
minhi = indmax + min(minima(minima>indmax) - indmax); %find min just above max
if isempty(minhi); minhi = length(yvals); end;

%Make the rand symmetric around indmax
if minhi-indmax > indmax - minlo
    minhi = indmax + (indmax - minlo);
else
    minlo = indmax - (minhi - indmax);
end
%POI is the intersection of the points in the range and the points above poi_ratio
poi_ind = intersect(minlo:minhi,find(yvals > poi_rat*max(yvals))); %find points to use for weighted avg

%calculate weighted average
poi_y = yvals(poi_ind);
poi_x = xvals(poi_ind);
wavg = sum(poi_y.*poi_x)/sum(poi_y);

%find furthest peak that is greater than 0.75 max power
max75 = maxima(maxima_rat > poi_rat);
if ~isempty(max75)
    furthest_peak_ind = max75(find(max(abs(max75 - indmax))==abs(max75 - indmax),1)); %find furthest peak (take lower value if more than 1)
    furthest_peak = xvals(furthest_peak_ind);
else furthest_peak = NaN; end
end


function [max_ind,max_ratio,min_ind,min_ratio] = findmaxmin(datavector)
% This function takes in a vector of data. It returns the indices (max_ind) of points
% that are local maxima, as well as the ratio of the magnitude of these
% points to the magnitude of the absolute maximum (max_ratio). The first and last points
% in the vector cannot be local extrema (eg, if a vector has length
% 10, the maxima can only have max_ind between 2 and 8). Indices of minima
% are returned in the same way.

dvapproach = datavector(2:end-1) - datavector(1:end-2);
dvdepart = datavector(2:end-1) - datavector(3:end);

dvincrease = sign(dvapproach);
dvdecrease = sign(dvdepart);

max_ind = find(dvincrease > 0 & dvdecrease > 0) + 1;
max_ratio = datavector(max_ind)/max(datavector(2:end-1));

min_ind = find(dvincrease < 0 & dvdecrease < 0) + 1;
min_ratio = datavector(min_ind)/max(datavector(2:end-1));

end
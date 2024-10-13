% Cite this code as:
% Penny, G., Daniels, K.E., Thompson, S.E. (2013). Local pattern properties. 
% Local properties of patterned vegetation: quantifying endogenous and 
% exogenous effects. Proceedings of the Royal Society A: Mathematical, 
% Physical and Engineering Sciences. doi:10.xxxx/rspa.xxxx.xxxx
% ************************************************************************

function [L_merge,D_merge] = LPPmerge(L,D,params,MinPower,MinPatchSize)
% LPPmerge.m
% This function merges the results from LocalPatternProps.m, which contain 
% overlapping windows. Additionally, this function includes filters to 
% remove windows that do not meet criteria specified by MinPower and 
% MinPatchSize. Once the filter is applied, the windows are merged such 
% that overlapping windows are averaged for each cell. This process is 
% applied to the pattern wavelength, orientation, and uniqueness metrics. 
% The edges of the results are then trimmed by half the width of a window 
% and the results are returned.
% 
% INPUTS:
% L: Structure output "L" from LocalPatternProps.m
% D: Structure output "D" from LocalPatternProps.m
% params: Structure output "params" from LocalPatternProps.m
% MinPower: The minimum power threshold. Any window with power below this
%      threshold will be discarded.
% MinPatchSize: The minimum patch size for a group of connected windows. 
%      In other words, if a group of 10 connected windows exist in 
%      isolation and MinPatchSize is 20, that group will be discarded.
% 
% OUTPUTS:
% L_merge (D_merge)
%    .final: wavelength field (direction field)
%    .uniqueness: uniqueness metric for wavelength (direction)

% WINDOW CHECK AND FILTER
WinCheckOrig = ~isnan(L.tile) & L.within_range==1; %establish baseline window check
WinCheck1 = L.tile < params.maxL & L.powmean > MinPower; %add min power criteria
WinCheck2 = bwmorph(WinCheck1,'close'); % remove holes
WinCheck3 = bwmorph(WinCheck2,'open'); % remove single points
CC = bwconncomp(WinCheck3,4); %calculate size of connected objects
WinCheck4 = WinCheck3;
for i = 1:CC.NumObjects %remove connected objects smaller than MinPatchSize
    if length(CC.PixelIdxList{i}) < MinPatchSize,
        WinCheck4(CC.PixelIdxList{i}) = 0;
    end
end
WinCheck_final = WinCheck1 & WinCheck4 & WinCheckOrig; %ensure that each window passes all tests

% WAVELENGTH MERGE
L_winCheck = L.tile;
L_winCheck(~(WinCheck_final)) = NaN;
%merge
fprintf('Wavelength window merge beginning\n')
Lwin_merge = nanmean(win_merge(L_winCheck,params.n),3); 

% DIRECTION MERGE
D_winCheck = D.tile;
D_winCheck(~(WinCheck_final)) = NaN;
%merge
fprintf('Wavelength window merge beginning\n')
Dwin_merge = DirectorAngleMatAvg(win_merge(D_winCheck,params.n),3);

% UNIQUENESS
% wavelength uniqueness
uniqueL_win = L.uniqueness;
uniqueL_win(~(WinCheck_final)) = NaN;
fprintf('Wavelength quality window merge beginning\n')
uniqueL_merge = nanmean(win_merge(uniqueL_win,params.n),3);
% direction uniqueness
uniqueD_win = D.uniqueness;
uniqueD_win(~(WinCheck_final)) = NaN;
fprintf('Direction quality window merge beginning\n')
uniqueD_merge = nanmean(win_merge(uniqueD_win,params.n),3);

% TRIM EDGES OF MERGED RESULTS BY HALF WINDOW WIDTH
% Get distance from boundaries
PattDistObjBoundaries = bwDistBoundaries(~isnan(Lwin_merge),'noholes')*params.dL;
% Trim edges 
L_merge.final = Lwin_merge;
L_merge.final(PattDistObjBoundaries < params.w/2) = NaN;
D_merge.final = Dwin_merge;
D_merge.final(PattDistObjBoundaries < params.w/2) = NaN;
L_merge.uniqueness = uniqueL_merge;
L_merge.uniqueness(PattDistObjBoundaries < params.w/2) = NaN;
D_merge.uniqueness = uniqueD_merge;
D_merge.uniqueness(PattDistObjBoundaries < params.w/2) = NaN;

end

%% DIRECTOR ARRAY ANGLE AVERAGE
function [danglematavg] = DirectorAngleMatAvg(danglemat,dim,weights)
% DirectorAngleMatAvg calculates the average of a set of matrices of
% director angles, relying on DirectorAngleAvg to average the director
% angles.

if nargin < 3
    weights = ones(size(danglemat)); end
if dim ~= 3
    error('Function only meant to average over 3rd dimension'); end

danglematavg = nan(size(danglemat,1),size(danglemat,2));
for i = 1:size(danglemat,1)
    for j = 1:size(danglemat,2)
        danglematavg(i,j) = DirectorAngleAvg(danglemat(i,j,:),weights(i,j,:));
    end
    fprintf(['Director Angle Averaging ',num2str(round(i/size(danglemat,1)*100*10)/10),'%% complete\n'])
end

end

%% DIRECTOR ANGLE AVERAGE
function [angleavg] = DirectorAngleAvg(d_anglevec,d_weights)

% DirectorAngleAvg takes a vector of director angles (either -pi/2 to pi/2 or
% 0 to pi). It then stretches the angles to a 2pi range, takes the vector
% average, and computes the new angle.

d_anglevec = d_anglevec(~isnan(d_anglevec));

if nargin < 2
    d_weights = ones(size(d_anglevec));
else
    d_weights = d_weights(~isnan(d_anglevec));
end

if isempty(d_anglevec); angleavg = NaN;
else
    

if max(d_anglevec) - min(d_anglevec) > pi
    error('Must be director angle'); end
if min(d_anglevec) < 0
    quads = -1;
else
    quads = 1; end

%stretch director to fill 0 to 2pi
a_stretch = 2*d_anglevec;
%take weighted average of x and y components
ax = sum(cos(a_stretch).*d_weights)/sum(d_weights);
ay = sum(sin(a_stretch).*d_weights)/sum(d_weights);

%determine quadrant of angle for each cell
q2 = ax < 0 & ay > 0;
q3 = ax < 0 & ay < 0;
q4 = ax > 0 & ay < 0;

%calculate angle
a_str_avg = atan(ay./ax); %good for quadrants 1 & 4

%shift angle into desired quadrants
if q2 %shift from -q4 to +q2
    a_str_avg = a_str_avg + pi;
elseif q3 & quads == 1 %shift from +q1 to +q3
    a_str_avg = a_str_avg + pi;
elseif q3 & quads == -1 %+q1 to -q3
    a_str_avg = a_str_avg - pi;
elseif q4 & quads == 1 %-q4 to q4
    a_str_avg = a_str_avg + 2*pi;
end

angleavg = a_str_avg/2;

end
end

%% WINDOW MERGE
function [merge] = win_merge(windowvar,n)
% This function combines the overlapping windows into a matrix.
% The matrix is 3D, and when averaged over the 3rd dimension, the resulting
% array is the result of averaging all the windows for each cell.

if nargin < 2
    error('Please input all variables'); end

[M,N] = size(windowvar);

O = M+n-1;
P = N+n-1;

merge = NaN(O,P,(n)^2);
for k = 1:n
    for l = 1:n
        ind = (k-1)*n+l;
        merge(k:end-n+k,l:end-n+l,ind) = windowvar;
    end
    fprintf(['Merge ', num2str(round(k/n*100*100)/100),'%% complete \n'])
end

end


%% DISTANCE TO BOUNDARIES
function [DistObjBoundaries] = bwDistBoundaries(obj,holeoption)
% This function takes in a 2D logical or numerical array of 1s and 0s,
% computes the boundaries around all objects (groups of 1s), then
% determines the distance to the nearest boundary at all points within objects.
%
% obj - object array, 1s represent pixels in objects
% holeoption - determines if hole boundaries are included in the distance
% calculation. 'holes' (default) includes hold boundaries. 'noholes'
% excludes hole boundaries.

if nargin < 2
    holeoption = 'holes'; end
if strcmp(holeoption,'noholes')
    obj = imfill(obj,'holes');
end

% obj = ~isnan(Lwin_merge);
[xx,yy] = meshgrid(1:size(obj,2),1:size(obj,1));

[B,L,N,A] = bwboundaries(obj,holeoption);
DistObjBoundaries = NaN(size(obj));

for i = 1:N
    %determine boundaries
    Bound = B{i}; %object boundary
    if strcmp(holeoption,'holes')
        for j = find(A(:,i))'
            Bound = [Bound;B{j}];
        end
    end
    
    Pxvec = xx(L==i);
    Pyvec = yy(L==i);
    
    numpts = length(Pyvec) * length(B{i});
    if numpts < 100e6 %use the entire matrix for multiplication
        Px = repmat(Pxvec,[1,size(Bound,1)]);
        Py = repmat(Pyvec,[1,size(Bound,1)]);
        Bx = transpose(repmat(Bound(:,2),[1,size(Px,1)])); %row ind
        By = transpose(repmat(Bound(:,1),[1,size(Py,1)])); %col ind

        dist_all = sqrt((Px-Bx).^2+(Py-By).^2);
        dist = min(dist_all,[],2);
    
    elseif numpts < 100e8 %use portions of the matrix
        int = 1:1000:length(Pxvec)-1;
        int = [int,length(Pxvec)+1];
        dist = [];
        for k = 1:length(int)-1
            Px = repmat(Pxvec(int(k):int(k+1)-1),[1,size(Bound,1)]);
            Py = repmat(Pyvec(int(k):int(k+1)-1),[1,size(Bound,1)]);
            Bx = transpose(repmat(Bound(:,2),[1,size(Px,1)])); %row ind
            By = transpose(repmat(Bound(:,1),[1,size(Py,1)])); %col ind

            dist_all = sqrt((Px-Bx).^2+(Py-By).^2);
            dist = [dist;min(dist_all,[],2)];
            fprintf(['Object ',num2str(i),' of ',num2str(N),' '...
                num2str(k/(length(int)-1)*100),'%% done\n'])
        end
    else
        error('length is too long'); end
    DistObjBoundaries(L==i) = dist;
    fprintf(['Object ',num2str(i),' of ',num2str(N),' done\n'])
end
end
%% 1. SET IMAGE AND CHOOSE WINDOW PARAMETERS
clc;
clear all;
cd('D:\recent_work_jupyter\DNA_work1110\Wavelet\PreparaData_To_Fig4')

%pathin='IMG1412';
%pathin='D:\recent_work_jupyter\DNA_work1110\Wavelet\Exp Data2\20230704&05 Droplet fusion\20+0.3 400R\20+0.3 400R FAM';
%load picture
folder = "D:\recent_work_jupyter\DNA_work1103\Exp Data2\20230704&05 Droplet fusion\10+0.3 400R\10+0.3 400R FAM";
files = dir(fullfile(folder, "*.png"));
filenames = {files.name};
Tmax=numel(filenames );
%im=dir(fullfile(pathin,'*.png'));
%[fsize z]=size(im);
%Dt=10; % min per frame
Dt=30; %timestep
dat=[];
for k = 10:Tmax;
    image = flipud(imread(fullfile(folder, filenames{k})));
    %this should be a binary image of vegetation/novegetation. convert to double to use pcolor
    %THESE VALUES SHOULD BE SET MANUALLY
    params.w =50; %set the width of the window (in pixels)
    params.dL = 10; %set the step length (px)
    params.minL = 1.0; %set the max wavelength (px)
    params.maxL = 50.0; %set the min wavelength (typ w/3) (px)
    %% 2. CALCULATE WAVELENGTH, DIRECTION, UNIQUENESS FOR EACH WINDOW
    [imcrop,L,D,params2] = LocalPattProps(image,params);

    %% 3. PLOT MEAN POWER (FROM POWER SPECTRUM) OF PATTERN-FORMING FREQUENCIES
    % resize imagecrop if larger that 2000^2 pixels
    numpx = size(imcrop,1) * size(imcrop,2);
    if numpx > 2000^2; imcrop = imresize(imcrop,2000/sqrt(numpx),'nearest'); end

%     hf = figure; subplot(1,2,1)
%     pcolor(imcrop),shading flat, colorbar,title('Patterned image')
%     subplot(1,2,2)
%     pcolor(L.powmean),shading flat,colorbar,title('Corresponding mean power of pattern forming frequencies')

    %% 4. DETERMINE MINIMUM POWER THRESHOLD AND MINIMUM PATCH SIZE THRESHOLD
    MinPower = 1.0e5; %THIS VALUE SET MANUALLY 2.5e10
    MinPatchSize = 10; %THIS VALUE SET MANUALLY

    %% 5. MERGE THE OVERLAPPING WINDOWS
    [L_merge,D_merge] = LPPmerge(L,D,params2,MinPower,MinPatchSize);

    %% 6. PLOT THE RESULTS
%     figure,pcolor(L_merge.final),shading flat,colorbar,title('Local pattern wavelength')
%     figure,pcolor(D_merge.final),shading flat,colorbar,title('Local pattern direction')
%     figure,pcolor(imcrop),shading flat, colorbar,title('Patterned image')
    L_merge.final=imresize(L_merge.final, 2, 'nearest');
    ML=L_merge.final(4:end-3,4:end-3);
    ML=ML*32/200; %400cm/3066 ML*0.05
%     nanmean(dd(:));
%     ML_crop = ML(5:end-5,5:end-5); % remove nan
%     以下代码在寻找每幅图的局部极大值点（波长）并计算均值
    i=1;
    for kw=1:min(size(ML))
        if ~isempty(max(findpeaks(ML(kw,:))))
            datT(i)=max(findpeaks(ML(kw,:)));
            i=i+1;
        end
        if ~isempty(max(findpeaks(ML(:,kw))))
            datT(i)=max(findpeaks(ML(:,kw)));
            i=i+1;
        end

    end
    time = k*Dt; 
    dat = [dat; k time mean(datT) std(datT)]; % dat=[dat; k*Dt nanmean(ML(:))]; 序列 时间 均值 标准差
end
%dlmwrite(strcat('Wave_Group6.dat'),dat, 'delimiter', '\t'); % pathin,
%dlmwrite('Wave_Group6.csv', dat, 'delimiter',',');
dlmwrite(strcat('DNA.dat'),dat, 'delimiter', '\t');
dlmwrite('DNA.csv', dat, 'delimiter',',');
% hf = figure; subplot(1,2,1)
% pcolor(imcrop),shading flat, colorbar,title('Patterned image')
% subplot(1,2,2)
% pcolor(L_merge.final),shading flat,colorbar,title('Corresponding mean power of pattern forming frequencies')

% TODO: cut a piece of image, compute radial averaged psd
clc; clear all;
fname='sklaw';
% fname='SDF';
dat1=load('Data_DNAnew2_30rho0.6/sklaw_data1.mat');
% dat2=load('Data/SDF_LSdata1.mat');

% load picture
% treename = dir('*.tif');
% treename = dir('vo*.jpg');
PSDTime=[];
for nfile=1:9
  
    dat=load(strcat('Data_DNAnew2_30rho0.6/',fname,'_data',num2str(nfile),'.mat'));
    Time=dat.Time;
    datDDA=dat.A;
    StartT=20; % 50 for DDA; 20 for SDF
    StartT0=StartT;
    for kk=StartT:1:size(datDDA,3); %size(treename,1); 
    % for images to analysis only
    % img = imread(treename(kk).name);
    % img_rgb = img;
    % img = rgb2gray(img);

    PT1=datDDA(:,:,kk);
    % mat2gray
    img = floor(255*mat2gray(PT1));

    img = 255 - uint8(img);
    mean_gray = mean(img(:));
    img = img - mean_gray;

    % perform 2D fft, compute psd
    [N,M] = size(img);
    % Compute power spectrum % with gpu
    img = single(img);
    % img = gpuArray(img);
    imgf = fftshift(fft2(img));
    % imgf = gather(imgf);
    % img = gather(img);
    imgfp = (abs(imgf)/(N*M)).^2;   

    %%% Adjust PSD size
    dimDiff = abs(N-M);
    dimMax = max(N,M);
    % Make square
    if N > M                                                                    % More rows than columns
        if ~mod(dimDiff,2)                                                      % Even difference
            imgfp = [NaN(N,dimDiff/2) imgfp NaN(N,dimDiff/2)];                  % Pad columns to match dimensions
        else                                                                    % Odd difference
            imgfp = [NaN(N,floor(dimDiff/2)) imgfp NaN(N,floor(dimDiff/2)+1)];
        end
    elseif N < M                                                                % More columns than rows
        if ~mod(dimDiff,2)                                                      % Even difference
            imgfp = [NaN(dimDiff/2,M); imgfp; NaN(dimDiff/2,M)];                % Pad rows to match dimensions
        else
            imgfp = [NaN(floor(dimDiff/2),M); imgfp; NaN(floor(dimDiff/2)+1,M)];% Pad rows to match dimensions
        end
    end
    halfDim = floor(dimMax/2) + 1;                                              % Only consider one half of spectrum (due to symmetry)

    %%% Compute radially average power spectrum
    [X Y] = meshgrid(-dimMax/2:dimMax/2-1, -dimMax/2:dimMax/2-1);               % Make Cartesian grid
    [theta rho] = cart2pol(X, Y);                                               % Convert to polar coordinate axes
    rho = round(rho)+1;


    rho = rho(:);
    imgfp= imgfp(:);

    Pf = zeros(max(rho),1);
    Counter = Pf;
%
    for index = 1:length(rho)
       Pf(rho(index)) = Pf(rho(index)) + imgfp(index); 
       Counter(rho(index)) = Counter(rho(index)) + 1; 
    end

    Pf = Pf./Counter;
    Pf(isnan(Pf)) = 0;

    % res = 0.273*2; %unit: m, change pixel to meter
    res = 1;
    Fs = 1/2*(1/res);
    L = length(Pf)-1;
    f = Fs*(0:L)'/L;

    % %note: only show wave length smaller than 1/10 of boxsize, 
    % %      to minimize finite size effect. Real picture do not have 
    % %      periodic boundary condition.
    % validindex = find(1./f<(dimMax*res/10));
    % f = f(validindex);
    % Pf = Pf(validindex);

    % setplot
    rdat=[f Pf];
    rdat(1,:)=[];
    rdat(rdat(:,1)>0.3,:)=[]; %note: only show wave length smaller than 1/3 of boxsize, 
    rdat(rdat(:,1)<3/size(img,1),:)=[];

    PSD=rdat(:,2);
    wavelength=1./rdat(:,1);
    
    [indx indy]=max(PSD);

    TimeWavelength(kk-StartT0+1)=wavelength(indy);

    % figure(1);
    % loglog(1./rdat(:,1),rdat(:,2),'bo-','LineWidth',2.0);
    % set(gca,'XColor','k','YColor','k','Xscale','log','Yscale','log');
    % set(gcf, 'Position', [0 0 700 600],'color','w','Visible','on');
    % set(gca, 'Position',[0.17 0.17 0.80 0.8],'LineWidth',2,'FontSize',20)
    % 
    % xlabel('wavelength',...
    %      'FontSize',24,'Interpreter','latex');
    % 
    % ylabel('power spectral density',...
    %     'FontSize',24,'Interpreter','latex','color','k');
    % saveas(gcf,strcat('pwd_',treename(kk).name(1:end-4),'.fig'));
    % saveas(gcf,strcat('pwd_',treename(kk).name(1:end-4),'.jpg'));
    PSDTime(:,kk-StartT0+1,nfile)=PSD; % xx, time, freqency
    
    end
    dataLS(1,:)=Time(StartT0+1:end);
    dataLS(nfile+1,:)=TimeWavelength;

end
dlmwrite(strcat(fname,'_wavelengthnew.csv'),dataLS,'delimiter','\t');
%% ==============

%% ==============
figure('Position', [10 10 600 500]);
timeseries=dataLS(1,:);
MeanWave=mean(dataLS(2:end,:));
SDWave=std(dataLS(2:end,:));
errorbar(dataLS(1,:),MeanWave,SDWave,'LineWidth',2);
% p1=shadedErrorBar(dataLS(1,:), dataLS(2:end,:), {@mean,@std}, 'lineprops', 'o-r');
% set(p1.edge,'LineWidth',2,'LineStyle','--')
% box on;
ylim([10 100])
yticks([10 20 40 100]);
xlabel('time, $t$','Interpreter','latex');
ylabel('wavelength, $\ell$','Interpreter','latex');
set(gca,'xscale','log','yscale','log','fontsize',18,'TickLength',[0.02 0.05],...
    'linewidth',2);

%% =====================================

dim=size(PSDTime);
COLOR=parula(dim(2)); COLOR(:,4)=0.7; FS=20;
temporal=mean(PSDTime,3);
figure('Position', [10 10 300 400]);
hold on
for kj=1:dim(2)
    loglog(wavelength,temporal(:,kj),'-','LineWidth',1.0,'color',COLOR(kj,:));
    set(gca,'XColor','k','YColor','k','Xscale','log','Yscale','log','TickLength',...
        [0.02 0.05],'FontSize',18,'linewidth',1.0,'layer','top');
%     set(gca, 'Position',[0.17 0.17 0.80 0.8]);
    xlabel('wavelength',...
         'FontSize',FS,'Interpreter','latex');
    ylabel('power spectral density',...
        'FontSize',FS,'Interpreter','latex','color','k');
end
box on;

xlim([4, 100]);
xticks([4 10 30 100])
% save2pdf(strcat('Fig3',fname,'_PSD'))
% dlmwrite(strcat(fname,'_LSPSD.csv'),temporal,'delimiter','\t');


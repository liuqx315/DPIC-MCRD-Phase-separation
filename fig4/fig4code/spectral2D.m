function [IrData,IomData]=spectral2D(Filename,Options);
% SPECTRAL2D - Two-dimensional spectral analysis of an image file
%    [Ir,Iom] = SPECTRAL2D('FILENAME',OPTIONS)
%    Spectral.m - Performes a spectral analysis as developed for ecology 
%    and geosciences by Renshaw and coworkers. FILENAME is the name of a
%    grayscale JPG or TIF image (TIF is better than JPG).
%
%    The image has to be square
%
%    The OPTIONS string contains character flags that do the following
%    'n' No graphs are displayed
%    'b' The bin distribution is shown
%    'c' Force a recalculation
%    't' Progress is given in text at the command line
%    'p' A progress bar shows progress
%
%    Ir  vectors of the radial spectrum and the two confidence limits
%    Iom  vectors of the radial spectrum and the two confidence limits
%
%    The function stores the periodogram in a file "filename.mat"
%    for quick future display of the results. You can force it to do a 
%    recalculation (for instance if you have different files with the 
%    same name) by giving 'p' as an option.
%    
%    For references about spectral analysis of 2D data see 
%    Renshaw and Ford 1984, Vegetatio 56:75-85
%    Muggulestone and Renshaw 1998, Computers & Geosciences 24
%    Couteron and Lejeune 2001, Journal of Ecology 89:616-628
%
%    Version 1. For an updated version, visit:
%    http://www.nioo.knaw.nl/homepages/koppel

global Ir Iom Chi_up_Ang Chi_dw_Ang Chi_up_Rad Chi_dw_Rad
fs=18;
on=1;off=0;
Graph=on;PlotBins=off;Calc=off;Pbar=off;ProgText=off;

if nargin>1,
    Graph=1-sum(Options=='n')>0;
    PlotBins=sum(Options=='b')>0;
    Calc=sum(Options=='c')>0;    
    Pbar=sum(Options=='p')>0;    
    ProgText=sum(Options=='t')>0;        
end;

% Get Screen dimensions and set Main Window Dimensions
x = get(0,'ScreenSize');
ScreenDim=x(3:4);
MainWindowDim=floor(ScreenDim.*[0.9 0.8]);

i=length(Filename);

Ext=lower(Filename([i-2:i])); % extracts the extension from the FileName

switch Ext
    case {'jpg','tif'},
       Filename=Filename([1:i-4]);
       Image=imread(Filename,Ext);            % Reads the image file if jpg or tif
       
    otherwise
       Image=imread(Filename,'tif');          % Reads the image file assuming tif
end;

Datafile=['constructed' Filename '.mat'];

AngleBinSize=15;            % The size of an angle bin in the angular spectrum
Resolution=0.45;             % 1 pixel is xxx umeters;

ShowSection=32;             % Size of the submatrix of I that is shown (a centre cutout)

ImInfo=imfinfo(Filename,Ext);

Colorpicture=1-strcmp(ImInfo.ColorType,'grayscale');

if Colorpicture,
   ImageBW=GrayImage(Image);
   [m,n,d]=size(Image); % Obtaining the size of the image array;
else
   ImageBW=Image;  % Obtaining the size of the image array;
   [m,n]=size(Image); 
   Image=zeros(m,n,3);
   Image(:,:,1)=ImageBW;
   Image(:,:,2)=ImageBW;
   Image(:,:,3)=ImageBW;
end;

% Checking if the image is square
if m~=n,
    disp('The image file needs to be square');
    beep;return;
end;

pmax=m/2;
qmax=m/2;

Y=zeros(m,n);
X=zeros(m,n);

a=zeros(pmax,qmax*2);
b=zeros(pmax,qmax*2);
I=zeros(pmax,qmax*2);
V=zeros(pmax,qmax*2);

Y(:)=ImageBW(:);                  % The image arraw is copied to a real-array that allows calculations

X = Y - mean(Y(:));             % Y is rescaled, so that the average equals zero;

% Below, the periodogram is calculated
if (exist(Datafile)) & (Calc==off),  % Checking if the periodogram has been calculated before
        load(Datafile);  % If so, the periodogram is loaded
else
    
    if ProgText, 
        disp('Processing image ...'); disp('Currently at   0%');
    end;
    if Pbar, h=waitbar(0,'Calculating periodogram');end;

    [t,s]=meshgrid(1:n,1:m);    
    
    for p=0:pmax,
       for q=0:qmax*2-1,
           a(p+1,q+1)=sum(sum(X.*cos(2*pi*(p*s/m+(q-qmax)*t/n))))/(m*n);
           b(p+1,q+1)=sum(sum(X.*sin(2*pi*(p*s/m+(q-qmax)*t/n))))/(m*n);
           I(p+1,q+1)=m*n*(a(p+1,q+1)^2+b(p+1,q+1)^2);

       end;
       if ProgText,
           disp(sprintf('\b\b\b\b\b%3.0f%%',p/pmax*100));
       end;
       if Pbar, waitbar(p/pmax,h);end;       
       drawnow;
    end;
    if Pbar, close(h);end;       
    save(Datafile,'I');                  
    
end;

% The variance is calculated

V=X.*X;                         % The variance array
Vtot=sum(V(:))/(m*n);           % The total variance

Dom=AngleBinSize;               % The AngleBinSize is copied into a variable with a shorter name
om_max=180/Dom;                 % The total number of bins

[pmax qmax]=size(I);            

rmax=min([pmax qmax]);

% Variable definitions
Ir  = zeros(1,rmax);            % The bins vor the radial spectrum is defined
Irc = zeros(1,rmax);            % The array with bin sizes for the radial spectrum is defined
Iom = zeros(1,om_max);          % The bins for the angular spectrum is defined
Iomc= zeros(1,om_max);          % The array with bin sizes for the angular spectrum is defined

Chi_up_Rad = zeros(1,rmax);     % The upper 95% confidence limits, Angular spectrum
Chi_dw_Rad = zeros(1,rmax);     % The lower 95% confidence limits, Angular spectrum
Chi_up_Ang = zeros(1,om_max);   % The upper 95% confidence limits, Angular spectrum
Chi_dw_Ang = zeros(1,om_max);   % The lower 95% confidence limits, Angular spectrum

Id=zeros(pmax,qmax);

% The bins are defined spatially
[q,p]=meshgrid(1:qmax,0:pmax-1);
qs=q-0.5*qmax-1;

% warning off MATLAB:divideByZero
om=ceil(atan(p./(qs))/pi*180/Dom)+om_max/2;
warning on;
om(1,(qmax/2+1):qmax)=0;

rs=round(sqrt(p.^2+(qs).^2));
rs(1,(qmax/2+1):qmax)=0;

% The angular bins are filled, and the number of cells counted
for i=1:om_max,
   Iom(i) = sum(sum((om==i).*I));
   Iomc(i)= sum(sum((om==i)));
end;    
    
% The radial bins are filled, and the number of cells counted
for i=1:rmax,
   Ir(i) = sum(sum((rs==i).*I));
   Irc(i)= sum(sum((rs==i)));
end;

% The bins are rescaled
Ir(:)=Ir(:)./Vtot./Irc(:);             
Iom(:)=Iom(:)./Vtot./Iomc(:); 

% Calculation of the Chi-squared critical values
Chi_up_Ang(:) = 1./(2*Iomc(:)).*chi2inv(0.975,Iomc(:)*2);
Chi_dw_Ang(:) = 1./(2*Iomc(:)).*chi2inv(0.025,Iomc(:)*2);
Chi_up_Rad(:) = 1./(2*Irc(:)).*chi2inv(0.975,Irc(:)*2);             
Chi_dw_Rad(:) = 1./(2*Irc(:)).*chi2inv(0.025,Irc(:)*2);            

% --------------------- Drawing the figures -------------------------------

% Doubling the original image for viewing, the lower part is rotated upwards
Itot = zeros(qmax,qmax);
Itot(pmax:2*pmax-1,:)=I(:,:);
Itot(1:pmax,:)=fliplr(flipud(I));
Itot(1:pmax,2:qmax)=Itot(1:pmax,1:qmax-1);

% Taking a subimage from the centre of the I image, with size ShowSection
s=ShowSection;
Ishow=zeros(ShowSection,ShowSection);
Ishow(:,:)=Itot((0.5*(m-s)+1):(0.5*(m+s)),(0.5*(n-s)+1):(0.5*(n+s)));

% Finally, the figure is drawn
if Graph==on,
%     Figure1=figure(... 
%            'Name','Spectral Analysis', ...
%            'NumberTitle','off', ...
%            'Position',[(ScreenDim-MainWindowDim)/2 MainWindowDim]...
%             ); 
    set(gcf,'position',[100 100 800 800],'color','w')
    tp=tiledlayout(2,2,'Tilespacing','Compact');
    tp.TileSpacing = 'compact';
    tp.Padding = 'compact';
    nexttile;
    
%     subplot(2,2,1);
    imagesc(flipud(uint8(Image))); 
    title(['Photograph -' Filename]);
    colormap('gray'); axis image; yticks([])%axis off;
    set(gca,'fontsize',fs,'Layer','top','linewidth',2,'TickDir','in','TickLength',[0.02 0.01]);          
           
%     subplot(2,2,2);
    nexttile;
    imagesc(max(Ishow(:))-Ishow); title('Periodogram');
    colormap('gray'); axis image; hold on;
    plot(s/2+1,s/2+1,'*');
    set(gca,'fontsize',fs,'Layer','top','linewidth',2,'TickDir','in','TickLength',[0.02 0.01]);

%     subplot(2,2,3); 
    nexttile;
    loglog(1:rmax,Ir(1:end),'-','markersize',10,'linewidth',2);hold on;
    kx=1:rmax; ky=Ir(1:end);
    save( ['rMax_' Filename  '.mat'], 'kx', 'ky');
%     yline([Chi_up_Rad(end) Chi_dw_Rad(end)],'--','linewidth',2,'fontsize',18);
%     yline(Chi_up_Rad(5:end),'b.','markersize',10);
%     yline(Chi_dw_Rad(5:end),'b.','markersize',10);
%     axis([0 ShowSection 0 max(Ir)*1.1]);
%     title('Radial spectrum');
    xlabel('Wavenumber');
    ylabel('Radial spectrum');
    set(gca,'fontsize',fs,'Layer','top','linewidth',2,'TickLength',[0.02 0.01]);
    
%     subplot(2,2,4); 
    nexttile;
    plot(Dom:Dom:180, Iom,'-s','markersize',10,'linewidth',2);hold on;
    datax=Dom:Dom:180;
    datay=Iom;
    save( ['Result_' Filename  '.mat'], 'datax', 'datay');
%     yline([Chi_up_Ang(end) Chi_dw_Ang(end)],'--',{'Max','Min'});
%     plot(Dom:Dom:180, Chi_up_Ang,'b.','markersize',10);
%     plot(Dom:Dom:180, Chi_dw_Ang,'b.','markersize',10);
    axis([-0.1 181 0 max([Iom 1.81]*1.1)]);
    ylabel('Angular spectrum');
    xlabel('Angle (degrees)');
    set(gca,'fontsize',fs,'Layer','top','linewidth',2,'TickLength',[0.02 0.01]);
    
end;
save2pdf(Filename,gcf,600)

% The plotting of the bins, (needed for checking)
if PlotBins==on,
    Figure1=figure( 'Name','The wavenumber and angle bins', ...
                    'NumberTitle','off', ...
                    'Position',[25 350 980 350] );
    subplot(1,2,1);
    imagesc(rs);axis image;
    subplot(1,2,2);
    imagesc(om);axis image;
    save2pdf(Filename,gcf,600)
end;

IomData=[Dom:Dom:180; Iom; Chi_up_Ang; Chi_dw_Ang];
IrData=[1:rmax; Ir; Chi_up_Rad; Chi_dw_Rad];

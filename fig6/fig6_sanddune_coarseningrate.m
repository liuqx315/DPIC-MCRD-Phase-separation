clc; clear all;

data=imread('fig1.png');
im2=imread('fig2.png');

data= rgb2gray(data);
im2= rgb2gray(im2);
im2= im2(:,1:end-5); % give up boundary column

dat2=data(:,6:end-2);
dat3(find(dat2==90))=1;

im2(find(im2==122))=1;

A=zeros(size(dat2));
A(find(dat2==90))=1;

dat3=A.*double(dat2);


Yind1=[1 23 55 107 171 237 325 411 507 611 703 811 925 1042];
% Yind1 = [10.7, 15.08615, 19.47231, 23.85846, 28.24462, 32.63077, 37.01693, 41.40308, 45.78924, 50.17539, 54.56155, 58.9477, 63.33386, 68.53];
% Yind1 = int64(Yind1);
Yind2=[1 91 164 250 354 460 544 646];
outdat1=zeros(size(dat3,2),14);

outdat1(:,1)=1:size(dat3,2);

for k=1:13
    datT=dat3(Yind1(k):Yind1(k+1),:);
    for kk=1:size(datT,2)
        outdat1(kk, k+1)=min(find(datT(:,kk)==90));
    end
%     imshow(datT)
%     imwrite(datT,['mycurve' num2str(k) '.png'])
end
save('outdat1v2.txt','outdat1','-ascii')

outdat2=zeros(size(im2,2),7);
outdat2(:,1)=1:size(im2,2);
for k=1:7
    datT_IM2=im2(Yind2(k):Yind2(k+1),:);
    for kk=1:size(datT_IM2,2)
        outdat2(kk, k+1)=min(find(datT_IM2(:,kk)==1));
    end
end
outdat2(:,1)=linspace(0,1*80,size(outdat2,1));
save('outdat2.txt','outdat2','-ascii')
% plot(outdat2(:,1),outdat2(:,7))
x=outdat2(:,1);
y=outdat2(:,6)-mean(outdat2(:,6));

%%
days_series=[1	16	36	81	146	174	207	224	267	345	388	462	566	754	813	910	971 	1109	1191	1330]
Periodlinear=[0.0941 0.0855 0.0753 0.3884  0.2690 0.2704 0.2736 ]
PeriodnonLinear=[0.2708 0.3920 0.3882 0.3731 0.3471 0.3398 0.3214 0.3178 0.3229 0.3203 0.3116 0.2985 0.2660 ]
Period=[0.2715 0.3921 0.3885 0.3736 0.3466 0.341 0.3219 0.3178 0.3229 0.3203 0.3126 0.2984 0.2661];
outdat1(:,1)=linspace(0,1*80,size(outdat1,1));
%%
for kk=4:4

x=outdat2(:,1);
%y=outdat2(:,5);
y=outdat2(:,2)-mean(outdat2(:,2));
yorigin=outdat2(:,2);
ffy = fft(y);
ffy(1) = [];
n = length(ffy);
power = abs(ffy(1:floor(n/2))).^2; % power of first half of transform data
maxfreq = 1/4;                   % maximum frequency
freq = (1:n/2)/(n/2)*maxfreq;    % equally spaced frequency grid
periodF=1./freq;

% plot(1./freq,power)
[ij TF2] = max(power);
plot(1./freq,power,periodF(TF2),power(TF2),'r*')
periodF(TF2);

end
%%
periodtotal=[Periodlinear,PeriodnonLinear]

%%
% cftool 
% plot dashed lines for power law
%x1=logspace(0.2,1.4,10);
x1=logspace(2.3,log10(1350),10)
%y1=13*x1.^0.2;
y1=5.5*x1.^0.19;
text(350,20,'$\ell\approx t^{0.19}$','Interpreter','latex', 'fontsize',20);
hold on
%plot(1:13,2*pi./Period,'.','markersize',18);
timeseries=linspace(10.7,68.53,13)
%plot(timeseries,2*pi./Period,'.','markersize',18);
plot(days_series(8:20),2*pi./PeriodnonLinear,'.','markersize',20);
wavelength_time=[days_series;2*pi./periodtotal]';
% csvwrite('wavelength_time.csv', wavelength_time);
plot(x1,y1,'--','linewidth',1.5);
box on
%xlim([0.9 22]);
xlim([210 1450]);
%ylim([15 25]);
ylim([14 25]);
%xticks([1 2 5 10 20]);
xticks([ 210 350 600 1400]);
%yticks([15 30 50 80]);
yticks([14 17 20 25]);
FS=18;
ylabel('Wavelength, $\ell$ [m]','Interpreter','latex');
xlabel('Time, $t$ [days]','Interpreter','latex');
set(gca,'fontsize',16,'xscal','log','yscal','log','fontsize',16,'linewidth',2);
set(gca,'fontsize',FS,'linewidth',2,'xminortick','off','yminortick','off',...
    'ticklength',[0.025 0.01]);
set(gca,'FontName','Times'); set(gcf,'Color',[1,1,1]);
%save2pdf('Coarsing_dune.pdf');

% plot(outdat(:,1),outdat(:,2:end))
% imshow(dat3)
% imwrite(dat3,'mycurve.png')
%%

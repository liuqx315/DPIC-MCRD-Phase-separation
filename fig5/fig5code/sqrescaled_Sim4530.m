
clc; clear all;
%fname='SDF'
num=5;
fname='sklaw';
%dat1=load('Data/DDA_LSdata1.mat');
dat1=load('/Users/liulab/Documents/DNA_work/Exp Data2/20230704&05 Droplet fusion/rho45_30_005/sklaw_data1.mat');
datDDA=dat1.B;
Time=dat1.Time;
%StartT=50;

StartT=40;
pick_num=round(linspace(160,194,5));
y = ones(5, 1);
pick_num=pick_num'.*y;

%%selected time corresponding qmax
Kmax = zeros(num,1) ;
%storage sq(r,t),qmax
Sq_rt = [];
SqQmax2 = [];
Q_Qmax = [];


%%
for i=1:num
    kk=pick_num(i);
    PT1=datDDA(:,:,kk);
    PT2=mat2gray(PT1);

    k =(1:1:floor(min(size(PT2))/3))';
    [SK] = Circularly_averaged_Sk_raster(PT2,k);
    Sq_rt(:,1)=SK(:,1);
    Sq_rt(:,i+1)=SK(:,2);
    qmax=sum(SK(:,1).*SK(:,2))./sum(SK(:,2));
    Kmax(i)=qmax;
   
end

for i = 1:num
    SqQmax2(:,i) = Sq_rt(:,i+1)*Kmax(i)^2 ;
    Q_Qmax(:,i) = Sq_rt(:,1)/Kmax(i) ;
end

T=round(Time(pick_num));

save('Sim4530v3.mat','T','Sq_rt','Kmax','SqQmax2','Q_Qmax') ;
%%
%%
data=load("Sim4530v3.mat")
%data2=load("Sim4530.mat")

markers = {'o','s','d','v','+','^','>','*','x','p','h','<'};
Times = data.T;
% 生成颜色序列
numMarkers = length(markers);
colors = lines(numMarkers); % 使用 lines 函数生成颜色

figure('Position', [10 10 600 500]);
hold on
FS=18;
Q_Qmax=data.Q_Qmax;
SqQmax2=data.SqQmax2;
for i = 1:length(Times)
    plot(Q_Qmax(:,i),SqQmax2(:,i),markers{i},'MarkerSize',8,'MarkerFaceColor',colors(i,:)) ; % 
    hold on
end
hold on
% Q_Qmax2=data2.Q_Qmax;
% SqQmax2_2=data2.SqQmax2;
% 
% for i = 1:length(Times)
%     plot(Q_Qmax2(:,i),SqQmax2_2(:,i),markers{i+5},'MarkerSize',8,'MarkerFaceColor',colors(i+5,:)) ; % 
%     hold on
% end

xxx = 1.0:1:30 ;
yyy = 3.3*10^3.5*xxx.^(-4.6) ;
h1 = plot(xxx,yyy,'--','linewidth',2,'color','b') ;


box on
ylim([0.1,1e5])
xlim([0.1,100])
yticks([10^(-1) 10^1 10^3 10^5]);

xlabel('$q/q_{\rm max}$','Interpreter','LaTex','Fontsize',FS) ;
ylabel('Structure function, $S(q)q_{\rm max}^2$','Interpreter','LaTex','Fontsize',FS) ;
% h = legend('$1514s$', '$2291s$', '$3311s$', '$5012s$', '$7244s$', '$1514s$', '$2291s$', '$3311s$', '$5012s$', '$7244s$') ;
h = legend('$1514s$', '$2291s$', '$3311s$', '$5012s$', '$7244s$');
set(h,'Interpreter','latex','Fontsize',FS,'Box','off','Location','northeast','Position',[0.084614944458008,0.115800001144409,0.322051747639974,0.279199998855592],'NumColumns',1);
%text(0.1,12,'Exp(10+0.3)','fontsize',FS,'Interpreter','LaTex') ;
text(4,90,'Slope $\sim -4.6$','fontsize',FS,'rotation',-53,'Interpreter','LaTex') ;
text(0.16,13,'Simulations:','fontsize',FS,'Interpreter','LaTex') ;
%title(h,'$\rho=15 \qquad \rho=22.5$','Interpreter','LaTex')
title(h,'$\quad \qquad v_0=30, \rho =22.5$','Interpreter','LaTex')
set(gca,'xscale','log','yscale','log','linewidth',1,'fontsize',FS,'TickLength',[0.02 0.025]);
set(gca,'fontsize',FS,'linewidth',2,'xminortick','off','yminortick','off',...
    'ticklength',[0.025 0.01]);

 save2pdf('rescaledSq_sim4530v3.pdf',gcf,600) ;
% x=Q_Qmax(:,2);
% x=log10(x);
% y=SqQmax2(:,2);
% y=log10(y);

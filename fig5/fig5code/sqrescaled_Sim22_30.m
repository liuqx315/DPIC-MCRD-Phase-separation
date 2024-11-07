
clc; clear all;

num=5;
fname='rho22.5_v30'; 
Nrun=1;
dat1=load(strcat('../fig5data/',fname,'/sklaw_data',num2str(Nrun),'.mat'));

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
Mean_SqQmax2=mean(SqQmax2,3);
Mean_Q_Qmax=mean(Q_Qmax,3);
T=round(Time(pick_num));

save(strcat(fname,'.mat'),'T','Sq_rt','Kmax','Mean_SqQmax2','Mean_Q_Qmax') ;
%%
%%
data=load(strcat(fname,'.mat'))

markers = {'o','s','d','v','+','^','>','*','x','p','h','<'};
Times = data.T;
% 生成颜色序列
numMarkers = length(markers);
colors = lines(numMarkers); % 使用 lines 函数生成颜色

figure('Position', [10 10 600 500]);
hold on
FS=18;
Q_Qmax=data.Mean_Q_Qmax;
SqQmax2=data.Mean_SqQmax2;
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
yyy = 3.3*10^3.5*xxx.^(-4.0) ;
h1 = plot(xxx,yyy,'--','linewidth',2,'color','b') ;


box on
ylim([0.1,1e5])
xlim([0.1,100])
yticks([10^(-1) 10^1 10^3 10^5]);

xlabel('$q/q_{\rm max}$','Interpreter','LaTex','Fontsize',FS) ;
ylabel('Structure function, $S(q)q_{\rm max}^2$','Interpreter','LaTex','Fontsize',FS) ;
h = legend('$1514$ s', '$2291$ s', '$3311$ s', '$5012$ s', '$7244$ s');
set(h,'Interpreter','latex','Fontsize',FS,'Box','off','Location','northeast','Position',[0.084,0.135,0.322,0.279],'NumColumns',1);
text(5,90,'Slope $\sim -4.0$','fontsize',FS,'rotation',-63,'Interpreter','LaTex') ;
% text(0.16,20,'Simulations:','fontsize',FS,'Interpreter','LaTex') ;
% title(h,'$\quad \qquad v_0=30, \rho =22.5$','Interpreter','LaTex')
set(gca,'xscale','log','yscale','log','linewidth',1,'fontsize',FS,'TickLength',[0.02 0.025],'xminortick','on','yminortick','on');

save2pdf(strcat(fname,'.pdf'),gcf,600) ;
% x=Q_Qmax(:,2);
% x=log10(x);
% y=SqQmax2(:,2);
% y=log10(y);

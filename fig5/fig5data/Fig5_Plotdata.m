clc; clear all;

num=5;
fname='rho22.5_v30'; 
for Nrun=1:10 % take average value for each time
    dat1=load(strcat('./',fname,'/sklaw_data',num2str(Nrun),'.mat'));
    Nrun
    datDDA=dat1.B;
    Time=dat1.Time;
    
    nplot=size(datDDA,3);
    pick_num=round(linspace(150,200,5));
    y = ones(5, 1);
    pick_num=pick_num'.*y;
    
    %%selected time corresponding qmax
    Kmax = zeros(num,1) ;    
    %storage sq(r,t),qmax
    Sq_rt = [];
    TSqQmax2 = [];
    TQ_Qmax = [];
    %%
    for i=1:num
        kk=pick_num(i);
        PT1=datDDA(:,:,kk);
        PT2=mat2gray(PT1);
        % PT2=PT1-mean((PT1(:)));
        k =(1:1:floor(min(size(PT2))/3))';
        [SK] = Circularly_averaged_Sk_raster(PT2,k);
        Sq_rt(:,1)=SK(:,1);
        Sq_rt(:,i+1)=SK(:,2);
        qmax=sum(SK(:,1).*SK(:,2))./sum(SK(:,2));
        Kmax(i)=qmax;
       
    end
    
    for i = 1:num
        TSqQmax2(:,i) = Sq_rt(:,i+1)*Kmax(i)^2 ;
        TQ_Qmax(:,i) = Sq_rt(:,1)/Kmax(i) ;
    end
    SqQmax2(:,:,Nrun)=TSqQmax2;
    Q_Qmax(:,:,Nrun)=TQ_Qmax;
end
Mean_SqQmax2=mean(SqQmax2,3);
Mean_Q_Qmax=mean(Q_Qmax,3);
T=round(Time(pick_num))';
FS=18;
loglog(Mean_Q_Qmax,Mean_SqQmax2)
xlabel('$q/q_{\rm max}$','Interpreter','LaTex','Fontsize',FS) ;
ylabel('Structure function, $S(q)q_{\rm max}^2$','Interpreter','LaTex','Fontsize',FS) ;
set(gca,'xscale','log','yscale','log','linewidth',1,'fontsize',FS,'TickLength',[0.02 0.025],'xminortick','on','yminortick','on');

save(strcat(fname,'.mat'),'T','Sq_rt','Kmax','Mean_SqQmax2','Mean_Q_Qmax') ;
%% 

data=load(strcat(fname,'.mat'));

markers = {'o','s','d','v','^','>','*','x','p','h','<'};
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
    plot(Q_Qmax(:,i),SqQmax2(:,i),strcat('-',markers{i}),'linewidth',1,'MarkerSize',8,'MarkerFaceColor',colors(i,:)) ; % 
    hold on
end
hold on

xxx = 0.8:1:20 ;
yyy = 5.0*10^3.5*xxx.^(-4.0) ;
h1 = plot(xxx,yyy,'--','linewidth',2,'color','b') ;


box on
ylim([0.1,1e5])
xlim([0.1,20])
yticks([10^(-1) 10^1 10^3 10^5]);

xlabel('$q/q_{\rm max}$','Interpreter','LaTex','Fontsize',FS) ;
ylabel('Structure function, $S(q)q_{\rm max}^2$','Interpreter','LaTex','Fontsize',FS) ;
 h = legend(num2str(T),'Interpreter','latex');
% h = legend('$1514$ s', '$2291$ s', '$3311$ s', '$5012$ s', '$7244$ s');
set(h,'Interpreter','latex','Fontsize',FS,'Box','off','Location','northeast','Position',[0.084,0.135,0.322,0.279],'NumColumns',1);
text(5,100,'Slope $\sim -4.0$','fontsize',FS,'rotation',-53,'Interpreter','LaTex') ;
% text(0.16,20,'Simulations:','fontsize',FS,'Interpreter','LaTex') ;
% title(h,{['$\quad \qquad v_0=25, \rho =15.0$']},'Interpreter','LaTex')
set(gca,'xscale','log','yscale','log','linewidth',1,'fontsize',FS,'TickLength',[0.02 0.025],'xminortick','on','yminortick','on');

save2pdf(strcat(fname,'.pdf'),gcf,600) ;


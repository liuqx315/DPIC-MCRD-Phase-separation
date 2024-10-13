clc; clear all;

% photo number
num=5;

%%group:20+0.6
folder = "/Users/liulab/Documents/DNA_work/Exp Data2/20230704&05 Droplet fusion/20+0.6 400R/20+0.6 400R FAM";
Time=120+30*linspace(0,239,240);

pick_num=round(linspace(40,240,5));
y = ones(5, 1);
pick_num=pick_num'.*y;

%%selected time corresponding qmax
Kmax = zeros(num,1) ;

%load data
files = dir(fullfile( folder,"*.png"));
filenames = {files.name};
Tmax = numel(filenames);
%storage sq(r,t),qmax
Sq_rt = [];
SqQmax2 = [];
Q_Qmax = [];


for i=1:num
    kk=pick_num(i);
    PT1 = imread(fullfile(folder, filenames{kk}));
    PT2 = mat2gray(PT1);
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

T=120+(pick_num-1)*30;

save('Exp2006v3.mat','T','Sq_rt','Kmax','SqQmax2','Q_Qmax') ;
%%
data=load("Exp2006v3.mat")
%data2=load("Exp2006v2.mat")
markers = {'o','s','d','v','+','^','>','*','x','p','h','<'};
% 生成颜色序列
numMarkers = length(markers);
colors = lines(numMarkers); % 
Times = data.T;
figure('Position', [10 10 600 500]);
hold on
FS=18;

% for i = 1:length(Times)
%     plot(Q_Qmax(:,i),SqQmax2(:,i),markers{i},'MarkerSize',8,'MarkerFaceColor',colors(i,:)) ; % 
%     hold on
% end
% hold on
Q_Qmax2=data.Q_Qmax;
SqQmax2_2=data.SqQmax2;

for i = 1:length(Times)
    plot(Q_Qmax2(:,i),SqQmax2_2(:,i),markers{i+5},'MarkerSize',8,'MarkerFaceColor',colors(i,:)) ; % 
    hold on
end

xxx = 1.0:1:30 ;
yyy = 7.0*10^4*xxx.^(-4.1) ;
h1 = plot(xxx,yyy,'--','linewidth',2,'color','b') ;


box on
ylim([0.1,1e5])
xlim([0.1,100])
yticks([10^(-1) 10^1 10^3 10^5]);

xlabel('$q/q_{\rm max}$','Interpreter','LaTex','Fontsize',FS) ;
ylabel('Structure function, $S(q)q_{\rm max}^2$','Interpreter','LaTex','Fontsize',FS) ;
%h = legend('$1290s$', '$2790s$', '$4290s$', '$5790s$', '$7290s$', '$1290s$', '$2790s$', '$4290s$','5790s','7290s') ;
h = legend('$1290s$', '$2790s$', '$4290s$', '$5790s$', '$7290s$') ;
set(h,'Interpreter','latex','Fontsize',FS,'Box','off','Location','southwest','NumColumns',1);
%text(0.1,12,'Exp(10+0.3)','fontsize',FS,'Interpreter','LaTex') ;
text(10,100,'Slope $\sim -4.1$','fontsize',FS,'rotation',-53,'Interpreter','LaTex') ;
text(0.13,25,'Experiments:','fontsize',FS,'Interpreter','LaTex') ;
% title(h,'\ 10+0.3 \quad \ 20+0.6')
title(h,'\ 20+0.6')
set(gca,'xscale','log','yscale','log','linewidth',1,'fontsize',FS,'TickLength',[0.02 0.025]);
set(gca,'fontsize',FS,'linewidth',2,'xminortick','off','yminortick','off',...
    'ticklength',[0.025 0.01]);

  save2pdf('rescaledSq2006v3.pdf',gcf,600) ;


clc; clear all;
% photo number
num=5;
%% group:20+0.6
folder = "../fig5data/experimentdata/10+0.3 400R/10+0.3 400R FAM";
Time=120+30*linspace(0,239,240);

pick_num=round(linspace(20,240,5));
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
    k =(1:1:floor(min(size(PT2))/4))';
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

% save('Exp1003v3.mat','T','Sq_rt','Kmax','SqQmax2','Q_Qmax') ;
%%
data=load("Exp1003v3.mat")
markers = {'o','s','d','v','^','>','*','x','p','h','<'};
% 生成颜色序列
numMarkers = length(markers);
colors = lines(numMarkers); % 
Times = data.T;
figure('Position', [10 10 600 500]);
hold on
FS=18;

Q_Qmax2=data.Q_Qmax;
SqQmax2_2=data.SqQmax2;

for i = 1:length(Times)
    plot(Q_Qmax2(1:end-100,i),SqQmax2_2(1:end-100,i),strcat('-',markers{i}),'linewidth',1,'MarkerSize',8,'MarkerFaceColor',colors(i,:)) ; % 
    hold on
end

xxx = 0.8:1:30 ;
yyy = 7.0*10^4*xxx.^(-4.0) ;
h1 = plot(xxx,yyy,'--','linewidth',2,'color','b') ;


box on
ylim([0.1,1e5])
xlim([0.1,20])
yticks([10^(-1) 10^1 10^3 10^5]);

xlabel('$q/q_{\rm max}$','Interpreter','LaTex','Fontsize',FS) ;
ylabel('Structure function, $S(q)q_{\rm max}^2$','Interpreter','LaTex','Fontsize',FS) ;
h = legend(num2str(T),'Interpreter','latex');

set(h,'Interpreter','latex','Fontsize',FS,'Box','off','Location','southwest','Position',[0.084,0.135,0.322,0.279],'NumColumns',1);
text(5,400,'Slope $\sim -4.0$','fontsize',FS,'rotation',-53,'Interpreter','LaTex') ;
% text(0.13,25,'Experiments:','fontsize',FS,'Interpreter','LaTex') ;
% title(h,'\ 10+0.3')
set(gca,'xscale','log','yscale','log','linewidth',1,'fontsize',FS,...
    'TickLength',[0.02 0.025],'xminortick','on','yminortick','on');

save2pdf('rescaledSq1003v3.pdf',gcf,600) ;


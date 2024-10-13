%dlmwrite(strcat('new10+0.3','_Sk1.csv'),datSk1,'delimiter','\t');
realsk = readtable('new10+0.3_Sk1.csv');
%y=110.*x.^(-0.28);% 0.25~0.26? the slope is between -0.26~-0.25 
% loglog(SK(:,1),SK(:,2),'-o')
%plot(datSk1(:,1),datSk1(:,2),'-o');
realskarr = table2array(realsk);

% 将数组转换为双精度格式
Realsk = double(realskarr);

%plot(x,y,'--'); 
x=data(:,1);
y=data(:,2);

data = readtable('Data_DNAnew2_30rho0.6_Sk1.csv');
%y=110.*x.^(-0.28);% 0.25~0.26? the slope is between -0.26~-0.25 
% loglog(SK(:,1),SK(:,2),'-o')
%plot(datSk1(:,1),datSk1(:,2),'-o');
data_array = table2array(data);
data = double(data_array);
% 将数组转换为双精度格式
% 初始化存储平均值和方差的数组
meanqmax = zeros(length(data), 1); % 存储平均值
stdqmax = zeros(length(data), 1); % 存储方差

for i=1:length(data)
    qmaxValues = data(i, 2:end);
    meanqmax(i) = mean(qmaxValues, 'omitnan'); % 计算平均值，忽略NaN值
    stdqmax(i) = std(qmaxValues, 'omitnan'); 
end
% 将平均值和方差添加到datSk2数组后面作为两列
data = [data, meanqmax, stdqmax];


%plot(x,y,'--'); 
x=data(:,1);
y=data(:,2);
% dlmwrite(strcat('DNA','_Sk2.csv'),datSk2,'delimiter','\t');
%%
figure('Position', [10 10 600 500]);
hold on
FS=18;
%plot(datSk1(:,1),datSk1(:,2),'-o');
scatter(Realsk(:,1),Realsk(:,2), '*','color','b');
hold on
plot(x,y,'-.','LineWidth',1,Color=[0.56,0,0.9]); 
hold on
errorbar(Realsk(:,1),Realsk(:,2),2*Realsk(:,3),'color','[0.79,0.96,0.79]')
hold on
plot(data(29:72,1),meanqmax(29:72)*0.4,'D','markersize',6,'color','magenta');
hold on
errorbar(data(29:72,1),meanqmax(29:72)*0.4,2*stdqmax(29:72),'color','[0.06,1.00,1.00]');
hold on
x=logspace(2,3.9,1000);
y=69.*x.^(-0.26);
hold on
% loglog(SK(:,1),SK(:,2),'-o')
%plot(datSk1(:,1),datSk1(:,2),'-o');
% scatter(Realsk(:,1),Realsk(:,2), '*','color','b');
% hold on
% plot(x,y,'-.','LineWidth',1,Color=[0.56,0,0.9]); 
% hold on
% errorbar(Realsk(:,1),Realsk(:,2),2*Realsk(:,3),'color','[0.79,0.96,0.79]')
text(200,14,'$q_{max}\propto t^{-0.26}$','fontsize',FS,'Interpreter','latex','rotation',-35);
text(150,7,'$q_{max}=\frac{\int qS(q)dq}{\int S(q)dq}$','fontsize',22,'Interpreter','latex','rotation',0);
box on;
box on;
%
%ylim([10,200])
ylim([5,25])
xlim([100,8500])
%xlim([500,8000])
yticks([6 10 15 25])
xticks([100 500  3000  8000])
xlabel('Time, $t\ (s)$', 'fontsize',22,'Interpreter','latex')
ylabel('$q_{\rm max}$','fontsize',22,'Interpreter','latex')
set(gca,'xscale','log','yscale','log','linewidth',1,'fontsize',FS,'TickLength',[0.02 0.025]);
set(gca,'fontsize',FS,'linewidth',2,'xminortick','off','yminortick','off',...
    'ticklength',[0.025 0.01]);
%save2pdf('SK_law')
%%save picture as pdf
% filename = 'SK_law_DNA_realdata_1204.pdf';
% saveas(gcf, filename);
clc; clear all;
clc; clear all;
%fname='SDF'
fname='sklaw';
%dat1=load('Data/DDA_LSdata1.mat');
dat1=load('Data_DNAnew2_30rho0.6/sklaw_data1.mat');
datDDA=dat1.A;
Time=dat1.Time;
%StartT=50;
StartT=40;
datSk1=[];
datSk2=[];
ij=1;
%%
for nfile=1:9
    dat=load(strcat('Data_DNAnew2_30rho0.6/',fname,'_data',num2str(nfile),'.mat'));
    %dat=load(strcat('Data/',fname,'_LSdata',num2str(nfile),'.mat'));
    Time=dat.Time;
    datDDA=dat.A;
    ij=1;
    StartT=40; % 50 for DDA; 20 for SDF
    StartT0=StartT;
    for kk=StartT:2:size(datDDA,3);
        PT1=datDDA(:,:,kk);
        PT2=mat2gray(PT1);

        k =(1:1:floor(min(size(PT2))/3))';
        [SK] = Circularly_averaged_Sk_raster(PT2,k);
        qmax=sum(SK(:,1).*SK(:,2))./sum(SK(:,2));
        datSk1(ij,1)=Time(kk);
        datSk1(ij,2)=qmax;
        datSk2(ij,1)=Time(kk);
        datSk2(ij,nfile+1)=qmax;

        ij=ij+1;
    end
end


%%

dlmwrite(strcat('Data_DNAnew2_30rho0.6','_Sk1.csv'),datSk2,'delimiter','\t');
%dlmwrite(strcat('LALI','_Sk2.csv'),datSk2,'delimiter','\t');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% ===
% figure('Position', [10 10 600 500]);
% hold on
% FS=18;
%x=logspace(1,5,50);
%y=400.*x.^(-0.25);
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

errorbar(data(:,1),meanqmax,stdqmax,'LineWidth',2);
%plot(x,y,'--'); 
x=data(:,1);
y=data(:,2);
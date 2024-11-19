clc; clear all;
clc; clear all;
%fname='SDF'
fname='sklaw';
%dat1=load('Data/DDA_LSdata1.mat');
dat1=load('rho20_32_3/sklaw_data1.mat');
datDDA=dat1.u;
Time=dat1.Time;
%StartT=50;

StartT=10;
datSk1=[];
datSk2=[];
ij=1;

%%
for nfile=1:5
    dat=load(strcat('rho20_32_3/',fname,'_data',num2str(nfile),'.mat'));
    %dat=load(strcat('Data/',fname,'_LSdata',num2str(nfile),'.mat'));
    Time=dat.Time;
    datDDA=dat.u;
    ij=1;
    StartT=40; % 50 for DDA; 20 for SDF
    StartT0=StartT;
    for kk=StartT:1:size(datDDA,3);
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


data = datSk2(:, 2:end); % 提取除时间列外的数据
meanqmax = mean(data, 2, 'omitnan'); % 计算每行数据的平均值
stdqmax = std(data, 0, 2, 'omitnan'); % 计算每行数据的标准差

data = [datSk2(:, 1), data, meanqmax, stdqmax]; % 添加平均值和标准差列
header = {'timeseries', 'data1', 'data2', 'data3', 'data4', 'data5','meanqmax', 'stdqmax'};
T = array2table(data, 'VariableNames', header);
writetable(T, 'rho_32_005.xls', 'Sheet', 1);
load handel.mat;
sound(y, Fs);

clear all
data2=readmatrix("rho_32_005.xls");
x=data2(:,1);
y=data2(:,end-1);
xlog=log10(x);
ylog=log10(y);


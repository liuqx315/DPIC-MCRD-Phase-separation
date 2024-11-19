%%%%if the scale from [15,40],note that we rescale it for clarity
clear all;
close all;
FS=18;
%dat=load('Wave_Group6.dat');
dat=load('DNA.dat');
% 读取CSV文件
data = readtable('sklaw_wavelengthnew.csv');

% 将数据保存为DAT文件
save('data_WL_0312.dat', 'data');
% 将表格数据转换为数组形式
data_array = table2array(data);

% 将数组转换为双精度格式
data_double = double(data_array);
timeseries=data_double(1,:);
MeanWave=mean(data_double(2:end,:));
SDWave=std(data_double(2:end,:));

%%
figure(1);
figure('Position', [10 10 600 500]);
set(gcf, 'position', [100 100 600 500],'color','w');

%loglog(dat(2:end,2)/60,5*dat(2:end,3),'*','markersize',8,'color','b')
loglog(dat(2:end,2),6.5*dat(2:end,3),'*','markersize',8,'color','b')
hold on
%errorbar(dat(2:end,2),dat(2:end,3),dat(2:end,4))
%loglog(dat(2:end,2),dat(2:end,3),'o','markersize',11)
%errorbar(dat(2:end,2),dat(2:end,3),dat(2:end,4))
errorbar(dat(2:end,2),6.5*dat(2:end,3),5*dat(2:end,4),'color','g')
hold on

%loglog((data_double(138:end,1)/1000),data_double(138:end,2),'D','markersize',6,'color','magenta')
loglog(timeseries(107:166),MeanWave(107:166),'D','markersize',6,'color','magenta')
hold on

hold on
%errorbar(data_double(138:end,1)/1000,data_double(138:end,2),data_double(138:end,3),'color','y')
errorbar(timeseries(107:166),MeanWave(107:166),SDWave(107:166),'color','y');
hold on
%% note that in order to combine the simulation results with lab data,we rescale the lab data by scaling factor 6 

x=linspace(330,4300,100);
x1=linspace(500,7000,100);
%y=13*x.^0.35;
 y=4.8*x.^0.28;
 y1=3.7*x.^0.28;
loglog(x,y,'r--','linewidth',2)
hold on
loglog(x1,y1,'b--','linewidth',2)
xlabel('time (hour)')
ylabel('spatial scale (mm)')

text(500,40,'$l\sim t^{0.28}$','Interpreter','latex','fontsize',18)

xlabel('Time, $t (s)$ ','Interpreter','latex');
ylabel('Spatial scale, $\ell$ [mm]','Interpreter','latex');
%ylim([10 35]);
 ylim([15 52]);
%xlim([0.9 400]);
xlim([300 10500]);
xticks([300 1000 3000 10000]);
yticks([ 15  20  30 50]);
%yticks([3 4 5 6 8]);


set(gca,'XScale','log','YScale','log');
set(gca,'fontsize',FS,'linewidth',2,'xminortick','off','yminortick','off',...
    'ticklength',[0.025 0.01]);
set(gca,'FontName','Times'); set(gcf,'Color',[1,1,1]);

% filename = 'DNA_local_wave_realdata_1118.pdf';
% 
% saveas(gcf, filename);
% filename2 = 'DNA_local_wave_realdata_1118.jpg';
% saveas(gcf, filename2);
 %% compute the growth rate through experiment data
clear all
% 定义文件夹路径
folderPath = '/Users/liulab/Documents/DNA_work/dispersion/1003/';
%data2
folderPath2 = '/Users/liulab/Documents/DNA_work/dispersion/dataset/';

%C = uisetcolor

% 定义文件名规则
files = dir(fullfile(folderPath, "*.mat"));
% data2
files2=dir(fullfile(folderPath2, "*.mat"));
filenames = {files.name};
Tmax = numel(filenames);
radius=[];
wavenumber=[];
h_k=[];
% load all the datasets
for i = 1:Tmax

    data = load(filenames{i});
    
    radius(:,i)=data.ky;
    
end
wavenumber=data.kx;
wavenumber=wavenumber';
h_k=[wavenumber,radius];
timeseries=120+30.*linspace(1,Tmax,Tmax);

h=[];
growth_rate=[];
am=[];
%%
%% group10+03
for j=1:length(wavenumber)
    h=log(radius(j,7:Tmax));
    linear_fit = fitlm(timeseries(7:Tmax), h);
   
    growth_rate(j)=linear_fit.Coefficients.Estimate(2);
    am(j)=linear_fit.Coefficients.Estimate(1);
end

%% group 20+0.6
%files2=dir(fullfile(folderPath2, "*.mat"));
addpath('/Users/liulab/Documents/DNA_work/dispersion/dataset/')  
filenames2 = {files2.name};
Tmax2 = numel(filenames2);
radius2=[];
wavenumber2=[];
h_k2=[];
% load all the datasets
for i = 1:Tmax2

    data2 = load(filenames2{i});
    
    radius2(:,i)=data.ky;
    
end
wavenumber2=data.kx;
wavenumber2=wavenumber2';
h_k2=[wavenumber2,radius2];
timeseries2=120+30.*linspace(1,Tmax2,Tmax2);

h2=[];
growth_rate2=[];
am2=[];
for j=1:length(wavenumber2)
    h2=log(radius2(j,1:7));
    linear_fit = fitlm(timeseries2(1:7), h2);
   
    growth_rate2(j)=linear_fit.Coefficients.Estimate(2);
    am2(j)=linear_fit.Coefficients.Estimate(1);
end
%% compute the theorectical lines of growth rate
% group1
beta1=3.75%3.75;
beta2=100%497%987;500
b1=25.0%25;
k2=2.0;
rho=4.0;%0.4
delta=0.003;

%%%
u=16.6%13%17.0;%17
v=3.9%7.5%11.0;%8.0
%%%
fu=beta1.*v-beta2.*rho./(k2+u+v)+rho.*beta2.*u./(k2+u+v)^2-1;
fv=b1+beta1*u+beta2*rho*u./(k2+u+v)^2;
%%%conserved form
k2=linspace(0,150,150);
tra=fu-fv-k2.^2.*(delta+1);
det_a=k2.^2.*(delta.*fv-fu+delta.*k2.^2);
discriminant = tra.^2-4.*det_a;

% 初始化 lambda1 和 lambda2
lambda1 = zeros(size(k2));
lambda2 = zeros(size(k2));

% 计算 lambda1 和 lambda2
for i = 1:length(k2)
    % if discriminant(i) >= 0
    %     lambda1(i) = real(0.5*(tra(i) + sqrt(discriminant(i))));
    %     lambda2(i) = real(0.5*(tra(i) - sqrt(discriminant(i))));
    % else
        lambda1(i) =real(0.5*(tra(i) + sqrt(discriminant(i)))); % 使用 NaN 来忽略复数根
        lambda2(i) = real(0.5*(tra(i) + sqrt(discriminant(i))));
    % end
end
%xlim([0.9 22]);
%
%%%group2
beta1=3.75%3.75;
beta2=100%497%987;
b1=25.0%25;
k2=2.0;
rho=16.0;%0.6
delta=0.0005;
u2=35.1%42.3%23.3 %17.0;%17
v2=2.7%2.4%11.0;%8.0
%%%
fu2=beta1.*v2-beta2.*rho./(k2+u2+v2)+beta2.*u2.*rho/(k2+u2+v2).^2-1;
fv2=b1+beta1*u2+beta2*rho*u2./(k2+u2+v2).^2;
%%%conserved form
k2=linspace(0,150,150);
y_zero = zeros(size(k2));

%k=linspace(0,100,100);
tra=fu2-fv2-k2.^2.*(delta+1);
det_a=k2.^2.*(delta.*fv2-fu2+delta.*k2.^2);
discriminant = tra.^2-4.*det_a;

% 初始化 lambda1 和 lambda2
lambda12 = zeros(size(k2));
lambda22 = zeros(size(k2));

% 计算 lambda1 和 lambda2
for i = 1:length(k2)
    % if discriminant(i) >= 0
        lambda12(i) = real(0.5*(tra(i) + sqrt(discriminant(i))));
        lambda22(i) = real(0.5*(tra(i) - sqrt(discriminant(i))));
    % else
    %     lambda12(i) = NaN; % 使用 NaN 来忽略复数根
    %     lambda22(i) = NaN;
    % end
end



%% normalized the growth rate by a scaling factor

% 绘制图像
max_lambda1=max(lambda1);
max_growth_rate=max(growth_rate);
scaling_factor=max_lambda1./max_growth_rate;
%data2
max_growth_rate=max(growth_rate2);
scaling_factor2=max_lambda1./max_growth_rate;
rescaled_gr1=growth_rate.*scaling_factor;
rescaled_gr2=growth_rate2.*scaling_factor2;
%%% calculate the mean
avg_wavenumber = [];
avg_growth_rate = [];
sd_growth_rate = [];
avg_growth_rate2 = [];
sd_growth_rate2 = [];
se_growth_rate=[];
se_growth_rate2=[];
% 每三个元素计算一次平均和方差
for i = 1:4:length(wavenumber)
    if i+3 <= length(wavenumber)
        % 计算当前组的平均 wavenumber
        current_avg_wavenumber = mean(wavenumber(i:i+3));
        avg_wavenumber = [avg_wavenumber, current_avg_wavenumber];

        % 计算当前组的平均 growth_rate
        current_avg_growth_rate = mean(rescaled_gr1(i:i+3));
        avg_growth_rate = [avg_growth_rate, current_avg_growth_rate];
        current_avg_growth_rate2 = mean(rescaled_gr2(i:i+3));
        avg_growth_rate2 = [avg_growth_rate2, current_avg_growth_rate2];

       % 计算当前组的 growth_rate 标准差
        current_sd_growth_rate = std(rescaled_gr1(i:i+3));
        sd_growth_rate = [sd_growth_rate, current_sd_growth_rate];
        current_sd_growth_rate2 = std(rescaled_gr2(i:i+3));
        sd_growth_rate2 = [sd_growth_rate2, current_sd_growth_rate2];

        % 计算当前组的 growth_rate 标准误差
        current_se_growth_rate = current_sd_growth_rate / sqrt(4);  % 因为每组有3个数据点
        se_growth_rate = [se_growth_rate, current_se_growth_rate];
        current_se_growth_rate2 = current_sd_growth_rate2 / sqrt(4);
        se_growth_rate2 = [se_growth_rate2, current_se_growth_rate2];
    end
end


%%
figure('Position', [10 10 900 320]);
%tiledlayout(1,1,TileSpacing="compact")


% 绘制 y=0 的直线'Color', [0.6, 0.6, 0.6], 'LineStyle', '--'
plot(k2./184, y_zero, 'Color', [0.6, 0.6, 0.6], 'LineStyle', '--','LineWidth', 1.5); % 'k--' 表示黑色虚线
hold on
FS=18;
% 绘制图像
% max_lambda1=max(lambda1);
% max_growth_rate=max(growth_rate);
% scaling_factor=max_lambda1./max_growth_rate;
% %data2
% max_growth_rate=max(growth_rate2);
% scaling_factor2=max_lambda1./max_growth_rate;
% rescaled_gr1=growth_rate.*scaling_factor;
% rescaled_gr2=growth_rate2.*scaling_factor2;
plot(k2./184, lambda1, 'color','[0.11 0.8 0.8]', 'LineWidth', 1.5);
%%%%%%%%%%%%%%%
hold on
plot(k2./184, lambda12, 'color','[1.0, 0.41, 0.16]', 'LineWidth', 1.5);


%ylim([-10 25]);
hold on
box on
hold on
% scatter(wavenumber./184, growth_rate .* scaling_factor, 'o', 'MarkerFaceColor', [0.68, 0.92, 1.00], 'MarkerEdgeColor', 'black');
% hold on
% scatter(wavenumber2./184, growth_rate2 .* scaling_factor2, 'D', 'MarkerFaceColor', [1.0, 0.41, 0.16], 'MarkerEdgeColor', 'black');
% 输出结果
% errorbar(avg_wavenumber./184,avg_growth_rate.*scaling_factor,var_growth_rate.*scaling_factor.*5,'o','MarkerFaceColor', [0.68, 0.92, 1.00], 'MarkerEdgeColor', 'black')
% errorbar(avg_wavenumber./184,avg_growth_rate2.*scaling_factor2,var_growth_rate2.*scaling_factor2.*5,'D', 'MarkerFaceColor', [1.0, 0.41, 0.16], 'MarkerEdgeColor', 'black')

errorbar(avg_wavenumber(:,1:30)./184,avg_growth_rate(:,1:30),sd_growth_rate(:,1:30),'o','MarkerFaceColor', [0.68, 0.92, 1.00], 'MarkerEdgeColor', 'black','LineWidth',1.5)
errorbar(avg_wavenumber./184,avg_growth_rate2,sd_growth_rate2,'D', 'MarkerFaceColor', [1.0, 0.41, 0.16], 'MarkerEdgeColor', 'black','LineWidth',1.5)


% 设置y轴和x轴的范围及刻度
ylim([-10 11]);
xlim([0 0.8]);
xticks(linspace(0, 0.8, 5));  % 设置5个刻度位置
xticklabels(sprintf('%.1f\n', get(gca, 'XTick')));
yticks([-10 0 10 20 30]);

% 设置图例和标签
legend('','Theoretical($\bar{u}=16.6,\bar{v}=3.9,\rho=4.0,\delta=0.005$)','Theoretical($\bar{u}=42.3,\bar{v}=2.4,\rho=20.0,\delta=0.0005$)', 'Experiment (10+0.3)', 'Experiment (20+0.6)','box', 'off','Interpreter', 'latex');
legend('','Theoretical($\bar{u}=16.6,\bar{v}=3.9$)','Theoretical($\bar{u}=35.1,\bar{v}=2.7$)', 'Experiment (10+0.3)', 'Experiment (20+0.6)','box', 'off','Interpreter', 'latex');
ylabel('Growth rate, $\sigma$ ','Interpreter','latex');
xlabel('Wave number, $q$ ','Interpreter','latex');
% % 设置图形属性
t=title( 'Wavelength, $\lambda$ ', 'Interpreter', 'latex', 'Color', 'magenta','fontsize', FS)


set(t, 'Position', [0.25, 20, 0]);
%set(t,'Position',get(t,'Position')+[0.0 1.5 0.5]);
set(gca, 'fontsize', FS, 'linewidth', 1.5, 'xminortick', 'off', 'yminortick', 'off', 'ticklength', [0.015 0.0025]);
set(gca, 'FontName', 'Times');
set(gcf, 'Color', [1, 1, 1]);
box off;
%%%set paper edge
set(gca, 'LooseInset', [0,0.1,0.03,0.1]);

ax1 = axes('Position', get(gca, 'Position'), 'XAxisLocation', 'top','YAxisLocation', 'right' ,'Color', 'none');

linkaxes([gca ax1], 'x');
ax1.XColor = 'red';  
ax1.YColor = 'k'% 设置上方x轴颜色为红色
%ax1.XLim = xlim./2;
ax1.XLim = xlim.*0.8;
ax1.XTick = linspace(0.1, 0.8, 5);  % 设置5个刻度位置
ax1.YTick = [];
%ax1.set_xticks([2,4,6,8,10])
ax1.XTickLabel = arrayfun(@(x) sprintf('%.1f', 2*pi./(x*1)), ax1.XTick, 'UniformOutput', false);  % 计算波长并设置标签
%xlabel(ax1, 'wavelength ($\lambda$)', 'Interpreter', 'latex', 'Color', 'red','Position','[450, 330]');
% 设置上边框为红色
ax1.FontName='Times';
ax1.FontSize=16;
ax1.LineWidth = 1.5;
ax1.XColor = 'magenta';
save2pdf('fig4v10.pdf');

 
% DRdata=[wavenumber,growth_rate',growth_rate2'];
% 
% header={'wavenum','1003growthrate','2006growthrate'};
% 
% file_name='DRdata.xls';
% writematrix(DRdata,file_name);
% datapr=readmatrix("DRdata.xls");
% 
% % 初始化新的数组
%%



% 
data1 = [avg_wavenumber.'/184, avg_growth_rate.', sd_growth_rate.', avg_growth_rate2.', sd_growth_rate2.'];
data2 = [k2.'/184, lambda1.', lambda12.'];
% 
writematrix(data1, 'fig4_experiment0704.xls');
 writematrix(data2, 'fig4_Theory0704.xls');
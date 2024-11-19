

clear all
data_f = xlsread('fn_data.xlsx',2);
data_g = readmatrix('gn_data.xls');

%% fitting f
x = data_f(:, 1);
y=data_f(:,2);
%% fitting g
y = data_g(:, end-2);
x=data_g(:,end);
%%
% linear model
lm_model = fitlm(x, y)
lm_bic = lm_model.ModelCriterion.BIC;
lm_aic = lm_model.ModelCriterion.AIC;
coefficients = lm_model.Coefficients.Estimate;

% 提取截距项和斜率项的估计值
a = coefficients(1);
b = coefficients(2);
scatter(x,y)
% 打印指数模型的数学形式和参数估计值
%fprintf('指数模型： y = %.4f * exp(%.4f * x)\n', a, b);
% quadratic模型
quad_model = fitlm(x, y, 'quadratic');
quad_bic = quad_model.ModelCriterion.BIC;
quad_aic = quad_model.ModelCriterion.AIC;

% logarithm
log_model = fitlm(log(x), y);
log_bic = log_model.ModelCriterion.BIC;
log_aic = log_model.ModelCriterion.AIC;
plot(x,y)

% exponetial
exp_model = fitlm(x, log(y));
exp_bic = exp_model.ModelCriterion.BIC;
exp_aic = exp_model.ModelCriterion.AIC;
% 获取指数模型的系数
coefficients = exp_model.Coefficients.Estimate;

% 提取截距项和斜率项的估计值
a = exp(coefficients(1));
b = coefficients(2);

% 打印指数模型的数学形式和参数估计值
%fprintf('指数模型： y = %.4f * exp(%.4f * x)\n', a, b);

% power law
power_model = fitlm(log(x), log(y));
power_bic = power_model.ModelCriterion.BIC;
power_aic = power_model.ModelCriterion.AIC;
coefficients1=power_model.Coefficients.Estimate;
c = exp(coefficients1(1));
d = coefficients1(2);

%  c*x/(x+a) 
model1 = fitnlm(x, y, @(b, x) b(1) .* x ./ (x + b(2)), [1, 1]);
model1_bic = model1.ModelCriterion.BIC;
model1_aic = model1.ModelCriterion.AIC;
coefficients1=model1.Coefficients.Estimate;
c = coefficients1(1);
d = coefficients1(2);
%%
%  c/x+a 
model2 = fitnlm(x, y, @(b, x) b(1) ./ (x + b(2)), [3, 3])
model2_bic = model2.ModelCriterion.BIC;
model2_aic = model2.ModelCriterion.AIC;
coefficients1=model2.Coefficients.Estimate;
c = coefficients1(1);
d = coefficients1(2);

%  c/x+c 
model3 = fitnlm(x, y, @(b, x) b(1) ./ (x + b(1)), [1])
model3_bic = model3.ModelCriterion.BIC;
model3_aic = model3.ModelCriterion.AIC;
coefficients1=model3.Coefficients.Estimate;
c2 = coefficients1(1);
%d = coefficients1(2);
%%

scatter(x,y);
hold on
x0=linspace(0,max(x)+1,100);
plot(x0,1.6./(1.6+x0))
ylim([min(y),max(y)]);
xlim([min(x)-0.3,max(x)+0.3]);
%xlim([0,1])
%xlim([50,100])
%%




% % 创建 AIC 和 BIC 的表格
% model_table = table({'Linear'; 'Quadratic'; 'Logarithmic'; 'Exponential'; 'Power-law'; 'c*x/(x+a)'; 'c/x+a'}, ...
%                     [lm_aic; quad_aic; log_aic; exp_aic; power_aic; modelx_aic; model2_aic], ...
%                     [lm_bic; quad_bic; log_bic; exp_bic; power_bic; model1_bic; model2_bic], ...
%                     'VariableNames', {'Model', 'AIC', 'BIC'});
% 
% % 打印表格
% disp(model_table)
%%%%%%%

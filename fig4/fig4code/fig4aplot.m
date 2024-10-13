
figure('Position', [1, 1, 900, 350])
t = tiledlayout(1,3,'TileSpacing','Compact','Padding','Compact');

h_k_plot=radius(:,7:Tmax);
%growh_rate_hk=growth_rate';
timeseries=120+30.*linspace(1,Tmax,Tmax);
Time=timeseries(7:Tmax);

index=[2,19,68]
wavelength_selected=[  578.0 ,60.9,17.0];

    % 初始化一个空表格
    data_table = table();

%%
for i =1:3
    %%yvalue

    ii=index(i)
    h_k_log=log(h_k_plot(ii,:));
    linear_fit = fitlm(Time, h_k_log);
    growth_rate=linear_fit.Coefficients.Estimate(2);
    am=linear_fit.Coefficients.Estimate(1);
    %%scatter plot
    %%%%%%%%%%%%%
     % 为当前index创建一个临时表格
    % temp_table = table(Time', h_k_plot(ii,:)', h_k_log', repmat(growth_rate, length(Time), 1), repmat(am, length(Time), 1), ...
    %     'VariableNames', {'Time', 'h_k_plot', 'h_k_log', 'GrowthRate', 'Am'});
    % 
    % % 将临时表格添加到主表格中
    % data_table = [data_table; temp_table];
    %%%%%%%%%%%%%%
    nexttile
    %scatter(Time,h_k_log,'o', 'MarkerFaceColor', [0.68, 0.92, 1.00], 'MarkerEdgeColor', 'black');
    scatter(Time,h_k_plot(ii,:),'o', 'MarkerFaceColor', [0.68, 0.92, 1.00], 'MarkerEdgeColor', 'black');
    hold on
    %plot(Time,growth_rate*Time+am,'color','[0.96,0.55,0.36]','LineStyle', '--','LineWidth', 1.5);
    plot(Time,exp(am).*exp(growth_rate.*Time),'color','[0.96,0.55,0.36]','LineStyle', '--','LineWidth', 1.5);
    %ylim([min(h_k_log)-0.5,max(h_k_log)+0.5])
    
    %legend('Experiment data','Theoretical line' , 'box', 'off');
    if i == 1
        ylim([min(h_k_plot(ii,:))-7,max(h_k_plot(ii,:))+7])
        ylabel('Amplitudes, $h_\lambda$', 'Interpreter', 'latex');
        yticks([60,80,110,140]);  % 设置5个刻度位置
        yticklabels(sprintf('%.0f\n', get(gca, 'YTick')));
    end
    if i == 2
        ylim([min(h_k_plot(ii,:))-7,max(h_k_plot(ii,:))+7])
        yticks([40,60,80,100]);  % 设置5个刻度位置
        yticklabels(sprintf('%.0f\n', get(gca, 'YTick')));
    end
    if i == 3
        ylim([15,30])
        yticks([15,20,25,30]);  % 设置5个刻度位置
        yticklabels(sprintf('%.0f\n', get(gca, 'YTick')));
    end
    xlabel('Time, $t\ (s)$ ','Interpreter','latex');
    % 根据y轴最小值和最大值，创建对数间隔的刻度
    % 例如，使用`logspace`设置对数间隔
    % 设置y轴为对数坐


    
    % yticks(linspace(min(h_k_plot(ii,:))-7,max(h_k_plot(ii,:))+7, 4));  % 设置5个刻度位置
    % yticklabels(sprintf('%.0f\n', get(gca, 'YTick')));
    %yticks(logspace(min(h_k_log)-0.5,max(h_k_log)+0.5, 4));
    % yticks(logspace(min(h_k_log)-0.5,max(h_k_log)+0.5, 4));
    % yticklabels(sprintf('%.1f\n', get(gca, 'YTick')));
    xlim([300,700])
    xticks([300,500,700]);
    xticklabels(sprintf('%.0f\n', get(gca, 'XTick')));

    % 设置边框和字体
    set(gca,'fontsize',FS,'yscale','log','linewidth',2,'xminortick','off','yminortick','off',...
    'ticklength',[0.02 0.0025]);

    % 设置背景色
    set(gcf, 'Color', [1, 1, 1]);

    % 设置标题
    title(sprintf('\\lambda=%.1f', wavelength_selected(i)));
    box on
end
save2pdf('growth_rate_sub3fig2v3.pdf');

% % 设置上边框为红色
% ax1.LineWidth = 1.5;
% ax1.XColor = 'red';
%writetable(data_table, 'data_output.xls');
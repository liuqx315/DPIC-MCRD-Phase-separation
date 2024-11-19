function [mean_qmax, var_qmax] = calculate_qmax_for_image(PT2)
    % 假设PT2是已经读取并转换为灰度的图像矩阵

    % 获取图像尺寸
    [height, width] = size(PT2);
    min_dimension = min(height, width);
    square_size = floor(min_dimension / 3);  % 正方形区域的大小

    % 定义五个正方形区域的起始点（左上角）
    starts = [1, 1;
              1, width - square_size + 1;
              height - square_size + 1, 1;
              height - square_size + 1, width - square_size + 1;
              floor((height - square_size) / 2) + 1, floor((width - square_size) / 2) + 1];
    
    qmax_values = zeros(5, 1);  % 存储每个区域的qmax值
    
    for idx = 1:5
        % 获取当前正方形区域
        row_start = starts(idx, 1);
        col_start = starts(idx, 2);
        square_region = PT2(row_start:row_start + square_size - 1, col_start:col_start + square_size - 1);
        
        % 计算当前正方形区域的qmax
        k =(1:1:floor(min(size(square_region))/3))';
        [SK] = Circularly_averaged_Sk_raster(square_region, k);
        qmax_values(idx) = sum(SK(:,1).*SK(:,2))./sum(SK(:,2));
    end
    
    % 计算平均值和方差
    mean_qmax = mean(qmax_values);
    var_qmax = var(qmax_values);
end

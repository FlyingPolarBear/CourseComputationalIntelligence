clear all
%% 数据预处理
data = csvread('iris.csv'); %读取图像
[I,J] = size(data); % I为数据量，J为数据维数
label = data(:,J);
data(:,J) = [];

[coeff,score,latent,tsquared] = pca(data);
data = score(:,1:2);
J = 2;

M = 20; %迭代次数
N = max(label); %聚类数
error = zeros(M,1);

% 计算每一维度的最大值和最小值,以及初始聚类中心
center = zeros(N,1); % 中心坐标
minmax = zeros(J,2); % 求最小值和最大值
for j = 1:J
    minmax(j,1) = min(data(:,j)); % 最小值
    minmax(j,2) = max(data(:,j)); % 最大值
    % 随机初始化聚类中心,每一行是一个初始化向量
    for n = 1:N
        center(n,j) = minmax(j,1) + (minmax(j,2) - minmax(j,1)) * rand();
    end
end

%% 迭代优化
for m = 1:M % 迭代M次
    %% 划分区域
    dis = zeros(I,N); % 距离矩阵
    sum_dis = zeros(I,N); % 计算距离过渡矩阵
    grade = zeros(I,N); % 隶属度矩阵
    for i = 1:I
        for n = 1:N
            % 计算距离
            dis(i,n) = exp(-norm(data(i,:)-center(n,:)));
        end
        Dis = sum(dis(i,:));
        for n = 1:N
            grade(i,n) = dis(i,n) / Dis;
        end
    end
    
    %% 更新聚类中心
    center_list = zeros(n,j);
    count = zeros(n,1);
    for n = 1:N
        center = sum(data*grade(:,n));
    end
    
    %% 计算误差
    err = 0;
    for i = 1:I
        if (near_center(i,1) - label(i,1) ~= 0)
            err = err + 1;
        end
    end
    error(m,1) = err / I;
end

scatter(data(:,1),data(:,2))
title('iris')

% %% 误差结果画图
% plot(error,'LineWidth',2)
% title('K-means错误率曲线');
% ylim([0,1]);
% xlabel('迭代次数');
% ylabel('错误率');
% grid on
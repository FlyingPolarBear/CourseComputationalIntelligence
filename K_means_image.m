 clear all
%% 数据预处理
image = imread('IMAGE1.png'); %读取图像
[I,J,K] = size(image); % I,J为图片大小，K为通道数
M = 20; %迭代次数
N = 3; %聚类数
% 计算每一维度的最大值和最小值,以及初始聚类中心
center = zeros(N,1); % 中心坐标
minmax = zeros(K,2); % 求最小值和最大值
for k = 1:K
    minmax(k,1) = min(min(image(:,:,1))); % 最小值
    minmax(k,2) = max(max(image(:,:,1))); % 最大值
    % 随机初始化聚类中心,每一行是一个初始化向量
    for n = 1:N
        center(n,k) = minmax(k,1) + (minmax(k,2) - minmax(k,1)) * rand();
    end
end
image = double(image); % 将图像矩阵转换为double类型

%% 迭代优化
for m = 1:M % 迭代M次
    m
%% 划分区域
dis = zeros(I*J,N); % 距离矩阵
near = zeros(I,J); % 近邻矩阵
near_center = zeros(I,J); % 聚类中心矩阵
near_center1 = zeros(I*J,1); % 编码后的聚类中心矩阵
for i = 1:I
    for j = 1:J
        % 每个点计算距离
        for n = 1:N
            dis((i-1)*J+j,n) = sqrt((image(i,j,1)-center(n,1))^2+...
            (image(i,j,2)-center(n,2))^2+(image(i,j,3)-center(n,3))^2); 
        end
        % 找到最近的聚类中心
        [near(i,j),near_center1((i-1)*J+j,1)] = min(dis((i-1)*J+j,:),[],2);
    end
end
for i = 1:I
    for j = 1:J
        near_center(i,j) = near_center1((i-1)*J+j,1); % 解码
    end
end

%% 更新聚类中心
center_list = zeros(n,k);
count = zeros(n,1);
for i = 1:I
    for j = 1:J
        for n = 1:N
            % 如果是近邻，则加入list
            if (near_center(i,j) == n)
                for k = 1:K
                    center_list(n,k) = center_list(n,k) + image(i,j,k);
                end
                count(n,1) = count(n,1) + 1;
            end
        end
    end
end
for n = 1:N
    % 如果存在亢余簇，则重新分配
    if (M - m > 2) && count(n,1) == 0
        for k = 1:K
            center(n,k) = minmax(k,1) + (minmax(k,2) - minmax(k,1)) * rand();
        end
    end
    % 计算均值点
    for k = 1:K
        center(n,k) = center_list(n,k) ./ count(n,1);
    end
end

end

%% 合成分割后的图片
result = uint8(zeros(I,J,K));
for i = 1:I
    for j = 1:J
        for n = 1:N
            if near_center(i,j) == n
                for k = 1:K
                    result(i,j,k) = center(n,k);
                end
            end
        end
    end
end
imshow(result) % 绘制分割后的图片

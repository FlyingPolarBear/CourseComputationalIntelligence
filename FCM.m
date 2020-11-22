clear all
%% 数据预处理
data = csvread('iris.csv'); % 读取数据
[I,J] = size(data); % I为数据量，J为数据维数
label = data(:,J);
J = J - 1;
data = data(:,1:J);

% PCA降至二维，方便可视化
[coeff,score] = pca(data);
data = score(:,1:2);
J = 2;

%% 设置参数
C = 20;         % 迭代次数
N = max(label); % 聚类数
m = 1.8;        % 模糊度

%% 初始化聚类中心
dis = zeros(I,J);
error = zeros(C,1);
center = min(data) + rand(1,N)' * (max(data) - min(data));

%% 迭代优化
U = zeros(I,N);
for c = 1:C % 迭代C次
    %% 计算隶属度
    for i = 1:I
        for n = 1:N
            dis(i,n) = norm(data(i,:) - center(n,:));
        end
        for n = 1:N
            den = 0;
            for nn = 1:N
                den = den + (dis(i,n)/dis(i,nn))^(2/(m-1));
            end
            U(i,n) = 1 / den;
        end
    end
    
    %% 更新聚类中心
    center = U' * data .* (3/I);
    
    %% 计算误差
    [~,res] = max(U');
    res = res';
    err = 0;
    for i = 1:I
        if (res(i,1) ~= label(i,1))
            err = err + 1;
        end
    end
    error(c,1) = err / I;
    
end

error(C,1)
%% 画图
% 收敛曲线图
figure(1)
plot(error,'Linewidth',2)
title('FCM收敛曲线图')
xlabel('迭代次数')
ylabel('平均误差')
grid on

% 聚类效果散点图
figure(2)
scatter(data(res==1,1),data(res==1,2),80,'r','marker','.')
hold on
scatter(data(res==2,1),data(res==2,2),80,'g','marker','.')
hold on
scatter(data(res==3,1),data(res==3,2),80,'b','marker','.')
hold on
scatter(center(1,1),center(1,2),120,'r','marker','P')
hold on
scatter(center(2,1),center(2,2),120,'g','marker','P')
hold on
scatter(center(3,1),center(3,2),120,'b','marker','P')
hold on
title('FCM聚类效果散点图')
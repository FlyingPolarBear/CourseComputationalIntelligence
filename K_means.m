clear all
%% ����Ԥ����
data = csvread('iris.csv'); %��ȡͼ��
[I,J] = size(data); % IΪ��������JΪ����ά��
label = data(:,J);
data(:,J) = [];

[coeff,score,latent,tsquared] = pca(data);
data = score(:,1:2);
J = 2;

M = 20; %��������
N = max(label); %������
error = zeros(M,1);

% ����ÿһά�ȵ����ֵ����Сֵ,�Լ���ʼ��������
center = zeros(N,1); % ��������
minmax = zeros(J,2); % ����Сֵ�����ֵ
for j = 1:J
    minmax(j,1) = min(data(:,j)); % ��Сֵ
    minmax(j,2) = max(data(:,j)); % ���ֵ
    % �����ʼ����������,ÿһ����һ����ʼ������
    for n = 1:N
        center(n,j) = minmax(j,1) + (minmax(j,2) - minmax(j,1)) * rand();
    end
end

%% �����Ż�
for m = 1:M % ����M��
    %% ��������
    dis = zeros(I,N); % �������
    sum_dis = zeros(I,N); % ���������ɾ���
    grade = zeros(I,N); % �����Ⱦ���
    for i = 1:I
        for n = 1:N
            % �������
            dis(i,n) = exp(-norm(data(i,:)-center(n,:)));
        end
        Dis = sum(dis(i,:));
        for n = 1:N
            grade(i,n) = dis(i,n) / Dis;
        end
    end
    
    %% ���¾�������
    center_list = zeros(n,j);
    count = zeros(n,1);
    for n = 1:N
        center = sum(data*grade(:,n));
    end
    
    %% �������
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

% %% �������ͼ
% plot(error,'LineWidth',2)
% title('K-means����������');
% ylim([0,1]);
% xlabel('��������');
% ylabel('������');
% grid on
clear all
%% ����Ԥ����
data = csvread('iris.csv'); % ��ȡ����
[I,J] = size(data); % IΪ��������JΪ����ά��
label = data(:,J);
J = J - 1;
data = data(:,1:J);

% PCA������ά��������ӻ�
[coeff,score] = pca(data);
data = score(:,1:2);
J = 2;

%% ���ò���
C = 20;         % ��������
N = max(label); % ������
m = 1.8;        % ģ����

%% ��ʼ����������
dis = zeros(I,J);
error = zeros(C,1);
center = min(data) + rand(1,N)' * (max(data) - min(data));

%% �����Ż�
U = zeros(I,N);
for c = 1:C % ����C��
    %% ����������
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
    
    %% ���¾�������
    center = U' * data .* (3/I);
    
    %% �������
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
%% ��ͼ
% ��������ͼ
figure(1)
plot(error,'Linewidth',2)
title('FCM��������ͼ')
xlabel('��������')
ylabel('ƽ�����')
grid on

% ����Ч��ɢ��ͼ
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
title('FCM����Ч��ɢ��ͼ')
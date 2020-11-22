 clear all
%% ����Ԥ����
image = imread('IMAGE1.png'); %��ȡͼ��
[I,J,K] = size(image); % I,JΪͼƬ��С��KΪͨ����
M = 20; %��������
N = 3; %������
% ����ÿһά�ȵ����ֵ����Сֵ,�Լ���ʼ��������
center = zeros(N,1); % ��������
minmax = zeros(K,2); % ����Сֵ�����ֵ
for k = 1:K
    minmax(k,1) = min(min(image(:,:,1))); % ��Сֵ
    minmax(k,2) = max(max(image(:,:,1))); % ���ֵ
    % �����ʼ����������,ÿһ����һ����ʼ������
    for n = 1:N
        center(n,k) = minmax(k,1) + (minmax(k,2) - minmax(k,1)) * rand();
    end
end
image = double(image); % ��ͼ�����ת��Ϊdouble����

%% �����Ż�
for m = 1:M % ����M��
    m
%% ��������
dis = zeros(I*J,N); % �������
near = zeros(I,J); % ���ھ���
near_center = zeros(I,J); % �������ľ���
near_center1 = zeros(I*J,1); % �����ľ������ľ���
for i = 1:I
    for j = 1:J
        % ÿ����������
        for n = 1:N
            dis((i-1)*J+j,n) = sqrt((image(i,j,1)-center(n,1))^2+...
            (image(i,j,2)-center(n,2))^2+(image(i,j,3)-center(n,3))^2); 
        end
        % �ҵ�����ľ�������
        [near(i,j),near_center1((i-1)*J+j,1)] = min(dis((i-1)*J+j,:),[],2);
    end
end
for i = 1:I
    for j = 1:J
        near_center(i,j) = near_center1((i-1)*J+j,1); % ����
    end
end

%% ���¾�������
center_list = zeros(n,k);
count = zeros(n,1);
for i = 1:I
    for j = 1:J
        for n = 1:N
            % ����ǽ��ڣ������list
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
    % ������ڿ���أ������·���
    if (M - m > 2) && count(n,1) == 0
        for k = 1:K
            center(n,k) = minmax(k,1) + (minmax(k,2) - minmax(k,1)) * rand();
        end
    end
    % �����ֵ��
    for k = 1:K
        center(n,k) = center_list(n,k) ./ count(n,1);
    end
end

end

%% �ϳɷָ���ͼƬ
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
imshow(result) % ���Ʒָ���ͼƬ

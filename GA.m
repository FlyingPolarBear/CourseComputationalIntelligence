clear all;
% �Ŵ��㷨����Ԫ�������ֵ

%% ��ʼ������
% GA����
T = 100;  % �������
N = 100;  % Ⱥ���ģ
pm = 0.05;% �������
pc = 0.8; % �������

% ����ȡֵ��Χ
umax1 = 12.1; % x1���ֵ
umin1 = -3.0; % x1��Сֵ
umax2 = 5.8; % x2���ֵ
umin2 = 4.1; % x2��Сֵ

% ��ʼ���������Ӧ��ֵ
L = 16;                     % ���������ִ����ȣ��ܱ��볤��2L
bval = round(rand(N,2*L));  % �����ʼ����Ⱥ
bestv = -inf;               % ������Ӧ�ȳ�ֵ

%% ������ʼ
for ii=1:T
    
    % ���룬������Ӧ��
    for i=1:N
        y1 = 0;
        y2 = 0;
        for j = 1:1:L
            y1 = y1+bval(i,L-j+1)*2^(j-1);
        end
        x1 = (umax1-umin1)*y1/(2^L-1)+umin1;            % x1����
        for j = 1:1:L
            y2 = y2+bval(i,2*L-j+1)*2^(j-1);
        end
        x2 = (umax2-umin2)*y2/(2^L-1)+umin2;            % x2����
        obj(i) = 21.5+x1*sin(4*pi*x1)+x2*sin(20*pi*x2); %Ŀ�꺯��
        xx(i,:) = [x1,x2];
    end
    func = obj;                 % Ŀ�꺯��ת��Ϊ��Ӧ�Ⱥ���
    p = func./sum(func);
    q = cumsum(p);              % �ۼ�
    [fmax,indmax] = max(func);	% �󵱴���Ѹ���
    if fmax >= bestv
        bestv = fmax;           % ��ĿǰΪֹ������Ӧ��ֵ
        bvalxx = bval(indmax,:);% ��ĿǰΪֹ���λ��
        optxx = xx(indmax,:);   % ��ĿǰΪֹ���Ų���
    end
    Bfit1(ii)=bestv;            % �洢ÿ����������Ӧ��
    
    % �Ŵ�������ʼ
    % 1.���̶�ѡ��
    for i=1:(N-1)
        r=rand;
        tmp=find(r<=q);
        newbval(i,:)=bval(tmp(1),:);
    end
    newbval(N,:)=bvalxx;  % ���ű���
    bval=newbval;
    % 2.���㽻��
    for i=1:2:(N-1)
        cc=rand;
        if cc<pc
            point=ceil(rand*(2*L-1)); % ȡ��һ��1��2L-1������
            ch=bval(i,:);
            bval(i,point+1:2*L)=bval(i+1,point+1:2*L);
            bval(i+1,point+1:2*L)=ch(1,point+1:2*L);
        end
    end
    bval(N,:)=bvalxx; % ���ű�������
    % 3.λ�����
    mm=rand(N,2*L)<pm; % С�ڱ�����ʸ�ֵΪ1����������Ϊ0
    mm(N,:)=zeros(1,2*L); % ���һ�в����죬ǿ�Ƹ�0
    bval(mm)=1-bval(mm);
end

%% ���
plot(Bfit1,'Linewidth',2)
title('�Ŵ��㷨����Ԫ�������ֵ����')
xlabel('�������')
ylabel('������Ӧ��ֵ')
bestv	% ���������Ӧ��ֵ
optxx	% ������ű���ȡֵ
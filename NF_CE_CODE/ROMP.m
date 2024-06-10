function [ theta ] =ROMP( y,A,K )
%CS_ROMP Summary of this function goes here
%Version: 1.0 written by jbb0523 @2015-04-24
%   Detailed explanation goes here
%   y = Phi * x
%   x = Psi * theta
%	y = Phi*Psi * theta
%   令 A = Phi*Psi, 则y=A*theta
%   现在已知y和A，求theta
%   Reference:Needell D，Vershynin R．Signal recovery from incomplete and
%   inaccurate measurements via regularized orthogonal matching pursuit[J]．
%   IEEE Journal on Selected Topics in Signal Processing，2010，4(2)：310—316.
    [y_rows,y_columns] = size(y);
    if y_rows<y_columns
        y = y';%y should be a column vector
    end
    [M,N] = size(A);%传感矩阵A为M*N矩阵
    theta = zeros(N,1);%用来存储恢复的theta(列向量)
    At = zeros(M,3*K);%用来迭代过程中存储A被选择的列
    Pos_theta = zeros(1,2*K);%用来迭代过程中存储A被选择的列序号
    Index = 0;
    r_n = y;%初始化残差(residual)为y
    %Repeat the following steps K times(or until |I|>=2K)
    for ii=1:K%迭代K次
        product = A'*r_n;%传感矩阵A各列与残差的内积
        %[val,pos] = max(abs(product));%找到最大内积绝对值，即与残差最相关的列
        [val,pos] = Regularize(product,K);%按正则化规则选择原子
        At(:,Index+1:Index+length(pos)) = A(:,pos);%存储这几列
        Pos_theta(Index+1:Index+length(pos)) = pos;%存储这几列的序号
        if Index+length(pos)<=M%At的行数大于列数，此为最小二乘的基础(列线性无关)
            Index = Index+length(pos);%更新Index，为下次循环做准备
        else%At的列数大于行数，列必为线性相关的,At(:,1:Index)'*At(:,1:Index)将不可逆
            break;%跳出for循环
        end
        A(:,pos) = zeros(M,length(pos));%清零A的这几列(其实此行可以不要因为它们与残差正交)
        %y=At(:,1:Index)*theta，以下求theta的最小二乘解(Least Square)
        theta_ls = (At(:,1:Index)'*At(:,1:Index))^(-1)*At(:,1:Index)'*y;%最小二乘解
        %At(:,1:Index)*theta_ls是y在At(:,1:Index)列空间上的正交投影
        r_n = y - At(:,1:Index)*theta_ls;%更新残差
        if norm(r_n)<1e-6%Repeat the steps until r=0
            break;%跳出for循环
        end
        if Index>=2*K%or until |I|>=2K
            break;%跳出for循环
        end
    end
    theta(Pos_theta(1:Index))=theta_ls;%恢复出的theta
end 
function [rs] = LAR_f(y,fi,Oinitial,Oupdata,k)
%Oupdata是最大的L个中的一个，外层循环每次带一个L中的一个进来检测

ss = length(Oinitial);
Oinitial(ss) = Oupdata;%和并初始下标索引集合和新选出的L个最大下标
Onew = Oinitial;
s = length(Onew);%新下标集合长度
[M,N] = size(fi);
fio = zeros(M,k);
sss = length(Onew);
fio(:,1:sss) = fi(:,Onew);
x = (fio(:,1:sss)'*fio(:,1:sss))^(-1)*fio(:,1:sss)'*y;%pinv(fio)求伪逆
r_n = y - fio(:,1:sss)*x;%初始残差
    for i = (s + 1) : k
        product = fi'*r_n;
        [val,pos] = max(abs(product));
        fio(:,i) = fi(:,pos);%存储这一列
        Onew(i) = pos;%存储这一列的序号
        theta_ls = (fio(:,1:i)'*fio(:,1:i))^(-1)*fio(:,1:i)'*y;
        r_n = y - fio(:,1:i)*theta_ls;
    end
    rs = norm(r_n);
end
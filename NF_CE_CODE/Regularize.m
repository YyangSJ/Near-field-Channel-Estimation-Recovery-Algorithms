function [val,pos] = Regularize(product,Kin)
%Regularize Summary of this function goes here
%   Detailed explanation goes here
%   product = A'*r_n;%传感矩阵A各列与残差的内积
%   K为稀疏度
%   pos为选出的各列序号
%   val为选出的各列与残差的内积值
%   Reference:Needell D，Vershynin R. Uniform uncertainty principle and
%   signal recovery via regularized orthogonal matching pursuit. 
%   Foundations of Computational Mathematics, 2009,9(3): 317-334.  
    productabs = abs(product);%取绝对值
    [productdes,indexproductdes] = sort(productabs,'descend');%降序排列
    for ii = length(productdes):-1:1
        if productdes(ii)>1e-6%判断productdes中非零值个数
            break;
        end
    end
    %Identify:Choose a set J of the K biggest coordinates
    if ii>=Kin
        J = indexproductdes(1:Kin);%集合J
        Jval = productdes(1:Kin);%集合J对应的序列值
        K = Kin;
    else%or all of its nonzero coordinates,whichever is smaller
        J = indexproductdes(1:ii);%集合J
        Jval = productdes(1:ii);%集合J对应的序列值
        K = ii;
    end
    %Regularize:Among all subsets J0∈J with comparable coordinates
    MaxE = -1;%循环过程中存储最大能量值
    for kk = 1:K
        J0_tmp = zeros(1,K);iJ0 = 1;
        J0_tmp(iJ0) = J(kk);%以J(kk)为本次寻找J0的基准(最大值)
        Energy = Jval(kk)^2;%本次寻找J0的能量
        for mm = kk+1:K
            if Jval(kk)<2*Jval(mm)%找到符合|u(i)|<=2|u(j)|的
                iJ0 = iJ0 + 1;%J0自变量增1
                J0_tmp(iJ0) = J(mm);%更新J0
                Energy = Energy + Jval(mm)^2;%更新能量
            else%不符合|u(i)|<=2|u(j)|的
                break;%跳出本轮寻找，因为后面更小的值也不会符合要求
            end
        end
        if Energy>MaxE%本次所得J0的能量大于前一组
            J0 = J0_tmp(1:iJ0);%更新J0
            MaxE = Energy;%更新MaxE，为下次循环做准备
        end
    end
    pos = J0;
    val = productabs(J0);
end 
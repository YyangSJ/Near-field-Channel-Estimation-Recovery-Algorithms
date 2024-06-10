function [pos_array,hat_x]=cs_somp(y,T_Mat,G,L)
 

[m,n]=size(y);
s=L; 

hat_x=zeros(G,n);                 
Aug_t=[];       
r_n=y;  

for n=1:s;  
    for i=1:G
        pro(i,:)=(T_Mat(:,i))'/abs(norm(T_Mat(:,i)))*r_n;
    end
   % pro=T_Mat'*r_n; 
    for i=1:G
       product(i)=sum(abs(pro(i,:)));
    end
   
    [val,pos]=max(product);   %最大投影系数对应的位置
        x_ind(n) = mod(pos - 1, sqrt(G)) + 1;
    z_ind(n) = ceil(pos / sqrt(G));
    Aug_t=[Aug_t,T_Mat(:,pos)];   %矩阵扩充
   % T_Mat(:,pos)=zeros(m,1); %选中的列置零
    aug_x=(Aug_t'*Aug_t+1e-7*eye(n))^(-1)*Aug_t'*y;   
  %   aug_x=pinv(Aug_t)*y;  
    r_n=y-Aug_t*aug_x;   %残差 
    pos_array(n)=pos;   %纪录最大投影系数的位置
    
end

hat_x(pos_array,:)=aug_x;  %  重构的向量 


end
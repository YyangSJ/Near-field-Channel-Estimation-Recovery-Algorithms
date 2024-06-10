function [ theta ] = LAOMP_f( y,A,t,L )
    [y_rows,y_columns] = size(y);
    if y_rows<y_columns
        y = y';
    end
    [M,N] = size(A);
    theta = zeros(N,1);
    At = zeros(M,t);
    Pos_theta = zeros(1,t);
    r_n = y;
    for ii=1:t
        sim = A'*r_n;
        [val1,pos1]=sort(abs(sim),'descend');%pos降序排列下标
        jv = pos1(1:L);
        [jj,ij,io] = intersect(jv,Pos_theta);
        jv(ij) = [];
        rs = zeros(1,length(jv));
        pos_theta = Pos_theta(:,1:ii);
        for i = 1 : length(jv)
            dex = jv(i);
            rs(i) = norm(LAR_f(y,A,pos_theta,dex,t));
        end
        [val,l] = min(rs);%l用来迭代过程中存储fi被选择的列序号
        pos = pos1(l);
        At(:,ii) = A(:,pos);
        Pos_theta(ii) = pos;
        theta_ls = (At(:,1:ii)'*At(:,1:ii))^(-1)*At(:,1:ii)'*y;%最小二乘解比伪逆好
        r_n = y - At(:,1:ii)*theta_ls;
    end
    theta(Pos_theta)=theta_ls;
end
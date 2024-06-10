function Result=FMFSBL(A,y,Max_iter)
    [M,N]=size(A);
    alpha=ones(N,1);
    mu=zeros(N,1);
    
 % mu=A'*inv(A*A')*y;  mu=A'*y;  bad performance
  
    a = 1e-6;
    b = 1e-6;
    c = 1e-6;
    d = 1e-6;
    S = svd(A,'econ');
    L = 2*S(1)^2;%+0*0.00001; % Lipschitz constant
    delta = 1;
    theta = mu;
    ATY = A'*y;
    ATA = A'*A;
    ba=diag(ATA);
    for iter = 1 : Max_iter
        mu_old = mu;
        mu=inv(eye(N)*L*delta/2+diag(alpha))*(L/2*delta*mu-delta*A'*(A*mu-y));
        for n=1:N
            sigma(n)=1/(delta*ba(n)+alpha(n));
            beta(n)=b+0.5*abs(mu(n))^2+0.5*sigma(n);
        end  
          
        lambda=d+0.5*norm(y-A*mu)^2+0.5*ba'*sigma.';
        
        for n=1:N
            alpha(n)=(0.5+a)/beta(n);
        end
        delta=(d+M/2)/lambda;
%         if norm(mu-mu_old)/norm(mu)<1e-8
%             break
%         end
    end
    Result = mu;
end
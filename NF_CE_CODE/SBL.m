function [mu] = SBL(y,PHI,II)

sig2=1;
G=size(PHI,2);
Lambda=eye(G);
rho=1*ones(G,1);

for iter=1:II
    SIG=inv(sig2^(-1)*PHI'*PHI+Lambda);
    mu=sig2^(-1)*SIG*PHI'*y;
    for i=1:G
        chi(i)=1-rho(i)*SIG(i,i);
        rho(i)=chi(i)/(abs(mu(i))+1e-10)^2;
    end
    sig2=norm(y-PHI*mu)^2/(G-sum(chi));
    Lambda=diag(rho);
end

end


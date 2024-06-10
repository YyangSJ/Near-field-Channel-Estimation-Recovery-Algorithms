function [mu] = V_SBL(y,PHI,II)

Esig2=1;
G=size(PHI,2);
Lambda=eye(G);
 

    a = 1e-6;
    b = 1e-6;
    c = 1e-6;
    d = 1e-6;
 
for iter=1:II
    SIG=inv(Esig2*PHI'*PHI+Lambda);
    mu=Esig2*SIG*PHI'*y;
    for i=1:G
       
        rho(i)=(1/2+a)/((abs(mu(i))^2+SIG(i,i))/2+b);
    end
    Esig2=(c+G/2)/((norm(y)^2-2*real(y'*PHI*mu)+trace(PHI'*PHI*(mu*mu'+SIG)))/2+d);
    Lambda=diag(rho);
end

end


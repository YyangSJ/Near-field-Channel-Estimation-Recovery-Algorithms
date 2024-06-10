function [mu] = MFV_SBL(y,PHI,PHI_co,II)

Esig2=1;
[M,G]=size(PHI);
Lambda=eye(G);


a = 1e-6;
b = 1e-6;
c = 1e-6;
d = 1e-6;

Py=PHI'*y;


 
 mu=zeros(G,1);
  mu=PHI'*inv(PHI*PHI')*y;

% for i = 1:G
%     SIG(i,i)=1/(Esig2*PHI_co(i,i)+Lambda(i,i));
%     mu(i,1)  = Esig2*PHI(:,i)'*(y)*SIG(i,i);
% end 
for iter=1:II
    
    for i=1:G
        %         PHI_t=PHI;
        %        PHI_t(:,i)=[];
        mu_t=mu;
        mu_t(i)=[];
        PHI_t=PHI_co(i,:);
        PHI_t(:,i)=[];
        SIG(i,i)=1/(Esig2*PHI_co(i,i)+Lambda(i,i));
        mu_t1(i)=Esig2*SIG(i,i)*PHI(:,i)'*(y-PHI*mu+PHI(:,i)*mu(i));
        cmu(i)=Esig2*SIG(i,i)*(Py(i)-PHI_t*mu_t);
        % mu(i)=Esig2*SIG(i,i)*PHI(:,i)'*(y-PHI_t*mu_t)
        rho(i)=(1/2+a)/((abs(mu(i))^2+SIG(i,i))/2+b);
        %Lambda(i,i)=(1/2+a)/((abs(mu(i))^2+SIG(i,i))/2+b);
    end
    
   % Esig2=(c+M/2)/((norm(y)^2-2*real(Py'*mu)+trace(PHI_co*(mu*mu'+SIG)))+d);
    
   Esig2=(c+M/2)/((norm(y)^2-2*real(Py'*mu)+mu'*PHI_co*mu+sum(diag(PHI_co).*diag(SIG)))/2+d);
    
    Lambda=diag(rho);
    
end

end


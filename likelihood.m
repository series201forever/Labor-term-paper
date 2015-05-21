function [ likelihood] = likelihood(para,data,n,exemp,tau)
%Compute all threshold epsilons using emax function
[epsstarv,epsstarstarv] = emax(para,data,exemp,tau);

%Construct likelihood
%Foreach observed education,race and experience, pick the two threshold
%computed before, compare them, and integrate the probability over the
%appropriate range.
    epsstarlik=zeros(1000,15);
lik=zeros(1000,15);
for i=1:n
s=data(1+15*(i-1),7);
r=data(1+15*(i-1),8);
    for t=1:15
        wobs=data(1+15*(i-1)+(15-t),4);
        eobs=data(1+15*(i-1)+(15-t),6);
        qobs=eobs+1-data(1+15*(i-1),6);
epsstarlik(i,(16-t))=epsstarstarv(i,qobs,(16-t))-epsstarv(i,qobs,(16-t));
        if (epsstarlik(i,(16-t))<0)&&(data(1+15*(i-1)+(15-t),3)==1)
           lik(i,(16-t))=log(max(10^(-100),(1-normcdf((epsstarv(i,qobs,(16-t))-(para(7)^2/(para(7)^2+para(8)^2))*(log(wobs)-para(1)-para(2)*s-para(3)*eobs-para(4)*eobs*r-para(9)*eobs^2))/(para(7)*sqrt(1-(para(7)^2/(para(7)^2+para(8)^2))))))...
            *normpdf(log(wobs)-para(1)-para(2)*s-para(3)*eobs-para(4)*eobs*r-para(9)*eobs^2,0,(sqrt(para(7)^2+para(8)^2)))));
        else if  (epsstarlik(i,(16-t))<0)&&(data(1+15*(i-1)+(15-t),3)==0)
              lik(i,(16-t))=log(max(10^(-100),normcdf(epsstarv(i,qobs,(16-t))/para(7))));
            else if (epsstarlik(i,(16-t))>0)&&(data(1+15*(i-1)+(15-t),3)==1)
               lik(i,(16-t))=log(max(10^(-100),(1-normcdf((epsstarstarv(i,qobs,(16-t))-(para(7)^2/(para(7)^2+para(8)^2))*(log(wobs)-para(1)-para(2)*s-para(3)*eobs-para(4)*eobs*r-para(9)*eobs^2))/(para(7)*sqrt(1-(para(7)^2/(para(7)^2+para(8)^2))))))...
           *normpdf(log(wobs)-para(1)-para(2)*s-para(3)*eobs-para(4)*eobs*r-para(9)*eobs^2,0,(sqrt(para(7)^2+para(8)^2)))));
                else if (epsstarlik(i,(16-t))>0)&&(data(1+15*(i-1)+(15-t),3)==0)
              lik(i,(16-t))=log(max(10^(-100),normcdf(epsstarstarv(i,qobs,(16-t))/para(7))));
                    else 
                        lik(i,(16-t))=-1000;                         
                    end
                end
             end
        end    
    
    end
end


likelihood=-sum(sum(lik(:,:)));

end


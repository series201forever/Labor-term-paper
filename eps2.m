function [ eps2] = eps2(s,e,r,yh,paramat,exemp,tau,t,emaxv)
eps2=log(max(10^(-100),(paramat{5}+paramat{6}.*yh+0.95*(emaxv(:,(1:15),(16-t)+1)-emaxv(:,(2:16),(16-t)+1))-(exemp.*tau))./(exp(paramat{1}+paramat{2}.*s+paramat{3}.*e+paramat{4}.*e.*r+paramat{9}.*e.^2).*(1-tau))));

end


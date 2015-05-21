function [ eps1] = eps1(s,e,r,yh,paramat,t,emaxv)
eps1=log(max(10^(-100),(paramat{5}+paramat{6}.*yh+0.95*(emaxv(:,(1:15),(16-t)+1)-emaxv(:,(2:16),(16-t)+1)))./exp(paramat{1}+paramat{2}.*s+paramat{3}.*e+paramat{4}.*e.*r+paramat{9}.*e.^2)));

end


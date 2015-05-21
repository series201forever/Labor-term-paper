function [epsstarv,epsstarstarv ] = emax(para,data,exemp,tau)
%To be revised: Matlab is column major. Use first argument to construct
%loop. i.e. set t as first argument.

%This function computes two threshold values determining working/not. The
%inputs are parameters to be estimated, data, amount of tax exemption and
%tax rate.

%For each period. I calculate both thresholds in number of
%individual*experience matrix.
%prepare parameters for matrix notations
paraa=para(1)*ones(1000,15);
parab=para(2)*ones(1000,15);
parac=para(3)*ones(1000,15);
parad=para(4)*ones(1000,15);
parae=para(5)*ones(1000,15);
paraf=para(6)*ones(1000,15);
parag=para(7)*ones(1000,15);
parah=para(8)*ones(1000,15);
parai=para(9)*ones(1000,15);
paramat=[{paraa},{parab},{parac},{parad},{parae},{paraf},{parag},{parah},{parai}];
exempmat=exemp*ones(1000,15);
taumat=tau*ones(1000,15);

%Construct matrices to store variables
%1000*15 (people*experience)*years
epsstarv=zeros(1000,15,15);
epsstarstarv=zeros(1000,15,15);
emaxv=zeros(1000,16,16);

%Pick an individual (set of initial condition i.e. race, education,
%husbands income. Make everything 1000*15 for each period.
s=data(1:15:end,7)*ones(1,15);
r=data(1:15:end,8)*ones(1,15);
%Solve DP from period 15. t=1 implies period 15.
for t=1:15
    yh=data(((16-t):15:end),5)*ones(1,15);
    %For each period, possible accumulation of experience is between zero
    %and 14, apparently if t>2 this is redundant though.
    q=ones(1000,1)*(1:15);
    e=data(1:15:end,6)*ones(1,15)+(q-1);
    wagefixed=exp(paramat{1}+paramat{2}.*s+paramat{3}.*e+paramat{4}.*e.*r+paramat{9}.*e.^2);
    %Compute thresholds.
    epsstar=eps1(s,e,r,yh,paramat,t,emaxv);
    epsstarstar=eps2(s,e,r,yh,paramat,exempmat,taumat,t,emaxv);
    epscut=log(exemp./wagefixed);
    
    epsstarv(:,:,(16-t))=epsstar;
    epsstarstarv(:,:,(16-t))=epsstarstar;


%To compute previous period, I need emax function corresponds to the
%current period. This is going to be stored in emaxv matrix, to be referred
%by eps1 and eps2 function in the previous period.
if epsstar>epsstarstar
    emaxv(:,(1:15),(16-t))=(yh.*(1+paramat{6})+paramat{5}).*normcdf(epsstar,0,para(7))...
        +(yh+wagefixed).*(normcdf(epscut,0,para(7))-normcdf(epsstar,0,para(7)))...
        +exp((paramat{7}.^2)/2).*(normcdf(paramat{7}-(epsstar./paramat{7}))-normcdf(paramat{7}-(epscut./paramat{7})))...
    +(yh+exemp.*tau+wagefixed.*(1-tau)).*(1-normcdf(epscut,0,para(7)))+(1-tau).*exp((paramat{7}.^2)/2).*normcdf(paramat{7}-(epsstarstar./paramat{7}));
elseif epsstar<epsstarstar
     emaxv(:,(1:15),(16-t))=(yh.*(1+paramat{6})+paramat{5}).*normcdf(epsstarstar,0,para(7))...
        +(yh+exemp.*tau+wagefixed.*(1-tau)).*(1-normcdf(epsstarstar,0,para(7)))+(1-tau).*exp((paramat{7}.^2)/2).*normcdf(paramat{7}-(epsstarstar./paramat{7}));
    
end
 
end
end


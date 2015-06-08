%This code is to estimate parameters using simulated method of moments

%Set initial parameters
ucurve = 0.6;
vfir   = 4;
vsec   = 0.5;
cint   = 25;
cslope = 5;    
cageint   = -0.3;
cageslope = 0.3;
conskid=0.7;
mum     = 2.83;
muf     = 2.503;
sigmam  = 0.14;
sigmaf  = 0.135;
Gamma = 0.09;

%These are the parameters estimated from SIPP data separately.
deltam  = 0.032;
deltaf  = 0.05;
lambdam = 0.4;
lambdaf = 0.3;


%Beta cannot be identified. So fix it to some value.
Beta    = 0.91;  

parainit = [ucurve,vfir,vsec,cint,cslope,cageint,cageslope,...
    conskid,mum,muf,sigmam,sigmaf,Gamma];
parafixed = [deltam,deltaf,lambdam,lambdaf,Beta];


%Minimization
obj=@(para)smmobjective(para,parafixed );

options = optimoptions(@fmincon,'Display','iter','MaxFunEvals',1e9,'MaxIter',1e9,'TolFun',1e-2,'TolX',1e-8); 
%A=[0 0 0 0 1 -1 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 1 -1];
%b=[0;0];
lb=[0.01, 0,   0.01, -50,  0,  -50, -20, 0,  0.01, 0.01, 0.01, 0.01, 0.001];
ub=[0.99, 100, 0.99,  50,  20, 50,  20,  10, 10,    10,   10,   10,  0.99];
[soltest4,fval,exitflag,output,lambda,grad,hessiantest4]=fmincon(obj,parainit,[],[],[],[],lb,ub,[],options);

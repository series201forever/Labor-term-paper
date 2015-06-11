%This code is to estimate parameters using simulated method of moments

%TEST, NO unemployment rate, no variance of wage
%Parameter of wage variance fixed for now.


%Set initial parameters
ucurve = 0.6;
vfir   = 4;
vsec   = 0.4;
cint   = 25;
cslope = 5;    
cageint   = -0.3;
cageslope = 0.3;
conskid=0.7;
mum     = 3.2140;
muf     =2.7134;
sigmam  = 0.14;
sigmaf  = 0.135;
Gamma = 0.05;

%These are the parameters estimated from SIPP data separately.
deltam  = 0.02;
deltaf  = 0.05;
lambdam = 0.4;
lambdaf = 0.3;


%Beta cannot be identified. So fix it to some value.
Beta    = 0.98;  

  parainit = [ucurve,vfir,vsec,cint,cslope,cageint,cageslope,...
    conskid,mum,muf,Gamma];
parafixed = [deltam,deltaf,lambdam,lambdaf,Beta,sigmam,sigmaf];


%Minimization
obj=@(para)smmobjective(para,parafixed );

options = optimoptions(@fmincon, 'Algorithm','interior-point','Display','iter','MaxFunEvals',1e9,'MaxIter',20,'TolFun',1e-2,'TolX',1e-8, 'Diffminchange',1e-2); 
%A=[0 0 0 0 1 -1 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 1 -1];
%b=[0;0];
lb=[0.01, 0,   0.01, -50,  0,  -50, -5, 0,  2.6, 1.8, 0.01];
ub=[0.99, 100, 0.99,  50,  20, 50,  5,  10, 3.7,    3.2,  0.7];

[soltest4,fval,exitflag,output,lambda,grad,hessiantest4]=fmincon(obj,parainit,[],[],[],[],lb,ub,[],options);


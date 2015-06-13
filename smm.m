%This code is to estimate parameters using simulated method of moments

%TEST, NO unemployment rate, no variance of wage
%Parameter of wage variance fixed for now.

%Set initial parameters
ucurve = 0.5018;
vfir   = 4.1059;
vsec   = 0.2988;
cint   = 35.9100;
cslope = 10.3759;    
cageint   = 0.1534;
cageslope = 0.2045;
conskid=0.6655;
mum     = 3.2778;
muf     =2.7616;
sigmam  = 0.8195;
sigmaf  = 0.6253;
Gamma = 0.2492;

%These are the parameters estimated from SIPP data separately.
deltam  = 0.1287;
deltaf  = 0.1416;
lambdam = 0.3804;
lambdaf = 0.2;


%Beta cannot be identified. So fix it to some value.
Beta    = 0.98; 

  parainit = [ucurve,vfir,vsec,cint,cslope,cageint,cageslope,...
    conskid,mum,muf,sigmam,sigmaf,Gamma];
parafixed = [deltam,deltaf,lambdam,lambdaf,Beta];


%Minimization
obj=@(para)smmobjective(para,parafixed );

options = optimoptions(@fmincon, 'Algorithm','interior-point','Display','iter','MaxFunEvals',1e9,'MaxIter',1,'TolFun',1e-1,'TolX',1e-4, 'Diffminchange',1e-2); 
%A=[0 0 0 0 1 -1 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 1 -1];
%b=[0;0];
lb=[0.1, 0,   0.01,  0,  0,  -50, -5, 0,  2.6, 1.8, 0.525, 0.6, 0.01];
ub=[0.8, 100, 0.99,  50,  20, 50,  5,  10, 3.7  3.2, 1.2, 1.2, 0.7];

[soltest4,fval,exitflag,output,lambda,grad,hessiantest4]=fmincon(obj,parainit,[],[],[],[],lb,ub,[],options);

%%
a=15620.316
b=125520428.8881
sigma=sqrt(log(b/(a^2)+1))


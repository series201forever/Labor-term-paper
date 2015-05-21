clear all;
rng(278);
data=importdata('PS3_data.csv');
%To be developed; tabulation can be done in stata.

%Set hypothetical para
exemp=3;
tau=0.5;
n=1000;
partest=[ 0.7513    0.0163    0.0165    0.0009    0.3847    0.3012    0.3982    0.1035   -0.0002]
%[g1,g2,g3,g4,a1,a2,sigmaeps,sigmaeta,g5];

%Check the function constructed
%[test2a,test2b]=emax(partest,data,exemp,tau)
%test2=likelihood(partest,data,n,exemp,tau)




%%

%Maximum likelihood
obj=@(para)likelihood(para,data,n,exemp,tau);
par0=partest;
lb=[-Inf,0.001,0.001,-Inf,0.01,0.01,0.001,0.001,-Inf];
ub=[Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf];
options=optimset('Display','iter','MaxFunEvals',1e9,'MaxIter',1e9,'TolFun',1e-04,'TolX',1e-10); 
[sol,fval,exitflag,output,lambda,grad,hessian]=fmincon(obj,par0,[],[],[],[],lb,ub,[],options)

%%
%Compute std error
var=inv(hessian);
stderror=sqrt(diag(var));


%%
%Simulating samples. i.e. draw 50 wage shock and 50 obs error for each
%person i, and for each draw compute optimal sequence of experience and
%wages. Here I computed sequence of experience using q (=experience
%accumulated in the sample period).

parahat=sol;
[epsstarhat,epsstarstarhat ] = emax(parahat,data,exemp,tau);
nind=1000;
nsim=100;
[q,works,wages]=simul(parahat,epsstarhat,epsstarstarhat,nind,nsim,data);
simdata =simdataset( data,q,works,wages,nsim );


%%
%Modify original dataset to add in-sample experience
origexp=kron(data(1:15:end,6),ones(15,1));
expinsample=data(:,6)-origexp;
obsdata=[data,expinsample];



%%
%Check the difference
%Average number of period overall/each educ level, and mean wage.
workmeandata=mean(obsdata(15:15:end,9))+mean(obsdata(15:15:end,3));
workmeansim=mean(simdata(15:15:end,9))+mean(simdata(15:15:end,3));
%To compute wage, I need to use only the part where I
%observe people working. For simulated data however, since for every period I observe at
%least one sequence where he works, I don't need to take subset of the
%data.
fulldatawork=obsdata(find(obsdata(:,3)>0),:);
meanwagedata=mean(fulldatawork(:,4));
meanwagesim=mean(simdata(:,4));


%For subset of samples, first I use all the data from each subset to
%compute length of working periods.
dataeducless12=obsdata(find(obsdata(:,7)<12),:);
simeducless12=simdata(find(simdata(:,7)<12),:);
less12workmeandata=mean(dataeducless12(15:15:end,9))+mean(dataeducless12(15:15:end,3));
less12workmeansim=mean(simeducless12(15:15:end,9))+mean(simeducless12(15:15:end,3));
%Then I reduce the sample further to pick up only working periods so that
%I can calculate mean wage.
dataeducless12full=dataeducless12(find(dataeducless12(:,3)>0),:);
less12wagemeandata=mean(dataeducless12full(:,4));
less12wagemeansim=mean(simeducless12(:,4));

%Do the same for other subsamples.
dataeducat12=obsdata(find(obsdata(:,7)==12),:);
simeducat12=simdata(find(simdata(:,7)==12),:);
at12workmeandata=mean(dataeducat12(15:15:end,9))+mean(dataeducat12(15:15:end,3));
at12workmeansim=mean(simeducat12(15:15:end,9))+mean(simeducat12(15:15:end,3));
dataeducat12full=dataeducat12(find(dataeducat12(:,3)>0),:);
at12wagemeandata=mean(dataeducat12full(:,4));
at12wagemeansim=mean(simeducat12(:,4));

dataeducbet1315=obsdata(find(obsdata(:,7)>12&obsdata(:,7)<16),:);
simeducbet1315=simdata(find(simdata(:,7)>12&simdata(:,7)<16),:);
bet1315workmeandata=mean(dataeducbet1315(15:15:end,9))+mean(dataeducbet1315(15:15:end,3));
bet1315workmeansim=mean(simeducbet1315(15:15:end,9))+mean(simeducbet1315(15:15:end,3));
dataeducbet1315full=dataeducbet1315(find(dataeducbet1315(:,3)>0),:);
bet1315wagemeandata=mean(dataeducbet1315full(:,4));
bet1315wagemeansim=mean(simeducbet1315(:,4));

dataeducmore16=obsdata(find(obsdata(:,7)>15),:);
simeducmore16=simdata(find(simdata(:,7)>15),:);
more16workmeandata=mean(dataeducmore16(15:15:end,9))+mean(dataeducmore16(15:15:end,3));
more16workmeansim=mean(simeducmore16(15:15:end,9))+mean(simeducmore16(15:15:end,3));
dataeducmore16full=dataeducmore16(find(dataeducmore16(:,3)>0),:);
more16wagemeandata=mean(dataeducmore16full(:,4));
more16wagemeansim=mean(simeducmore16(:,4));

%Fraction of women working at each age, and mean wage
fracworktdata=zeros(1,15);
fracworktsim=zeros(1,15);
meanwagetdata=zeros(1,15);
meanwagetsim=zeros(1,15);
for t=1:15
    workt=obsdata(t:15:end,3);
    fracworktdata(t)=mean(workt);
    wagetaux=obsdata(t:15:end,4);
    waget=wagetaux(wagetaux>0);
    meanwagetdata(t)=mean(waget);
    worktsim=simdata(t:15:end,3);
    fracworktsim(t)=mean(worktsim);
    wagetsimaux=simdata(t:15:end,4);
    wagetsim=wagetsimaux(wagetsimaux>0);
    meanwagetsim(t)=mean(wagetsim);
end

%Fraction of women working by experience levels
worke10data=obsdata(find(obsdata(:,6)<=10),:);
worke1120data=obsdata(find(obsdata(:,6)>10&obsdata(:,6)<=20),:);
worke21data=obsdata(find(obsdata(:,6)>20),:);
worke10sim=simdata(find(simdata(:,6)<=10),:);
worke1120sim=simdata(find(simdata(:,6)>10&simdata(:,6)<=20),:);
worke21sim=simdata(find(simdata(:,6)>20),:);

fracworke10data=mean(worke10data(:,3));
fracworke1120data=mean(worke1120data(:,3));
fracworke21data=mean(worke21data(:,3));
fracworke10sim=mean(worke10sim(:,3));
fracworke1120sim=mean(worke1120sim(:,3));
fracworke21sim=mean(worke21sim(:,3));

%%

input3.data = [fracworke10data fracworke1120data fracworke21data;fracworke10sim fracworke1120sim fracworke21sim]
input3.transposeTable = 1;
latex1 = latexTable(input3)
%%
%Q3
%Solve for counterfactual
%Baseline: original specification
nind=1000;
nsim=100;
%Under the new threshold, simulate the economy again, this time without measurement error.
[q0,works0,wages0]=simul2(parahat,epsstarhat,epsstarstarhat,nind,nsim,data);

%Construct simulated samples as before.
simdata0 =simdataset( data,q0,works0,wages0,nsim );

%Calculate employment and wage by age
%Now I need to make sure I just take the mean over accepted wages (not
%including zeros) because in some period we never observe a person works across
%simulations. In that case, in the data wage=0 for that guy at that period.
fracworktsim0=zeros(1,15);
meanwagetsim0=zeros(1,15);
for t=1:15
    worktsim0=simdata0(t:15:end,3);
    fracworktsim0(t)=mean(worktsim0);
    wagetsimaux0=simdata0(t:15:end,4);
    wagetsim0=wagetsimaux0(wagetsimaux0>0);
    meanwagetsim0(t)=mean(wagetsim0);
end

%%
%Counterfactual 1, new policy from the beginning
newexemp=2;
newtau=0.7;
[epsstarnew,epsstarstarnew ] = emax(parahat,data,newexemp,newtau);
nind=1000;
nsim=1000;
%Under the new threshold, simulate the economy.
[q1,works1,wages1]=simul2(parahat,epsstarnew,epsstarstarnew,nind,nsim,data);

%Construct simulated samples as before.
simdata1 =simdataset( data,q1,works1,wages1,nsim );

%Calculate employment and wage by age
%Now I need to make sure I just take the mean over accepted wages (not
%including zeros) because in some period we never observe a person works across
%simulations. In that case, in the data wage=0 for that guy at that period.
fracworktsim1=zeros(1,15);
meanwagetsim1=zeros(1,15);
for t=1:15
    worktsim1=simdata1(t:15:end,3);
    fracworktsim1(t)=mean(worktsim1);
    wagetsimaux1=simdata1(t:15:end,4);
    wagetsim1=wagetsimaux1(wagetsimaux1>0);
    meanwagetsim1(t)=mean(wagetsim1);
end


%%
%Counterfactual 2.
%I need to modify emax function to accomodate the change in tax rate.
exemp=3;
tau=0.5;
newexemp=2;
newtau=0.7;
[epsstarnew2,epsstarstarnew2 ] = emax2(parahat,data,exemp,tau,newexemp,newtau);
nind=1000;
nsim=100;
%Under the new threshold, simulate the economy.
[q2,works2,wages2]=simul2(parahat,epsstarnew2,epsstarstarnew2,nind,nsim,data);
simdata2 =simdataset( data,q2,works2,wages2,nsim );

fracworktsim2=zeros(1,15);
meanwagetsim2=zeros(1,15);
for t=1:15
    worktsim2=simdata2(t:15:end,3);
    fracworktsim2(t)=mean(worktsim2);
    wagetsimaux2=simdata2(t:15:end,4);
    wagetsim2=wagetsimaux2(wagetsimaux2>0);
    meanwagetsim2(t)=mean(wagetsim2);
end

%%
%Counterfactual 3
%Since they don't expect the policy change, their threshold comes from the
%original specification till period 7, while it comes from new
%specification from period 8.
[epsstartest,epsstarstartest ] = emax(parahat,data,3,0.5);
[epsstartest2,epsstarstartest2 ] = emax(parahat,data,2,0.7);
epsstarnew3=zeros(nind,15,15);
epsstarstarnew3=zeros(nind,15,15);
for t=1:7
    epsstarnew3(:,:,t)=epsstartest(:,:,t);
    epsstarstarnew3(:,:,t)=epsstarstartest(:,:,t);
end
for t=8:15
    epsstarnew3(:,:,t)=epsstartest2(:,:,t);
    epsstarstarnew3(:,:,t)=epsstarstartest2(:,:,t);
end
nind=1000;
nsim=100;
%Under the new threshold, simulate the economy.
[q3,works3,wages3]=simul2(parahat,epsstarnew3,epsstarstarnew3,nind,nsim,data);
simdata3 =simdataset( data,q3,works3,wages3,nsim );

fracworktsim3=zeros(1,15);
meanwagetsim3=zeros(1,15);
for t=1:15
    worktsim3=simdata3(t:15:end,3);
    fracworktsim3(t)=mean(worktsim3);
    wagetsimaux3=simdata3(t:15:end,4);
    wagetsim3=wagetsimaux3(wagetsimaux3>0);
    meanwagetsim3(t)=mean(wagetsim3);
end

%%
input3.data = [fracworktsim0;fracworktsim1;fracworktsim2;fracworktsim3]
input4.data=[[1:15];meanwagetsim0;meanwagetsim1;meanwagetsim2;meanwagetsim3 ]
input4.transposeTable = 1;
input4.tableBorders = 0;
latex1 = latexTable(input4)
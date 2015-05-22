clear all;
rng(278);

%Set model
%Consumption utility u(wm+wf)=(wm+wf)^(1-ucurve)/(1-ucurve). 
%Utility from child v(nt)=vfir*nt-vsec*nt^2.
%Childcare cost c_t(a_t)=-max(0,cint-cslope*at).
%parameters are para=[ucurve,vcurve,cint,cslope].

%%%%%%%%%%%%%%
%NO STOCHASTIC UTILITY INTRODUCED YET. ALL THE DECISIONS ARE DETERMINISTIC.
%%%%%%%%%%%%%%

%Flow utility stored in "flowu.m".

%Set hypothetical para
%flow utility para
ucurve=0.3;
vfir=0.5;
vsec=3;
cint=2.5;
cslope=0.5;
cage=0.01;

%dynamic para, borrowed from Dey Flinn
beta=0.9;
deltam=0.032;
deltaf=0.05;
lambdam=0.279;
lambdaf=0.243;
mum=2.514;
muf=2.403;
sigmam=0.139;
sigmaf=0.135;
%Child arrival probability
gamma=0.3;

%%
%Set state space

%Set demographic (outside income)
y=0.5;

%%%%%%%%%%
%Set terminal value to be zero (for now).
%%%%%%%%%%
term=0;

%State space:
%wm,wf: continuous, following log-normal. Discretize them into 10 points (for
%now)
wm=linspace(5,20,10);
wf=linspace(5,20,10);

%Number of children: up to 5. Child age up to 5 (for now, because at age 5
%c(at)=0.
nt=0:5;
at=1:5;

%%
%The number of state space is 10*10*6*5=3000. 
%Make it a one-dimensional vector.
[a,b,c,d]=ndgrid(wm,wf,nt,at);
state=[a(:),b(:),c(:),d(:)];
wmstate=state(:,1);
wfstate=state(:,2);
ntstate=state(:,3);
atstate=state(:,4);
%Order of state: wm first, then wf, then nt, and finally at changes.
sizestate=numel(state(:,1));


%For each state, they have up to four choices (either dft= 0 or 1 
%and dmt= 0 or 1). The number of choices depend on the number of offers
%and/or their working status in the previous period.

%Depending on the choices made in the previous period, transition matrix
%differs. Compute conditional transition matrix (3000*3000) for each of
%four choices.




%%
%For each state, they have four choices (either dft= 0 or 1 and dmt= 0 or 1).
%Compute utility corresponding to the four

%To use array fun, make "para" 3000*1 dimensional.
simucurve=ones(sizestate,1)*ucurve;
simvfir=ones(sizestate,1)*vfir;
simvsec=ones(sizestate,1)*vsec;
simcint=ones(sizestate,1)*cint;
simcslope=ones(sizestate,1)*cslope;
simcage=ones(sizestate,1)*cage;

yy=ones(sizestate,1)*y;
%Also, make working decision 3000*1 dimensional
work=ones(sizestate,1);
nowork=zeros(sizestate,1);
%%
%Simulate the model for 20 periods.

%Terminal period
t=20;
tt=ones(sizestate,1)*t;

%Conditional utilities corresponding to four cases.
%Man work, woman work
VeeT=arrayfun(@flowu,simucurve,simvfir,simvsec,simcint,simcslope,simcage,wmstate,wfstate,ntstate,atstate,work,work,yy,tt)+term;
%Man work, woman not
VeuT=arrayfun(@flowu,simucurve,simvfir,simvsec,simcint,simcslope,simcage,wmstate,wfstate,ntstate,atstate,work,nowork,yy,tt)+term;
%Man not, woman work
VueT=arrayfun(@flowu,simucurve,simvfir,simvsec,simcint,simcslope,simcage,wmstate,wfstate,ntstate,atstate,nowork,work,yy,tt)+term;
%Man not, woman not
VuuT=arrayfun(@flowu,simucurve,simvfir,simvsec,simcint,simcslope,simcage,wmstate,wfstate,ntstate,atstate,nowork,nowork,yy,tt)+term;

%%%%%%%%%%%%%%%%
%Might worth considering%
%Full state space only needed for ee. For uu, it only depends on nt and at.
%i.e. 30 dimensions are enough.
%Does this lose speed? I don't know. In calculating Emax later, I need to
%boost VuuT to 3000*1 by taking Kronecker product anyway.
%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%
%ITERATION OVER t TO BE INSERTED HERE.
%%%%%%%%%%%%%%


%Now we are in between Period 19 and Period 20.

%First solve "Period 20's optimal choice". i.e. value and policy function
%conditional on realized wage offer at period 20.

%Two choices: work and fertility. First consider working decision.

%WORKING DECISION
%This can be calculated by two-step. First given everything observed (i.e.
%period 20 state already realized), calculate optimal choice CONDITIONAL on
%whether each spouse can work or not.

%If two of them can choose whether to work (two offers, one already working 
%and one offer,etc)
%Pick the choice that maximizes value
aux=[VeeT,VeuT,VueT,VuuT];
[VmfT,policymfT]=max(aux,[],2);
%mf implies both of them can work.

%If only  male can choose (female not getting offer, job just
%destroyed etc...)
aux2=[VeuT,VuuT];
[VmT,policymT]=max(aux2,[],2);

%If only female can choose
aux3=[VeuT,VuuT];
[VfT,policyfT]=max(aux3,[],2);

%If noone can work (both unemployed previously and no offers, or one
%unemployed previously and the other job destroyed, etc)
VT=VuuT;
policyT=4*ones(sizestate,1);

clear aux aux2 aux3
%This yields policy function and conditional value functions (conditional
%on wages. 

%%

%To compute Emax (evaluated at period 19), take expectation over
%wage. 


%Assume log-normal wage.
distwm=lognpdf(wm,mum,sigmam);
distwf=lognpdf(wf,muf,sigmaf);


%Note that the variable we are taking expectation over differs depending on
%which state we come from.
%Consider equation (2) in the paper (see updated version 1.2 on May 21).
%This expectation is taken over wm. Not over wf because wf value is
%predetermined at previous period (note state ue). If we consider Vt^eu, then we need to
%take expectation over wf instead.

%So Emax is path-dependent. Need for computing them each by each.
%Start from equation number (1)-(6) in the paper.

%Equation 2 requires taking expectation of Vmf over wm
%Calculate 300*1 vector corresponding to expectation over wm.
aux=(distwm*reshape(VmfT,numel(wm),numel(wf)*numel(nt)*numel(at))).';
%Align the order of rows to match them (i.e. LHS state (wf,nt,at) =RHS
%state (wf,nt,at+1) for each row).
%Delete top rows with at=0.
aux(1:numel(wf)*numel(nt))=[];
%%%%%%%%%%%Treatment at corner%%%%%%%%%%%%%%%%%
%If a19=5, then set a20=5 as well (it doesn't matter anyway).
%Do the same for all cases below.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aux((1+numel(wf)*numel(nt)*(numel(at)-1)):numel(wf)*numel(nt)*numel(at))=...
    aux((1+numel(wf)*numel(nt)*(numel(at)-2)):numel(wf)*numel(nt)*(numel(at)-1));
%Finally, make it 3000*1 again.
eq2=kron(aux,ones(numel(wm)));
clear aux

%equation 3 requires expectation of Vm over wm.
aux=(distwm*reshape(VmT,numel(wm),numel(wf)*numel(nt)*numel(at))).';
%Same procedure as above.
aux(1:numel(wf)*numel(nt))=[];
aux((1+numel(wf)*numel(nt)*(numel(at)-1)):numel(wf)*numel(nt)*numel(at))=...
    aux((1+numel(wf)*numel(nt)*(numel(at)-2)):numel(wf)*numel(nt)*(numel(at)-1));

eq3=kron(aux,ones(numel(wm)));
clear aux

%Equation 4 just requires alignment of at.
aux=VT;
aux(1:numel(wm)*numel(wf)*numel(nt))=[];
aux((1+numel(wm)*numel(wf)*numel(nt)*(numel(at)-1)):numel(wm)*numel(wf)*numel(nt)*numel(at))=...
    aux((1+numel(wm)*numel(wf)*numel(nt)*(numel(at)-2)):numel(wm)*numel(wf)*numel(nt)*(numel(at)-1));

eq4=aux;
clear aux

%Same for equation 6
aux=VfT;
aux(1:numel(wm)*numel(wf)*numel(nt))=[];
aux((1+numel(wm)*numel(wf)*numel(nt)*(numel(at)-1)):numel(wm)*numel(wf)*numel(nt)*numel(at))=...
    aux((1+numel(wm)*numel(wf)*numel(nt)*(numel(at)-2)):numel(wm)*numel(wf)*numel(nt)*(numel(at)-1));

eq6=aux;
clear aux

%Now, for equation 5, I need to compare Vue(nt,at+1) and Vue(nt+1,1).
%Note that I need VueT. NO Need for VfT. This is because their job status hasn't
%changed from the previous period (essentially no choice today. Last
%period's optimal choice=today's optimal choice).

%If give a birth,
%add 1 to Nt, at goes back to 1.
%First, slash Nt=0 region and at>1 region. 
aux=VueT((1+numel(wm)*numel(wf)):numel(wm)*numel(wf)*numel(nt));
%%%%%%%%%%%%%%%Treatment at corner%%%%%%%%%%%%%%%%
%Once N=5, then by giving a birth N stays 5. Only at is re-set to 1
%This makes it suboptimal to give a birth when N=5 by construction.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This gives me Nt \in [1,2,3,4,5,5] and at=1 region.
aux((1+numel(wm)*numel(wf)*(numel(nt)-1)):numel(wm)*numel(wf)*numel(nt))=aux((1+numel(wm)*numel(wf)*(numel(nt)-2)):numel(wm)*numel(wf)*(numel(nt)-1));
%Repeat them numel(at) times. That gives me the appropriately updated state
%space.
aux2=repmat(aux,numel(at),1);

%If not give a birth, same as before.
aux3=VueT;
aux3(1:numel(wm)*numel(wf)*numel(nt))=[];
aux3((1+numel(wm)*numel(wf)*numel(nt)*(numel(at)-1)):numel(wm)*numel(wf)*numel(nt)*numel(at))=...
    aux3((1+numel(wm)*numel(wf)*numel(nt)*(numel(at)-2)):numel(wm)*numel(wf)*numel(nt)*(numel(at)-1));
aux4=[aux2,aux3];

[eq5,Pbue]=max(aux4,[],2);
clear aux2 aux3 aux4
%%




clear all;
rng(278);

% Consumption utility is u(wm + wf) = (wm + wf)^(1 - ucurve) / (1 - ucurve)
% Utility from child     v(nt)      = vfir * nt - vsec * nt^2
% Childcare cost         c_t(a_t)   = -max(0, cint - cslope * at)
% Parameters are         para       = [ucurve, vcurve, cint, cslope]

%%%%%%%%%%%%%%
%NO STOCHASTIC UTILITY INTRODUCED YET. ALL THE DECISIONS ARE DETERMINISTIC.
%%%%%%%%%%%%%%

% Flow utility is stored in "flowu.m"

% Set hypothetical para
% Flow utility para
ucurve = 0.3;
vfir   = 0.5;
vsec   = 3;
cint   = 2.5;
cslope = 0.5;
cage   = 0.01;

% Dynamic parameters, borrowed from Dey Flinn
Beta    = 0.9;  % Capitalized to prevent confusion with beta function
deltam  = 0.032;
deltaf  = 0.05;
lambdam = 0.279;
lambdaf = 0.243;
mum     = 2.514;
muf     = 2.403;
sigmam  = 0.139;
sigmaf  = 0.135;

%Child arrival probability
Gamma = 0.3;  % Capitalized to prevent confusion with the gamma function

%%
% Set state space

% Set demographic (outside income)
y = 0.5;

%%%%%%%%%%
% Set terminal value to be zero (for now).
%%%%%%%%%%
term = 0;

% State space:
% wm, wf: continuous, following log-normal. 
% Discretize them into 10 points (for now)
wm = linspace(5, 20, 10);
wf = linspace(5, 20, 10);

% Number of children: up to 5. Child age up to 5 (for now, because at age 5 
% c(at)=0.

% Update 5/22 (model v1.2 update)
% Slight modification from the paper made. Used to be: if nt=0, at=0 and no
% evolution.
% Now: no matter what nt is, at evolves in the same way, but if nt=0 then at
% does not show up in utility.
% This makes the evolution of state space simpler.
nt = 0:5;
at = 1:5;

%%
% The size of the state space is 10 * 10 * 6 * 5 = 3000. 
% Make it a one-dimensional vector.
[a, b, c, d] = ndgrid(wm, wf, nt, at);
state        = [a(:), b(:), c(:), d(:)];
wmstate      = state(:, 1);
wfstate      = state(:, 2);
ntstate      = state(:, 3);
atstate      = state(:, 4);
% Order of state: wm first, then wf, then nt, and finally at changes.
sizestate    = numel(state(:, 1));

clear a b c d



%%
% To use arrayfun later, make parameters 3000*1 dimensional object as well.
simucurve = ones(sizestate, 1) * ucurve;
simvfir   = ones(sizestate, 1) * vfir;
simvsec   = ones(sizestate, 1) * vsec;
simcint   = ones(sizestate, 1) * cint;
simcslope = ones(sizestate, 1) * cslope;
simcage   = ones(sizestate, 1) * cage;

yy = ones(sizestate, 1) * y;

% Also, make working decision 3000*1 dimensional
work   = ones(sizestate, 1);
nowork = zeros(sizestate, 1);

% Now all the inputs to the function "flowu" are set as 3000*1 vectors.

%%
% Simulate the model for 20 periods. 
period = 20;

% To start backward induction, we need to have utility corresponding to four
% choices (male/female work/not work) at the terminal period.

% Utility is nothing but just flow utility. So calculate flow utility at all
% the state points.

% Now we are in terminal period
t  = period;
tt = ones(sizestate,1)*t;

% Conditional utilities corresponding to four cases.
% Man work, woman work
VeeT = arrayfun(@flowu,simucurve,simvfir,simvsec,simcint,simcslope,simcage,...
                wmstate,wfstate,ntstate,atstate,work,work,yy,tt)+term;
            
% Man work, woman not
VeuT = arrayfun(@flowu,simucurve,simvfir,simvsec,simcint,simcslope,simcage,...
                wmstate,wfstate,ntstate,atstate,work,nowork,yy,tt)+term;
            
% Man not, woman work
VueT = arrayfun(@flowu,simucurve,simvfir,simvsec,simcint,simcslope,simcage,...
                wmstate,wfstate,ntstate,atstate,nowork,work,yy,tt)+term;
          
% Man not, woman not
VuuT = arrayfun(@flowu,simucurve,simvfir,simvsec,simcint,simcslope,simcage,...
                wmstate,wfstate,ntstate,atstate,nowork,nowork,yy,tt)+term;
            
% This is conditional on realized wage AND CHOICES. i.e. THIS IS NOT A VALUE
% FUNCTION AT PERIOD 20. 

%%%%%%%%%%%%%%%%
% Might worth considering%
% Full state space only needed for ee. For uu, it only depends on nt and at.
% i.e. 30 dimensions are enough.
% Does this lose speed? I don't know. In calculating Emax later, I need to
% boost VuuT to 3000*1 by taking Kronecker product anyway.
%%%%%%%%%%%%%%%%


%Construct matrices to store values and policies. The dimension is number of state space *
%number of periods
Veemat  = zeros(sizestate,period);
Veumat  = zeros(sizestate,period);
Vuemat  = zeros(sizestate,period);
Vuumat  = zeros(sizestate,period);
Vmfmat  = zeros(sizestate,period);
Vmmat   = zeros(sizestate,period);
Vfmat   = zeros(sizestate,period);
Vumat   = zeros(sizestate,period);
Pmfmat  = zeros(sizestate,period);
Pmmat   = zeros(sizestate,period);
Pfmat   = zeros(sizestate,period);
Pumat   = zeros(sizestate,period);
Pbeemat = zeros(sizestate,period);
Pbeumat = zeros(sizestate,period);
Pbuemat = zeros(sizestate,period);
Pbuumat = zeros(sizestate,period);

% Store terminal period value that I calculated above.
Veemat(:,t) = VeeT;
Veumat(:,t) = VeuT;
Vuemat(:,t) = VueT;
Vuumat(:,t) = VuuT;


%%
%%%%%%
% Instruction of iteration
%%%%%%

% Below, I iterate the economy from period 19 to period 1.
% Each step proceeds as follows. First, given conditional utilities at
% period t+1 (conditional on realized wage and choice), calculate conditional
% value function at period t+1. i.e. calculate optimal choice at period t+1
% given offer/job destruction.
% Second, using conditional value function obtained, calculate interim value function
% between period t and t+1 (interim=offer realizes but not the value of
% wage, called EMAX in the lecture). i.e. RHS of value function equations (1)-(24) in the paper.
% Finally, calculate conditional utility function at period t (LHS=RHS, conditional
% on realized wage at t and choice at t). This
% completes period t's iteration.


%Iteration from period 19 to priod 1.
%Set period 20's conditional utility function as initial value.
Veeup = VeeT;
Veuup = VeuT;
Vueup = VueT;
Vuuup = VuuT;

for tau = 1:(period-1)
t = period-tau;

%Use the updated utility functions at period t+1.
Vee = Veeup;
Veu = Veuup;
Vuu = Vuuup;
Vue = Vueup;


% Step 1: Solve period t+1's optimal choice.

% Two choices: work and fertility. First consider working decision.

% Working decision
% Calculate optimal choice conditional on whether each spouse can make a choice or not.

% If two of them can choose whether to work (two offers, one already working 
% and one offer,etc)
% Pick the choice that maximizes value
aux       = [Vee,Veu,Vue,Vuu];
[Vmf,Pmf] = max(aux,[],2);
% mf implies both of them can work.

% If only  male can choose (female not getting offer, job just
% destroyed etc...)
aux2    = [Veu,Vuu];
[Vm,Pm] = max(aux2,[],2);

% If only female can choose
aux3    = [Vue,Vuu];
[Vf,Pf] = max(aux3,[],2);

% If noone can work (both unemployed previously and no offers, or one
% unemployed previously and the other job destroyed, etc)
Vu = Vuu;
Pu = 4*ones(sizestate,1);

clear aux aux2 aux3
% This yields policy function and conditional value functions (conditional
% on wages). 


% Also, I need to solve for policy function for fertility. But I will solve
% it jointly with calculating Emax functions, which is computationally
% simpler.



%%
% Second step: Compute interim value functions (EMAX).
% Given the optimal choice (for work) at period t+1, calculate RHS of equation 
% (1)-(24).

% To compute Emax taking expectation over wage needed. 
% Assume log-normal wage.
distwm = lognpdf(wm,mum,sigmam);
distwf = lognpdf(wf,muf,sigmaf);


%Note that the variable we are taking expectation over differs depending on
%which state we come from. Emax is path-dependent. Need for computing them each by each.
%Start from equation number (1)-(6) in the paper.

%Equation 2 requires taking expectation of Vmf over wm
%Calculate 300*1 vector corresponding to expectation over wm.
aux=(distwm*reshape(Vmf,numel(wm),numel(wf)*numel(nt)*numel(at))).';
%Align the order of rows to match them (i.e. LHS state (wf,nt,at) =RHS
%state (wf,nt,at+1) for each row).
%Delete top rows with at=1.
aux(1:numel(wf)*numel(nt))=[];
%%%%%%%%%%%Treatment at corner%%%%%%%%%%%%%%%%%
%If at=5, then set a{t+1}=5 as well (it doesn't matter anyway).
%Do the same for all cases below.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aux((1+numel(wf)*numel(nt)*(numel(at)-1)):numel(wf)*numel(nt)*numel(at))=...
    aux((1+numel(wf)*numel(nt)*(numel(at)-2)):numel(wf)*numel(nt)*(numel(at)-1));
%Finally, make it 3000*1 again.
eq2=kron(aux,ones(numel(wm),1));
clear aux
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Equation 2 to be replaced by simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%equation 3 requires expectation of Vm over wm.
aux=(distwm*reshape(Vm,numel(wm),numel(wf)*numel(nt)*numel(at))).';
%Same procedure as above.
aux(1:numel(wf)*numel(nt))=[];
aux((1+numel(wf)*numel(nt)*(numel(at)-1)):numel(wf)*numel(nt)*numel(at))=...
    aux((1+numel(wf)*numel(nt)*(numel(at)-2)):numel(wf)*numel(nt)*(numel(at)-1));

eq3=kron(aux,ones(numel(wm),1));
clear aux

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Equation 3 to be replaced by simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Equation 4 just requires alignment of at.
aux=Vu;
aux(1:numel(wm)*numel(wf)*numel(nt))=[];
aux((1+numel(wm)*numel(wf)*numel(nt)*(numel(at)-1)):numel(wm)*numel(wf)*numel(nt)*numel(at))=...
    aux((1+numel(wm)*numel(wf)*numel(nt)*(numel(at)-2)):numel(wm)*numel(wf)*numel(nt)*(numel(at)-1));

eq4=aux;
clear aux

%Same for equation 6. No choices in this case, so use Vue instead of Vf.
aux=Vue;
aux(1:numel(wm)*numel(wf)*numel(nt))=[];
aux((1+numel(wm)*numel(wf)*numel(nt)*(numel(at)-1)):numel(wm)*numel(wf)*numel(nt)*numel(at))=...
    aux((1+numel(wm)*numel(wf)*numel(nt)*(numel(at)-2)):numel(wm)*numel(wf)*numel(nt)*(numel(at)-1));

eq6=aux;
clear aux

%Now, for equation 5, I need to two things\ solve optimization at t+1 concerning 
%childbirth (which I didn't do above). i.e. compare Vue(nt,at+1) and Vue(nt+1,1).
%Since their job status hasn't changed from the previous period (essentially
%no choice today. Last period's optimal choice=today's optimal choice), I
%only need Vue and not Vf.

%Then, I need to align the state space values to make equation 5 hold.
%Do these simultaneously.

%If give a birth,
%add 1 to Nt, at goes back to 1.
%First, slash Nt=0 region and at>1 region. 
aux=Vue((1+numel(wm)*numel(wf)):numel(wm)*numel(wf)*numel(nt));
%%%%%%%%%%%%%%%Treatment at corner%%%%%%%%%%%%%%%%
%Once N=5, then by giving a birth N stays 5. Only at is re-set to 1
%This makes it suboptimal to give a birth when N=5 by construction.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This gives me Nt \in [1,2,3,4,5,5] and at=1 region.
aux((1+numel(wm)*numel(wf)*(numel(nt)-1)):numel(wm)*numel(wf)*numel(nt))=aux((1+numel(wm)*numel(wf)*(numel(nt)-2)):numel(wm)*numel(wf)*(numel(nt)-1));
%Repeat them numel(at) times. That gives me the appropriately aligned state
%space values.
aux2=repmat(aux,numel(at),1);

%If not give a birth, same as before.
aux3=Vue;
aux3(1:numel(wm)*numel(wf)*numel(nt))=[];
aux3((1+numel(wm)*numel(wf)*numel(nt)*(numel(at)-1)):numel(wm)*numel(wf)*numel(nt)*numel(at))=...
    aux3((1+numel(wm)*numel(wf)*numel(nt)*(numel(at)-2)):numel(wm)*numel(wf)*numel(nt)*(numel(at)-1));
aux4=[aux2,aux3];

%Policy function denoted as Pbue.
[eq5,Pbue]=max(aux4,[],2);
%Pbue=1 give birth, =2 not give birth.

clear aux aux2 aux3 aux4


%Equation 8=Equation 6.
eq8=eq6;

%Equation 9
aux=Veu;
aux(1:numel(wm)*numel(wf)*numel(nt))=[];
aux((1+numel(wm)*numel(wf)*numel(nt)*(numel(at)-1)):numel(wm)*numel(wf)*numel(nt)*numel(at))=...
    aux((1+numel(wm)*numel(wf)*numel(nt)*(numel(at)-2)):numel(wm)*numel(wf)*numel(nt)*(numel(at)-1));
eq9=aux;
clear aux

%Equation 10=Equation 4.
eq10=eq4;

%Eguation 11=similar to equation 5.
%If not giving birth,
aux=Vee;
aux(1:numel(wm)*numel(wf)*numel(nt))=[];
aux((1+numel(wm)*numel(wf)*numel(nt)*(numel(at)-1)):numel(wm)*numel(wf)*numel(nt)*numel(at))=...
    aux((1+numel(wm)*numel(wf)*numel(nt)*(numel(at)-2)):numel(wm)*numel(wf)*numel(nt)*(numel(at)-1));
%If giving birth without quitting, same procedure as the first part of
%equation 5.
aux2=Vee((1+numel(wm)*numel(wf)):numel(wm)*numel(wf)*numel(nt));
aux2((1+numel(wm)*numel(wf)*(numel(nt)-1)):numel(wm)*numel(wf)*numel(nt))=aux2((1+numel(wm)*numel(wf)*(numel(nt)-2)):numel(wm)*numel(wf)*(numel(nt)-1));
%Repeat them numel(at) times. That gives me the appropriately updated state
%space.
aux3=repmat(aux2,numel(at),1);

%If give birth and male quit, same as above using Vue.
aux4=Vue((1+numel(wm)*numel(wf)):numel(wm)*numel(wf)*numel(nt));
aux4((1+numel(wm)*numel(wf)*(numel(nt)-1)):numel(wm)*numel(wf)*numel(nt))=aux4((1+numel(wm)*numel(wf)*(numel(nt)-2)):numel(wm)*numel(wf)*(numel(nt)-1));
%Repeat them numel(at) times. That gives me the appropriately updated state
%space.
aux5=repmat(aux4,numel(at),1);

%If give birth and female quit
aux6=Veu((1+numel(wm)*numel(wf)):numel(wm)*numel(wf)*numel(nt));
aux6((1+numel(wm)*numel(wf)*(numel(nt)-1)):numel(wm)*numel(wf)*numel(nt))=aux6((1+numel(wm)*numel(wf)*(numel(nt)-2)):numel(wm)*numel(wf)*(numel(nt)-1));
%Repeat them numel(at) times. That gives me the appropriately updated state
%space.
aux7=repmat(aux6,numel(at),1);
aux8=[aux3,aux,aux5,aux7];

%Policy function denoted as Pbue.
[eq11,Pbee]=max(aux8,[],2);
%Policy: 1=give birth no quit, 2 no give birth, 3=give birth male quit, 
%4=give birth female quit.

clear aux aux2 aux3 aux4 aux5 aux6 aux7 aux8


%Equation 12
aux=Vee;
aux(1:numel(wm)*numel(wf)*numel(nt))=[];
aux((1+numel(wm)*numel(wf)*numel(nt)*(numel(at)-1)):numel(wm)*numel(wf)*numel(nt)*numel(at))=...
    aux((1+numel(wm)*numel(wf)*numel(nt)*(numel(at)-2)):numel(wm)*numel(wf)*numel(nt)*(numel(at)-1));
eq12=aux;
clear aux


%Equation 14
%Take expectation of Vmf over both wm and wf.
%First, take expectation over wm. Exactly the same procedure as in equation
%2.
aux=(distwm*reshape(Vmf,numel(wm),numel(wf)*numel(nt)*numel(at))).';
aux(1:numel(wf)*numel(nt))=[];
aux((1+numel(wf)*numel(nt)*(numel(at)-1)):numel(wf)*numel(nt)*numel(at))=...
    aux((1+numel(wf)*numel(nt)*(numel(at)-2)):numel(wf)*numel(nt)*(numel(at)-1));
aux2=kron(aux,ones(numel(wm),1));

%Now, take expectation over wf. No need for adjustment of nt anymore. So
%just take expectation, and restore the size by Kronecker.
%Sort rows with respect to wf, so that I can use simple multiplicative form.
aux3=[state,aux2];
aux4=sortrows(aux3,[4,3,1]);
aux5=aux4(:,5);
%Do the same as above.
aux6=(distwf*reshape(aux5,numel(wf),numel(wm)*numel(nt)*numel(at))).';
aux7=kron(aux6,ones(numel(wf),1));
%Re-sort them to original order
aux8=[state,aux7];
aux9=sortrows(aux8,[4,3,2]);

eq14=aux9(:,5);
clear aux aux2 aux3 aux4 aux5 aux6 aux7 aux8 aux9

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Equation 14 to be replaced by simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Equation 15=Equation 3
eq15=eq3;

%Equation 16
%Sort rows with respect to wf.
aux=[state,Vf];
aux2=sortrows(aux,[4,3,1]);
aux3=aux2(:,5);
%Do the same thing as equation 3.
aux4=(distwf*reshape(aux3,numel(wf),numel(wm)*numel(nt)*numel(at))).';
aux4(1:numel(wm)*numel(nt))=[];
aux4((1+numel(wm)*numel(nt)*(numel(at)-1)):numel(wm)*numel(nt)*numel(at))=...
    aux4((1+numel(wm)*numel(nt)*(numel(at)-2)):numel(wm)*numel(nt)*(numel(at)-1));
aux5=kron(aux4,ones(numel(wf),1));
%Re-sort them to original order
aux6=[state,aux5];
aux7=sortrows(aux6,[4,3,2]);

eq16=aux7(:,5);
clear aux aux2 aux3 aux4 aux5 aux6 aux7

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Equation 16 to be replaced by simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Equation 17
%Similar to equation 5
%If give a birth,
%add 1 to Nt, at goes back to 1.
%First, slash Nt=0 region and at>1 region. 
aux=Vuu((1+numel(wm)*numel(wf)):numel(wm)*numel(wf)*numel(nt));
aux((1+numel(wm)*numel(wf)*(numel(nt)-1)):numel(wm)*numel(wf)*numel(nt))=aux((1+numel(wm)*numel(wf)*(numel(nt)-2)):numel(wm)*numel(wf)*(numel(nt)-1));
%Repeat them numel(at) times. That gives me the appropriately updated state
%space.
aux2=repmat(aux,numel(at),1);

%If not give a birth, same as before.
aux3=Vuu;
aux3(1:numel(wm)*numel(wf)*numel(nt))=[];
aux3((1+numel(wm)*numel(wf)*numel(nt)*(numel(at)-1)):numel(wm)*numel(wf)*numel(nt)*numel(at))=...
    aux3((1+numel(wm)*numel(wf)*numel(nt)*(numel(at)-2)):numel(wm)*numel(wf)*numel(nt)*(numel(at)-1));
aux4=[aux2,aux3];

%Policy function denoted as Pbuu.
[eq17,Pbuu]=max(aux4,[],2);

clear aux aux2 aux3 aux4


%Equation 18=equation 4
eq18=eq4;


%Case eu is the mirror image of case ue, with just re-ordering required.
%Equation 20 (see model version 1.3 updated on May 22)
%Sort rows with respect to wf.
aux=[state,Vmf];
aux2=sortrows(aux,[4,3,1]);
aux3=aux2(:,5);
%Do the same thing as equation 2.
aux4=(distwf*reshape(aux3,numel(wf),numel(wm)*numel(nt)*numel(at))).';
%Align the order of rows to match them (i.e. LHS state (wf,nt,at) =RHS
%state (wf,nt,at+1) for each row).
%Delete top rows with at=0.
aux4(1:numel(wm)*numel(nt))=[];
aux4((1+numel(wm)*numel(nt)*(numel(at)-1)):numel(wm)*numel(nt)*numel(at))=...
    aux4((1+numel(wm)*numel(nt)*(numel(at)-2)):numel(wm)*numel(nt)*(numel(at)-1));
aux5=kron(aux4,ones(numel(wf),1));
%Re-sort them to original order
aux6=[state,aux5];
aux7=sortrows(aux6,[4,3,2]);


eq20=aux7(:,5);
clear aux aux2 aux3 aux4 aux5 aux6 aux7

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Equation 20 to be replaced by simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%equation 21 =Equation 16
eq21=eq16;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Equation 21 to be replaced by simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Equation 22 =Equation 4
eq22=eq4;

%Equation 23 similar to equation 5.
%If give a birth,
%add 1 to Nt, at goes back to 1.
%First, slash Nt=0 region and at>1 region. 
aux=Veu((1+numel(wm)*numel(wf)):numel(wm)*numel(wf)*numel(nt));
%This gives me Nt \in [1,2,3,4,5,5] and at=1 region.
aux((1+numel(wm)*numel(wf)*(numel(nt)-1)):numel(wm)*numel(wf)*numel(nt))=aux((1+numel(wm)*numel(wf)*(numel(nt)-2)):numel(wm)*numel(wf)*(numel(nt)-1));
%Repeat them numel(at) times. That gives me the appropriately updated state
%space.
aux2=repmat(aux,numel(at),1);

%If not give a birth, same as before.
aux3=Veu;
aux3(1:numel(wm)*numel(wf)*numel(nt))=[];
aux3((1+numel(wm)*numel(wf)*numel(nt)*(numel(at)-1)):numel(wm)*numel(wf)*numel(nt)*numel(at))=...
    aux3((1+numel(wm)*numel(wf)*numel(nt)*(numel(at)-2)):numel(wm)*numel(wf)*numel(nt)*(numel(at)-1));
aux4=[aux2,aux3];

%Policy function denoted as Pbeu.
[eq23,Pbeu]=max(aux4,[],2);

clear aux aux2 aux3 aux4

%Equation 24 =Equation 9.
eq24=eq9;

%All equations on RHS derived.
%Interim value function (EMAX) derivation done.

%Step 3
%Conditional utility at period t (conditional on choice and wage).
%Just calculate the four equations on the paper.
tt=ones(sizestate,1)*t;

%State ue
Vueup=arrayfun(@flowu,simucurve,simvfir,simvsec,simcint,simcslope,simcage,wmstate,wfstate,ntstate,atstate,nowork,work,yy,tt)+Beta*lambdam*(1-deltaf)*eq2+Beta*lambdam*deltaf*eq3+Beta*(1-lambdam)*deltaf*eq4+Beta*(1-lambdam)*Gamma*(1-deltaf)*eq5+Beta*(1-lambdam)*(1-Gamma)*(1-deltaf)*eq6;
%State ee
Veeup=arrayfun(@flowu,simucurve,simvfir,simvsec,simcint,simcslope,simcage,wmstate,wfstate,ntstate,atstate,work,work,yy,tt)+Beta*deltam*(1-deltaf)*eq8+Beta*(1-deltam)*deltaf*eq9+Beta*deltam*deltaf*eq10+Beta*Gamma*(1-deltam)*(1-deltaf)*eq11+Beta*(1-Gamma)*(1-deltam)*(1-deltaf)*eq12;
%State uu
Vuuup=arrayfun(@flowu,simucurve,simvfir,simvsec,simcint,simcslope,simcage,wmstate,wfstate,ntstate,atstate,nowork,nowork,yy,tt)+Beta*lambdam*lambdaf*eq14+Beta*lambdam*(1-lambdaf)*eq15+Beta*(1-lambdam)*lambdaf*eq16+Beta*(1-lambdam)*(1-lambdaf)*Gamma*eq17+Beta*(1-lambdam)*(1-lambdaf)*(1-Gamma)*eq18;
%State eu
Veuup=arrayfun(@flowu,simucurve,simvfir,simvsec,simcint,simcslope,simcage,wmstate,wfstate,ntstate,atstate,work,nowork,yy,tt)+Beta*lambdaf*(1-deltam)*eq20+Beta*lambdaf*deltam*eq21+Beta*(1-lambdaf)*deltam*eq22+Beta*(1-lambdaf)*(1-deltam)*Gamma*eq23+Beta*(1-lambdaf)*(1-deltam)*(1-Gamma)*eq24;


%Everything is done for period t. Store them in the matrices created
%before. Note that calculated policy is at period t+1. So they are
%allocated t+1 column. For interim value function (EMAX), I put it in t+1
%column as well (doesn't matter).
Veemat(:,t)=Veeup;
Veumat(:,t)=Veuup;
Vuemat(:,t)=Vueup;
Vuumat(:,t)=Vuuup;
Vmfmat(:,t+1)=Vmf;
Vmmat(:,t+1)=Vm;
Vfmat(:,t+1)=Vf;
Vumat(:,t+1)=Vu;
Pmfmat(:,t+1)=Pmf;
Pmmat(:,t+1)=Pm;
Pfmat(:,t+1)=Pf;
Pumat(:,t+1)=Pu;
Pbeemat(:,t+1)=Pbee;
Pbeumat(:,t+1)=Pbeu;
Pbuemat(:,t+1)=Pbue;
Pbuumat(:,t+1)=Pbuu;

%4 conditional utilities and 4 EMAX.
%4 Policy functions concerning working choice (conditional on offer arrival
%and wage) and 4 policies concerning childbirth (one each for four working
%status)
end


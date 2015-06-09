function [ smmobjective ] = smmobjective(para,parafixed)

rng(278);

%Parameters
ucurve = para(1);
vfir   = para(2);
vsec   = para(3);
cint   = para(4);
cslope = para(5);    
cageint   = para(6);
cageslope = para(7);
conskid= para(8);
mum     = para(9);
muf     = para(10);
sigmam  = para(11);
sigmaf  = para(12);
Gamma   = para(13);
deltam  = parafixed(1);
deltaf  = parafixed(2);
lambdam = parafixed(3);
lambdaf = parafixed(4);
Beta    = parafixed(5);


% Set state space

% Set demographic (outside income)
y = 0.5;  %ADJUST INDIVIDUALLY


% State space:
% wm, wf: continuous, following log-normal. 
botm=5;
topm=30;
botf=3;
topf=25;
wm = linspace(botm, topm, 70);
wf = linspace(botf, topf, 70);

% Number of children: up to 5. Child age up to 5 (for now, because at age 5 
% c(at)=0.

% Update 5/22 (model v1.2 update)
% Slight modification from the paper made. Used to be: if nt=0, at=0 and no
% evolution.
% Now: no matter what nt is, at evolves in the same way, but if nt=0 then at
% does not show up in utility.
% This makes the evolution of state space simpler.
nt = 0:5;  %ADJUST INDIVIDUALLY
at = 1:5;  %ADJUST INDIVIDUALLY

%DP periods
period = 20; %ADJUST INDIVIDUALLY

%%%%%%%
% In conducting SMM, adjust state space (y,bot,top,wm,wf,nt,at) size etc on
% this function.
%%%%%%%Preparation done



% Make it a one-dimensional vector.
[aux1, aux2, aux3, aux4] = ndgrid(wm, wf, nt, at);
state        = [aux1(:), aux2(:), aux3(:), aux4(:)];
wmstate      = state(:, 1);
wfstate      = state(:, 2);
ntstate      = state(:, 3);
atstate      = state(:, 4);
% Order of state: wm first, then wf, then nt, and finally at changes.
sizestate    = numel(state(:, 1));

clear aux1 aux2 aux3 aux4


% To use arrayfun later, make parameters 3000*1 dimensional object as well.
simucurve = ones(sizestate, 1) * ucurve;
simvfir   = ones(sizestate, 1) * vfir;
simvsec   = ones(sizestate, 1) * vsec;
simcint   = ones(sizestate, 1) * cint;
simcslope = ones(sizestate, 1) * cslope;
simcageint   = ones(sizestate, 1) * cageint;
simcageslope   = ones(sizestate, 1) * cageslope;
simconskid   = ones(sizestate, 1) * conskid;
simBeta   = ones(sizestate, 1) * Beta;

yy = ones(sizestate, 1) * y;

% Also, make working decision 3000*1 dimensional
work   = ones(sizestate, 1);
nowork = zeros(sizestate, 1);

% Now all the inputs to the function "flowu" are set as 3000*1 vectors.


%%%%%%%%%%
% From here, I denote "choice-based value functions" as Vee,Veu,Vue,Vuu.
% These are functions of state space AND CHOICES.
% I denote "value functions" as Vmf, Vm, Vf, Vu. These are functions only
% of state space.

% To start backward induction, we need to have utility corresponding to four
% choices (male/female work/not work) at the terminal period.

% Utility is nothing but just flow utility. So calculate flow utility at all
% the state points.

% Now we are in terminal period
t  = period;
tt = ones(sizestate,1)*t;

%Terminal period value as a continuation payoff terminal state
termee =  arrayfun(@term, simucurve, simvfir, simvsec, simcint,...
              simcslope, simcageint, simcageslope, simconskid, simBeta, wmstate, wfstate, ntstate,...
              atstate, work, work, yy, tt);
termeu = arrayfun(@term, simucurve, simvfir, simvsec, simcint,...
             simcslope, simcageint, simcageslope, simconskid, simBeta, wmstate, wfstate, ntstate,...
              atstate, work, nowork, yy, tt);
termue = arrayfun(@term, simucurve, simvfir, simvsec, simcint,...
            simcslope, simcageint, simcageslope, simconskid, simBeta, wmstate, wfstate, ntstate,...
            atstate, nowork, work, yy, tt);
termuu = arrayfun(@term, simucurve, simvfir, simvsec, simcint,...
             simcslope, simcageint, simcageslope, simconskid, simBeta, wmstate, wfstate, ntstate,...
              atstate, nowork, nowork, yy, tt);


% Choice-based value functions corresponding to four cases.
% Man work, woman work
VeeT = arrayfun(@flowu_2, simucurve, simvfir, simvsec, simcint,...
                simcslope, simcageint, simcageslope, simconskid, wmstate, wfstate, ntstate,...
                atstate, work, work, yy, tt) + termee;
            
% Man work, woman not
VeuT = arrayfun(@flowu_2, simucurve, simvfir, simvsec, simcint,...
                simcslope, simcageint, simcageslope, simconskid, wmstate, wfstate, ntstate,...
                atstate, work, nowork, yy, tt) + termeu;      
            
% Man not, woman work
VueT = arrayfun(@flowu_2, simucurve, simvfir, simvsec, simcint,...
                simcslope, simcageint, simcageslope, simconskid, wmstate, wfstate, ntstate,...
                atstate, nowork, work, yy, tt) + termue;
            
% Man not, woman not
VuuT = arrayfun(@flowu_2, simucurve, simvfir, simvsec, simcint,...
                simcslope, simcageint, simcageslope, simconskid, wmstate, wfstate, ntstate,...
                atstate, nowork, nowork, yy, tt) + termuu;
            

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
Pwbeemat = zeros(sizestate,period);
Pwbeumat = zeros(sizestate,period);
Pwbuemat = zeros(sizestate,period);
Pwbuumat = zeros(sizestate,period);

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
% Each step proceeds as follows. First, given choice based value functions at
% period t+1 (conditional on realized wage and choice), calculate 
% value functions at period t+1. i.e. calculate optimal choice at period t+1
% given offer/job destruction.
% Second, using value function at t+1 obtained, calculate  EMAX by taking
% appropriate integrals. i.e. RHS of equations (1)-(24) in the paper.
% Finally, calculate choice-based value functions at period t (LHS=RHS, conditional
% on realized wage at t and choice at t). This
% completes period t's iteration.


%Iteration from period 19 to priod 1.
%Set period 20's choice based value functions as initial value.
Veeup = VeeT;
Veuup = VeuT;
Vueup = VueT;
Vuuup = VuuT;

for tau = 1:(period-1)
t = period-tau;

%Use the updated choice based value functions at period t+1.
Vee = Veeup;
Veu = Veuup;
Vuu = Vuuup;
Vue = Vueup;


% Step 1: Solve period t+1's optimal choice. i.e. Derive value function at t+1.

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
[Vm,aux4] = max(aux2,[],2);
Pm = 2*aux4; % To make the indeces of choices aligned

% If only female can choose
aux3    = [Vue,Vuu];
[Vf,aux5] = max(aux3,[],2);
Pf = 2+aux5; % To make the indeces of choices aligned

% If noone can work (both unemployed previously and no offers, or one
% unemployed previously and the other job destroyed, etc)
Vu = Vuu;
Pu = 4*ones(sizestate,1); % To make the indeces of choices aligned

clear aux aux2 aux3 aux4 aux5
% This yields policy function and value functions (conditional only on
% state values and NOT on choices). 


% Also, I need to solve for policy function for fertility. But I will solve
% it jointly with calculating Emax functions, which is computationally
% simpler.




% Second step: Compute EMAX.
% i.e. Given the optimal choice (for work) at period t+1, calculate RHS of equation 
% (1)-(24).

% To compute Emax taking expectation over wage needed. 
% Assume log-normal wage.
distwm = ((top-bot)/numel(wm))*lognpdf(wm,mum,sigmam);
distwf = ((top-bot)/numel(wf))*lognpdf(wf,muf,sigmaf);


% Note that the variable we are taking expectation over differs depending on
% which state we come from. Emax is path-dependent. Need for computing them each by each.
% Start from equation number (1)-(6) in the paper.

% Equation 2 requires taking expectation of Vmf over wm
% Calculate 300*1 vector corresponding to expectation over wm.
aux = (distwm * reshape(Vmf, numel(wm), numel(wf)*numel(nt)*numel(at))).';
% Align the order of rows to match them (i.e. LHS state (wf,nt,at) =RHS
% state (wf,nt,at+1) for each row).
% Delete top rows with at=1.
aux(1 : numel(wf)*numel(nt)) = [];
%%%%%%%%%%%Treatment at corner%%%%%%%%%%%%%%%%%
% If at=5, then set a{t+1}=5 as well (it doesn't matter anyway).
% Do the same for all cases below.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aux((1+numel(wf)*numel(nt)*(numel(at)-1)) : numel(wf)*numel(nt)*numel(at))=...
    aux((1+numel(wf)*numel(nt)*(numel(at)-2)) : numel(wf)*numel(nt)*(numel(at)-1));
%Finally, make it 3000*1 again.
eq2 = kron(aux, ones(numel(wm),1));
clear aux
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Equation 2 to be replaced by simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%equation 3 requires expectation of Vm over wm.
aux = (distwm * reshape(Vm,numel(wm), numel(wf)*numel(nt)*numel(at))).';
%Same procedure as above.
aux(1 : numel(wf)*numel(nt)) = [];
aux((1+numel(wf)*numel(nt)*(numel(at)-1)) : numel(wf)*numel(nt)*numel(at))=...
    aux((1+numel(wf)*numel(nt)*(numel(at)-2)) : numel(wf)*numel(nt)*(numel(at)-1));

eq3 = kron(aux, ones(numel(wm),1));
clear aux

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Equation 3 to be replaced by simulation
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
aux4=[aux3,aux2];

%Value and policy
[eq5,Pbue]=max(aux4,[],2);
Pwbue=3*ones(sizestate,1);

%Pbue=1 not give birth, =2 give birth. I separate birth choice (stored in
%Pbue) from corresponding work choice (stored in Pwbue). In this case,
%there's no work choice involved, so Pwbue contains "only female keeps
%working".

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
aux8=[aux,aux3,aux5,aux7];

%Policy function denoted as Pbee.
[eq11,aux9]=max(aux8,[],2);
%Policy: 1 =no give birth, 2=give birth no quit, 3=give birth male quit, 
%4=give birth female quit.

% Split this policy into two: work and fertility
Pbee = ones(sizestate,1);
Pbee(aux9>1) = 2;
Pwbee = zeros(sizestate,1);
Pwbee(aux9<3) = 1;
Pwbee(aux9==3) = 3;
Pwbee(aux9==4) = 2;
% If Pbee =2, give birth, =1 no give birth.
% If Pwbee=1, both work, =2 only male work, =3 only female work.
% In this state, it makes a difference to separate the two choices.

clear aux aux2 aux3 aux4 aux5 aux6 aux7 aux8 aux9 aux10


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
aux4=[aux3, aux2];

%Policy function denoted as Pbuu.
[eq17,Pbuu]=max(aux4,[],2);
Pwbuu=4*ones(sizestate,1);

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
aux4=[aux3,aux2];

%Policy function denoted as Pbeu.
[eq23,Pbeu]=max(aux4,[],2);
Pwbeu = 2*ones(sizestate,1);

clear aux aux2 aux3 aux4

%Equation 24 =Equation 9.
eq24=eq9;

%All equations on RHS derived.
%EMAX derivation done.

%Step 3
%Choice-based value function at period t (conditional on choice and state at t).
%Just calculate the four equations on the paper.
tt=ones(sizestate,1)*t;

%State ue
Vueup=arrayfun(@flowu_2,simucurve,simvfir,simvsec,simcint,simcslope, simcageint, simcageslope,simconskid,wmstate,wfstate,ntstate,atstate,nowork,work,yy,tt)+Beta*lambdam*(1-deltaf)*eq2+Beta*lambdam*deltaf*eq3+Beta*(1-lambdam)*deltaf*eq4+Beta*(1-lambdam)*Gamma*(1-deltaf)*eq5+Beta*(1-lambdam)*(1-Gamma)*(1-deltaf)*eq6;
%State ee
Veeup=arrayfun(@flowu_2,simucurve,simvfir,simvsec,simcint,simcslope, simcageint, simcageslope,simconskid,wmstate,wfstate,ntstate,atstate,work,work,yy,tt)+Beta*deltam*(1-deltaf)*eq8+Beta*(1-deltam)*deltaf*eq9+Beta*deltam*deltaf*eq10+Beta*Gamma*(1-deltam)*(1-deltaf)*eq11+Beta*(1-Gamma)*(1-deltam)*(1-deltaf)*eq12;
%State uu
Vuuup=arrayfun(@flowu_2,simucurve,simvfir,simvsec,simcint,simcslope, simcageint, simcageslope,simconskid,wmstate,wfstate,ntstate,atstate,nowork,nowork,yy,tt)+Beta*lambdam*lambdaf*eq14+Beta*lambdam*(1-lambdaf)*eq15+Beta*(1-lambdam)*lambdaf*eq16+Beta*(1-lambdam)*(1-lambdaf)*Gamma*eq17+Beta*(1-lambdam)*(1-lambdaf)*(1-Gamma)*eq18;
%State eu
Veuup=arrayfun(@flowu_2,simucurve,simvfir,simvsec,simcint,simcslope, simcageint, simcageslope,simconskid,wmstate,wfstate,ntstate,atstate,work,nowork,yy,tt)+Beta*lambdaf*(1-deltam)*eq20+Beta*lambdaf*deltam*eq21+Beta*(1-lambdaf)*deltam*eq22+Beta*(1-lambdaf)*(1-deltam)*Gamma*eq23+Beta*(1-lambdaf)*(1-deltam)*(1-Gamma)*eq24;


%Everything is done for period t. Store them in the matrices created
%before. Note that calculated policies are at period t+1. So they are
%allocated to t+1 column. For value functions (not choice-based), I put it in t+1
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
Pwbeemat(:,t+1)=Pwbee;
Pwbeumat(:,t+1)=Pwbeu;
Pwbuemat(:,t+1)=Pwbue;
Pwbuumat(:,t+1)=Pwbuu;

% 4 value functions (mf,m,f,u) and 4 choice based value functions (ee,eu,ue,uu).
% 4 Policy functions concerning working choice (Pmf, Pm, Pf, Pu), 4 policy
% functions concerning working choice when child arrives (Pwbee, Pwbeu, Pwbue,
% Pwbuu) and 4 policies concerning childbirth (Pbee, Pbeu, Pbue, Pbuu).
% 
% Pmfmat 1 = both work, 2 = male work, 3 = female work, 4 = none work
% Pfmat 3 = female work, 4 = none work
% Pmmat 2 = male work, 4= none work
% Pumat 4 = none work
% Pwbeemat 1= both work, 2 = male work, 3 = female work,
% Pwbeumat 2= male work
% Pwbuemat 3= female work
% Pwbuumat 4= none work

% In other words, 1= both work, 2 = male work, 3= female work and 4 =none
% work no matter what the states are. This is why I aligned the number
% above (e.g. lines 228,233,238).

% Pbeemat, Pbuemat, Pbeumat, Pbuumat 1 = no give birth, 2 = give birth.

%%%%%%%%%%%%%%%%%%%%%%%%
% Note; by construction working choice at period t depends on period t's realized 
% state value, while childbirth choice at period t depend on period t-1's
% state value.
%%%%%%%%%%%%%%%%%%%%%%%%
end


%%
% I noticed that I need Policy at period 1, so let me calculate it here in
% addition.
aux       = [Veemat(:,1),Veumat(:,1),Vuemat(:,1),Vuumat(:,1)];
[~,Pmfmat(:,1)] = max(aux,[],2);

aux2      = [Veumat(:,1),Vuumat(:,1)];
[~,aux4] = max(aux2,[],2);
Pmmat(:,1) = 2*aux4;

aux3      = [Vuemat(:,1),Vuumat(:,1)];
[~, aux5] = max(aux3,[],2);
Pfmat(:,1) = 2+aux5;

%Vumat(:,1) = Vuumat(:,1);
Pumat(:,1) = 4*ones(sizestate,1);

clear aux aux2 aux3 aux4 aux5



%%%%%

%DP solved. Calculate moments by simulation.

%%%%%

% Simulation of 10000 individuals
nsim = 5000;

% Draw shocks: job arrival, wage offer, job destruction (for both spouses)
% and child arrival per period. In total, number of period * 7 shocks.

% Job arrival
simofferm = reshape(binornd(1,lambdam,period*nsim,1),nsim,period);
simofferf = reshape(binornd(1,lambdaf,period*nsim,1),nsim,period);

% Wage offer is drawn from log-normal.
% Generate random indeces 
aux = rand(nsim,period);
aux2 = rand(nsim,period);
aux3 = zeros(numel(wm),1);
aux4 = zeros(numel(wf),1);
for k=1:numel(wm)
aux3(k)=sum(distwm(1:k));
aux4(k)=sum(distwf(1:k));
end

simwofferm=zeros(nsim,period);
simwofferf=zeros(nsim,period);

for i=1:nsim
    for t=1:period
    simwofferm(i,t) = find(aux3>aux(i,t),1,'first');
    simwofferf(i,t) = find(aux4>aux2(i,t),1,'first');
    end
end


clear aux aux2 aux3 aux4

%Job destruction
simdestm = reshape(binornd(1,deltam,period*nsim,1),nsim,period);
simdestf = reshape(binornd(1,deltaf,period*nsim,1),nsim,period);


% Child arrival
simchildarr = reshape(binornd(1,Gamma,period*nsim,1),nsim,period);


%Create matrix to store behaviors and realized states
simwork = zeros(nsim,period);
simoffer= zeros(nsim,period);
simnkids = zeros(nsim,period);
simat    = ones(nsim,period);
simwm    = zeros(nsim,period);
simwf    = zeros(nsim,period);
simreswagem =zeros(nsim,period);
simreswagef =zeros(nsim,period);

%Just to align notation as before
simstate = zeros(nsim,period);


%Things used in the next section
quitbirthm=zeros(nsim,period);
quitbirthf=zeros(nsim,period);
comebackf=zeros(nsim,period);
comebackm=zeros(nsim,period);
quitdurm=zeros(nsim,8);
quitdurf=zeros(nsim,8);
quitagem=zeros(nsim,8);
quitagef=zeros(nsim,8);
wagedifm=zeros(nsim,16);
wagediff=zeros(nsim,16);

%For each iteration, seven variables to update.+things used in the next
%section.

% Calculate initial period
% Initial arrival rate set higher so that we have more employed young
% workers
initofferm = reshape(binornd(1,0.8,nsim,1),nsim,1);
initofferf = reshape(binornd(1,0.7,nsim,1),nsim,1);

% Set initial working status and state
for i = 1:nsim
    if initofferm(i) + initofferf(i)==2 
        simoffer(i,1) = 1;%Two offers
        simwm(i,1) = simwofferm(i,1);
        simwf(i,1) = simwofferf(i,1);
        simstate(i,1) = simwm(i,1)+(simwf(i,1)-1)*numel(wm);
           %Assign the location of state (combination of drawn wm and wf)      
        simwork(i,1) = Pmfmat(simstate(i,1),1);
        
    elseif initofferm(i)==1             
        simoffer(i,1) = 2;%Only male offer
        simwm(i,1) = simwofferm(i,1);
        simwf(i,1) = simwofferf(i,1);  % Irrelevant. Assign random value
        simstate(i,1) = simwm(i,1)+(simwf(i,1)-1)*numel(wm);
        simwork(i,1) = Pmmat(simstate(i,1),1);
        
    elseif initofferf(i)==1 
        simoffer(i,1) = 3;   %Only female offer
        simwm(i,1) = simwofferm(i,1);  % Irrelevant. Assign random value
        simwf(i,1) = simwofferf(i,1);  
        simstate(i,1) = simwm(i,1)+(simwf(i,1)-1)*numel(wm);
        simoffer(i,1) = 3;
        simwork(i,1) = Pfmat(simstate(i,1),1);
        
    else                                    
        simoffer(i,1) = 4;                  %No offer
        simwm(i,1) = simwofferm(i,1);  % Irrelevant. Assign random value
        simwf(i,1) = simwofferf(i,1);  % Irrelevant. Assign random value
        simstate(i,1) = simwm(i,1)+(simwf(i,1)-1)*numel(wm);
        simwork(i,1) = Pumat(simstate(i,1),1);
        
    end
end
%%
%Now solve for the sequence of choices
for i = 1:nsim
    for t = 2:period
        
        if simwork(i,t-1)==1 % If both worked in the previous period
            
            if simdestm(i,t)+simdestf(i,t)==2 % If both jobs destroyed
                simoffer(i,t) = 4;            %No offers
                simwm(i,t) = simwofferm(i,t);  % Irrelevant. Assign random value
                simwf(i,t) = simwofferf(i,t);  % Irrelevant. Assign random value   
                simnkids(i,t) =simnkids(i,t-1); %Kid doesn't change
                simat(i,t) = min(simat(i,t-1)+1,numel(at));    %Age added by 1
                
                simstate(i,t) = (simat(i,t)-1)*numel(nt)*numel(wf)*numel(wm)+...
                  (simnkids(i,t))*numel(wf)*numel(wm)+(simwf(i,t)-1)*numel(wm)+simwm(i,t);
                   %Assign state value
                simwork(i,t) = Pumat(simstate(i,t),t); %Working choice
                simreswagem(i,t) = reswagemuu(simstate(i,t),t);
                simreswagef(i,t) = reswagefuu(simstate(i,t),t);
            
            elseif simdestm(i,t)==1 %Only male job destroyed
                simoffer(i,t) = 3;             % Only female offer
                simwm(i,t) = simwofferm(i,t);  % Irrelevant. Assign random value
                simwf(i,t) = simwf(i,t-1);     % She is offered the same wage as before
                simnkids(i,t) =simnkids(i,t-1); %Kid doesn't change
                simat(i,t) = min(simat(i,t-1)+1,numel(at));    %Age added by 1
                
                simstate(i,t) = (simat(i,t)-1)*numel(nt)*numel(wf)*numel(wm)+...
                  (simnkids(i,t))*numel(wf)*numel(wm)+(simwf(i,t)-1)*numel(wm)+simwm(i,t);
                   
                simwork(i,t) = Pfmat(simstate(i,t),t);
                simreswagem(i,t) = reswagemue(simstate(i,t),t);
                simreswagef(i,t) = simwf(i,t);
                
            elseif simdestf(i,t)==1 % Only female job destroyed
                simoffer(i,t) = 2;          %Only male offer
                simwm(i,t) = simwm(i,t-1);  % He is offered the same wage
                simwf(i,t) = simwofferf(i,t);  % Irrelevant. Assign random value   
                simnkids(i,t) =simnkids(i,t-1);
                simat(i,t) = min(simat(i,t-1)+1,numel(at));
                
                simstate(i,t) = (simat(i,t)-1)*numel(nt)*numel(wf)*numel(wm)+...
                  (simnkids(i,t))*numel(wf)*numel(wm)+(simwf(i,t)-1)*numel(wm)+simwm(i,t);
                   
                simwork(i,t) = Pmmat(simstate(i,t),t);
                simreswagem(i,t) = simwm(i,t);
                simreswagef(i,t) = reswagefeu(simstate(i,t),t);
                
            elseif simchildarr(i,t)==1 %Child arrives
                simoffer(i,t) = 1;         %Offer status remains the same.
                simwm(i,t) = simwm(i,t-1);  % He is offered the same wage
                simwf(i,t) = simwf(i,t-1);  % She also.
                simnkids(i,t) = min(simnkids(i,t-1)+Pbeemat(simstate(i,t-1),t)-1,numel(nt));
                     %If give birth, add 1 to "simnkids".
                     
                if Pbeemat(simstate(i,t-1),t)== 1
                    % State evolution depends on choices. If not giving a birth
                    simat(i,t) = min(simat(i,t-1)+1,numel(at));
                    simstate(i,t) = (simat(i,t)-1)*numel(nt)*numel(wf)*numel(wm)+...
                  (simnkids(i,t))*numel(wf)*numel(wm)+(simwf(i,t)-1)*numel(wm)+simwm(i,t);
                       
                    simwork(i,t) = simwork(i,t-1); %Work status not changing
                    simreswagem(i,t) = simwm(i,t);
                    simreswagef(i,t) = simwf(i,t);

                elseif Pbeemat(simstate(i,t-1),t)== 2
                    % If giving a birth,
                 simat(i,t) = 1; %at re-set to 1.
                simstate(i,t) = (simat(i,t)-1)*numel(nt)*numel(wf)*numel(wm)+...
                  (simnkids(i,t))*numel(wf)*numel(wm)+(simwf(i,t)-1)*numel(wm)+simwm(i,t);
                
                 simwork(i,t)=Pwbeemat(simstate(i,t-1),t); %Work choice involved.
                      %Choice based on state at t-1. This is how Pwbee is
                      %defined.
                  simreswagem(i,t) = reswagemue(simstate(i,t),t);
                  simreswagef(i,t) = reswagefeu(simstate(i,t),t);
              
                end
                
            elseif simchildarr(i,t)==0 %Nothing happens
                simoffer(i,t) = 1;
                simwm(i,t) = simwm(i,t-1);  % He is offered the same wage
                simwf(i,t) = simwf(i,t-1);  % She also.
                simnkids(i,t) = simnkids(i,t-1);
                simat(i,t) = min(simat(i,t-1)+1,numel(at));
                 simstate(i,t) = (simat(i,t)-1)*numel(nt)*numel(wf)*numel(wm)+...
                  (simnkids(i,t))*numel(wf)*numel(wm)+(simwf(i,t)-1)*numel(wm)+simwm(i,t);
                simwork(i,t) = simwork(i,t-1);
                simreswagem(i,t) = simwm(i,t);
                simreswagef(i,t) = simwf(i,t);
            else
                disp('error')
                
            end
            
        elseif simwork(i,t-1)==2 %Only male worked in the previous period
            
            if simofferf(i,t)==1 && simdestm(i,t)==0 %Female got offer
                 simoffer(i,t) = 1;
                 simwm(i,t) = simwm(i,t-1);  % He is offered the same wage
                 simwf(i,t) = simwofferf(i,t);  % She got a new offer
                 simnkids(i,t) =simnkids(i,t-1);
                 simat(i,t) = min(simat(i,t-1)+1,numel(at));
                 simstate(i,t) = (simat(i,t)-1)*numel(nt)*numel(wf)*numel(wm)+...
                  (simnkids(i,t))*numel(wf)*numel(wm)+(simwf(i,t)-1)*numel(wm)+simwm(i,t);
                %State at evolves by one. Everything else stays the same.

                 simwork(i,t) = Pmfmat(simstate(i,t),t); %Choose who to work
                 simreswagem(i,t) = reswagemue(simstate(i,t),t);
                 simreswagef(i,t) = reswagefeu(simstate(i,t),t);

                 
            elseif simofferf(i,t)==1 && simdestm(i,t)==1 %Female offer & male destroy
                 simoffer(i,t) = 3;    %Only female offer
                 simwm(i,t) = simwofferm(i,t);  % Irrelevant:
                 simwf(i,t) = simwofferf(i,t);  % She got a new offer
                 simnkids(i,t) =simnkids(i,t-1);
                 simat(i,t) = min(simat(i,t-1)+1,numel(at));
                 simstate(i,t) = (simat(i,t)-1)*numel(nt)*numel(wf)*numel(wm)+...
                  (simnkids(i,t))*numel(wf)*numel(wm)+(simwf(i,t)-1)*numel(wm)+simwm(i,t);
                %State at evolves by one. Everything else stays the same.

                 simwork(i,t) = Pfmat(simstate(i,t),t); %Femal choose whether to work
                 simreswagem(i,t) = reswagemuu(simstate(i,t),t);
                 simreswagef(i,t) = reswagefuu(simstate(i,t),t);
                 
            elseif simofferf(i,t)==0 && simdestm(i,t)==1 %Male destroy
                 simoffer(i,t) = 4;    %No offer
                 simwm(i,t) = simwofferm(i,t);  % Irrelevant
                 simwf(i,t) = simwofferf(i,t);  % Irrelevant
                 simnkids(i,t) =simnkids(i,t-1);
                 simat(i,t) = min(simat(i,t-1)+1,numel(at));
                 simstate(i,t) = (simat(i,t)-1)*numel(nt)*numel(wf)*numel(wm)+...
                  (simnkids(i,t))*numel(wf)*numel(wm)+(simwf(i,t)-1)*numel(wm)+simwm(i,t);
                %State at evolves by one. Everything else stays the same.

                 simwork(i,t) = Pumat(simstate(i,t),t); 
                 simreswagem(i,t) = reswagemuu(simstate(i,t),t);
                 simreswagef(i,t) = reswagefuu(simstate(i,t),t);
                 
            elseif   simchildarr(i,t)==1 %Child arrives
                simoffer(i,t) = 2;         %Offer status remains the same.
                simwm(i,t) = simwm(i,t-1);  % He is offered the same wage
                simwf(i,t) = simwofferf(i,t);  % Irrelevant.
                simnkids(i,t) =  min(simnkids(i,t-1)+Pbeumat(simstate(i,t-1),t)-1,numel(nt));
                     %If give birth, add 1 to "simnkids".
                     
                if Pbeumat(simstate(i,t-1),t)== 1
                    % State evolution depends on choices. If not giving a birth
                    simat(i,t) = min(simat(i,t-1)+1,numel(at));
                    simstate(i,t) = (simat(i,t)-1)*numel(nt)*numel(wf)*numel(wm)+...
                  (simnkids(i,t))*numel(wf)*numel(wm)+(simwf(i,t)-1)*numel(wm)+simwm(i,t);
                       
                    simwork(i,t) = simwork(i,t-1); %Work status not changing
                    simreswagem(i,t) = simwm(i,t);
                   simreswagef(i,t) = reswagefeu(simstate(i,t),t);

                elseif Pbeumat(simstate(i,t-1),t)== 2
                    % If giving a birth,
                 simat(i,t) = 1; %at re-set to 1.
                simstate(i,t) = (simat(i,t)-1)*numel(nt)*numel(wf)*numel(wm)+...
                  (simnkids(i,t))*numel(wf)*numel(wm)+(simwf(i,t)-1)*numel(wm)+simwm(i,t);
                
                 simwork(i,t)=simwork(i,t-1);%Pwbeumat(simstate(i,t-1),t); %Work choice not 
                 %involved. So simwork(i,t)=simwork(i,t-1) should work.
                 %Here I let them choose to check if the actually choose to
                 %stay (validity check).
                 simreswagem(i,t) = simwm(i,t);
                  simreswagef(i,t) = reswagefeu(simstate(i,t),t);
              
                end
                
            elseif simchildarr(i,t)==0 %Nothing happens
                simoffer(i,t) = 2;
                simwm(i,t) = simwm(i,t-1);  % He is offered the same wage
                simwf(i,t) = simwofferf(i,t);  % Irrelevant
                simnkids(i,t) = simnkids(i,t-1);
                simat(i,t) = min(simat(i,t-1)+1,numel(at));
                 simstate(i,t) = (simat(i,t)-1)*numel(nt)*numel(wf)*numel(wm)+...
                  (simnkids(i,t))*numel(wf)*numel(wm)+(simwf(i,t)-1)*numel(wm)+simwm(i,t);
                simwork(i,t) = simwork(i,t-1);
                simreswagem(i,t) = simwm(i,t);
                simreswagef(i,t) = reswagefeu(simstate(i,t),t);
                
            else
                disp('error')
            end
            
        elseif simwork(i,t-1)==3 %Only female worked in the previous period
            
            if simofferm(i,t)==1 && simdestf(i,t)==0 %Male got offer
                 simoffer(i,t) = 1;
                 simwm(i,t) = simwofferm(i,t);  % He got an offer
                 simwf(i,t) = simwf(i,t-1);  % She is offered the same wage
                 simnkids(i,t) =simnkids(i,t-1);
                 simat(i,t) = min(simat(i,t-1)+1,numel(at));
                 simstate(i,t) = (simat(i,t)-1)*numel(nt)*numel(wf)*numel(wm)+...
                  (simnkids(i,t))*numel(wf)*numel(wm)+(simwf(i,t)-1)*numel(wm)+simwm(i,t);
                %State at evolves by one. Everything else stays the same.

                 simwork(i,t) = Pmfmat(simstate(i,t),t); %Choose who to work
                 simreswagem(i,t) = reswagemue(simstate(i,t),t);
                 simreswagef(i,t) = reswagefeu(simstate(i,t),t);
                 
            elseif simofferm(i,t)==1 && simdestf(i,t)==1 %Male offer & Female destroy
                 simoffer(i,t) = 2;    %Only Male offer
                 simwm(i,t) = simwofferm(i,t);  % He got a new offer:
                 simwf(i,t) = simwofferf(i,t);  % Irrelevant
                 simnkids(i,t) =simnkids(i,t-1);
                 simat(i,t) = min(simat(i,t-1)+1,numel(at));
                 simstate(i,t) = (simat(i,t)-1)*numel(nt)*numel(wf)*numel(wm)+...
                  (simnkids(i,t))*numel(wf)*numel(wm)+(simwf(i,t)-1)*numel(wm)+simwm(i,t);
                %State at evolves by one. Everything else stays the same.

                 simwork(i,t) = Pmmat(simstate(i,t),t); %Femal choose whether to work
                 simreswagem(i,t) = reswagemuu(simstate(i,t),t);
                 simreswagef(i,t) = reswagefuu(simstate(i,t),t);
                 
                 
            elseif simofferm(i,t)==0 && simdestf(i,t)==1 %Female destroy
                 simoffer(i,t) = 4;    %No offer
                 simwm(i,t) = simwofferm(i,t);  % Irrelevant
                 simwf(i,t) = simwofferf(i,t);  % Irrelevant
                 simnkids(i,t) =simnkids(i,t-1);
                 simat(i,t) = min(simat(i,t-1)+1,numel(at));
                 simstate(i,t) = (simat(i,t)-1)*numel(nt)*numel(wf)*numel(wm)+...
                  (simnkids(i,t))*numel(wf)*numel(wm)+(simwf(i,t)-1)*numel(wm)+simwm(i,t);
                %State at evolves by one. Everything else stays the same.

                 simwork(i,t) = Pumat(simstate(i,t),t); 
                 simreswagem(i,t) = reswagemuu(simstate(i,t),t);
                 simreswagef(i,t) = reswagefuu(simstate(i,t),t);
                 
            elseif   simchildarr(i,t)==1 %Child arrives
                simoffer(i,t) = 3;         %Offer status remains the same.
                simwm(i,t) = simwofferm(i,t);  % Irrelevant.
                simwf(i,t) = simwf(i,t-1);  % She is offered the same wage.
                simnkids(i,t) =  min(simnkids(i,t-1)+Pbuemat(simstate(i,t-1),t)-1,numel(nt));
                     %If give birth, add 1 to "simnkids".
                     
                if Pbuemat(simstate(i,t-1),t)== 1
                    % State evolution depends on choices. If not giving a birth
                    simat(i,t) = min(simat(i,t-1)+1,numel(at));
                    simstate(i,t) = (simat(i,t)-1)*numel(nt)*numel(wf)*numel(wm)+...
                  (simnkids(i,t))*numel(wf)*numel(wm)+(simwf(i,t)-1)*numel(wm)+simwm(i,t);
                       
                    simwork(i,t) = simwork(i,t-1); %Work status not changing
                    simreswagem(i,t) = reswagemuu(simstate(i,t),t);
                    simreswagef(i,t) = simwf(i,t);

                elseif Pbuemat(simstate(i,t-1),t)== 2
                    % If giving a birth,
                 simat(i,t) = 1; %at re-set to 1.
                simstate(i,t) = (simat(i,t)-1)*numel(nt)*numel(wf)*numel(wm)+...
                  (simnkids(i,t))*numel(wf)*numel(wm)+(simwf(i,t)-1)*numel(wm)+simwm(i,t);
                
                 simwork(i,t)=simwork(i,t-1);%Pwbuemat(simstate(i,t-1),t); %Work choice not 
                 %involved. So simwork(i,t)=simwork(i,t-1) should work.
                 %Here I let them choose to check if the actually choose to
                 %stay (validity check).
                 simreswagem(i,t) = reswagemuu(simstate(i,t),t);
                 simreswagef(i,t) = simwf(i,t);
              
                end
                
            elseif simchildarr(i,t)==0 %Nothing happens
                simoffer(i,t) = 3;
                simwm(i,t) = simwofferm(i,t);  % Irrelevant
                simwf(i,t) = simwf(i,t-1);  % same wage
                simnkids(i,t) = simnkids(i,t-1);
                simat(i,t) = min(simat(i,t-1)+1,numel(at));
                 simstate(i,t) = (simat(i,t)-1)*numel(nt)*numel(wf)*numel(wm)+...
                  (simnkids(i,t))*numel(wf)*numel(wm)+(simwf(i,t)-1)*numel(wm)+simwm(i,t);
                simwork(i,t) = simwork(i,t-1);
                simreswagem(i,t) = reswagemuu(simstate(i,t),t);
                simreswagef(i,t) = simwf(i,t);
            else
                disp('error')
            end
            
        elseif simwork(i,t-1)==4 % None worked in the previous period
            
            if simofferm(i,t)==1 && simofferf(i,t)==1 %Both got offer
                 simoffer(i,t) = 1;
                 simwm(i,t) = simwofferm(i,t);  % He got an offer
                 simwf(i,t) = simwofferf(i,t);  % She also
                 simnkids(i,t) =simnkids(i,t-1);
                 simat(i,t) = min(simat(i,t-1)+1,numel(at));
                 simstate(i,t) = (simat(i,t)-1)*numel(nt)*numel(wf)*numel(wm)+...
                  (simnkids(i,t))*numel(wf)*numel(wm)+(simwf(i,t)-1)*numel(wm)+simwm(i,t);
               

                 simwork(i,t) = Pmfmat(simstate(i,t),t); %Choose who to work
                 simreswagem(i,t) = reswagemuu(simstate(i,t),t);
                 simreswagef(i,t) = reswagefuu(simstate(i,t),t);
                 
            elseif simofferm(i,t)==1 && simofferf(i,t)==0 %Male offer
                 simoffer(i,t) = 2;    %Only Male offer
                 simwm(i,t) = simwofferm(i,t);  % He got a new offer:
                 simwf(i,t) = simwofferf(i,t);  % Irrelevant
                 simnkids(i,t) =simnkids(i,t-1);
                 simat(i,t) = min(simat(i,t-1)+1,numel(at));
                 simstate(i,t) = (simat(i,t)-1)*numel(nt)*numel(wf)*numel(wm)+...
                  (simnkids(i,t))*numel(wf)*numel(wm)+(simwf(i,t)-1)*numel(wm)+simwm(i,t);
                %State at evolves by one. Everything else stays the same.

                 simwork(i,t) = Pmmat(simstate(i,t),t); %Femal choose whether to work
                 simreswagem(i,t) = reswagemuu(simstate(i,t),t);
                 simreswagef(i,t) = reswagefuu(simstate(i,t),t);
                 
            elseif simofferm(i,t)==0 && simofferf(i,t)==1 %Female offer
                 simoffer(i,t) = 3;    %No offer
                 simwm(i,t) = simwofferm(i,t);  % Irrelevant
                 simwf(i,t) = simwofferf(i,t);  % She got an offer
                 simnkids(i,t) =simnkids(i,t-1);
                 simat(i,t) = min(simat(i,t-1)+1,numel(at));
                 simstate(i,t) = (simat(i,t)-1)*numel(nt)*numel(wf)*numel(wm)+...
                  (simnkids(i,t))*numel(wf)*numel(wm)+(simwf(i,t)-1)*numel(wm)+simwm(i,t);
                %State at evolves by one. Everything else stays the same.

                 simwork(i,t) = Pfmat(simstate(i,t),t); 
                 simreswagem(i,t) = reswagemuu(simstate(i,t),t);
                 simreswagef(i,t) = reswagefuu(simstate(i,t),t);
                 
            elseif   simchildarr(i,t)==1 %Child arrives
                simoffer(i,t) = 4;         %Offer status remains the same.
                simwm(i,t) = simwofferm(i,t);  % Irrelevant.
                simwf(i,t) = simwofferf(i,t-1);  % Irrelevant.
                simnkids(i,t) =  min(simnkids(i,t-1)+Pbuumat(simstate(i,t-1),t)-1,numel(nt));
                     %If give birth, add 1 to "simnkids".
                     
                if Pbuumat(simstate(i,t-1),t)== 1
                    % State evolution depends on choices. If not giving a birth
                    simat(i,t) = min(simat(i,t-1)+1,numel(at));
                    simstate(i,t) = (simat(i,t)-1)*numel(nt)*numel(wf)*numel(wm)+...
                  (simnkids(i,t))*numel(wf)*numel(wm)+(simwf(i,t)-1)*numel(wm)+simwm(i,t);
                       
                    simwork(i,t) = simwork(i,t-1); %Work status not changing
                    simreswagem(i,t) = reswagemuu(simstate(i,t),t);
                    simreswagef(i,t) = reswagefuu(simstate(i,t),t);
                 
                elseif Pbuumat(simstate(i,t-1),t)== 2
                    % If giving a birth,
                 simat(i,t) = 1; %at re-set to 1.
                 simstate(i,t) = (simat(i,t)-1)*numel(nt)*numel(wf)*numel(wm)+...
                  (simnkids(i,t))*numel(wf)*numel(wm)+(simwf(i,t)-1)*numel(wm)+simwm(i,t);
                
                 simwork(i,t)=simwork(i,t-1);
                 simreswagem(i,t) = reswagemuu(simstate(i,t),t);
                 simreswagef(i,t) = reswagefuu(simstate(i,t),t);
              
                end
                
            elseif simchildarr(i,t)==0 %Nothing happens
                simoffer(i,t) = 4;
                simwm(i,t) = simwofferm(i,t);  % Irrelevant
                simwf(i,t) = simwofferf(i,t-1);  % Irrelevant
                simnkids(i,t) = simnkids(i,t-1);
                simat(i,t) = min(simat(i,t-1)+1,numel(at));
                simstate(i,t) = (simat(i,t)-1)*numel(nt)*numel(wf)*numel(wm)+...
                  (simnkids(i,t))*numel(wf)*numel(wm)+(simwf(i,t)-1)*numel(wm)+simwm(i,t);
                simwork(i,t) = simwork(i,t-1);
                simreswagem(i,t) = reswagemuu(simstate(i,t),t);
                simreswagef(i,t) = reswagefuu(simstate(i,t),t);
                
            else
                disp('error')
            end
        end
        
        %Things used in the next section
        if simnkids(i,t-1)<simnkids(i,t) && simwork(i,t-1)==1 && simwork(i,t)==3
            quitbirthm(i,t)=1;
        elseif simnkids(i,t-1)<simnkids(i,t) && simwork(i,t-1)==1 && simwork(i,t)==2
            quitbirthf(i,t)=1;
        end
        
        if (simwork(i,t-1)==2||simwork(i,t-1)==4)&&(simwork(i,t)==1||simwork(i,t)==3)
            comebackf(i,t)=1;
        elseif (simwork(i,t-1)==3||simwork(i,t-1)==4)&&(simwork(i,t)==1||simwork(i,t)==2)
            comebackm(i,t)=1;
        end
        
    end
    
    if sum(comebackm(i,:))>0&&sum(quitbirthm(i,:))>0&&sum(quitbirthm(i,:))<9
        aux=find(quitbirthm(i,:));
        for j=1:numel(aux)
            aux2=aux(j):period;
            aux3=min(intersect(find(comebackm(i,:)),aux2));
            if isempty(aux3)==0
                quitdurm(i,j)=aux3-aux(j);
                wagedifm(i,j)=simwm(i,aux(j))-simwm(i,aux3);
                quitagem(i,j)=aux(j);
            else
                quitdurm(i,j)=99;
                wagedifm(i,j)=99;
                quitagem(i,j)=aux(j);
            end
        end
    elseif sum(quitbirthm(i,:))>8
        i %If one family experiences more than one child-quit, this doesn't work.
    end
    clear aux aux2 aux3
     if sum(comebackf(i,:))>0&&sum(quitbirthf(i,:))>0&&sum(quitbirthf(i,:))<9
        aux=find(quitbirthf(i,:));
        for j=1:numel(aux)
            aux2=aux(j):period;
            aux3=min(intersect(find(comebackf(i,:)),aux2));
            if isempty(aux3)==0
                quitdurf(i,j)=aux3-aux(j);
                wagediff(i,j)=simwf(i,aux(j))-simwf(i,aux3);
                quitagef(i,j)=aux(j);
            else
                quitdurf(i,j)=99;
                wagediff(i,j)=99;
                quitagef(i,j)=aux(j);
            end
        end
    elseif sum(quitbirthf(i,:))>8
        i %If one family experiences more than one child-quit, this doesn't work.
     end
    clear aux aux2 aux3

end

        
%%%%%       

%Simulation done. Calculate moments.
%"Conditional on age" means conditional on "young" (20-29) and "old"
%(30-39).

%%%%%

%Observed moments
birthrateyobs=
birthrateoobs=

kids0obs=0.2636;
kids1obs=0.1889;
kids2obs=0.2873;
kids3obs=0.1567;
kids4obs=0.0661;

pquitmyobs=0.13532;
pquitmoobs=0.06684;
pquitfyobs=0.4862;
pquitfoobs=0.426;

meanquitdurmyobs=0.9173; %=47.7/52
meanquitdurmoobs=0.7669; %=39.88/52
meanquitdurfyobs=1.2906; %=67.11/52
meanquitdurfoobs=1.0067; %=52.35/52

empmkidageobs=[ , , , , ];
empfkidageobs=[ , , , , ];

meanwagemykidobs=[ , , ];
meanwagefykidobs=[ , , ];
meanwagemokidobs=[ , , ];
meanwagefokidobs=[ , , ];

varwagemykidobs=[ , , ];
varwagefykidobs=[ , , ];
varwagemokidobs=[ , , ];
varwagefokidobs=[ , , ];

%%%%%

%(1) Childbirth rate conditional on age
%(2) Probability of quitting conditional on childbirth, parents age and
%gender


%Childbirth rate
birthratey=sum(simnkids(:,(period/2)))./(nsim*period/2);
birthrateo=(sum(simnkids(:,(period)))-sum(simnkids(:,(period/2))))./(nsim*period/2);
moment1=[birthratey,birthrateo]-[birthrateyobs,birthrateoobs];

%Probability of quitting with birth conditional on parents' age and gender
pquitmy=sum(sum(quitbirthm(:,1:(period/2))))/sum(simnkids(:,(period/2)));
pquitmo=sum(sum(quitbirthm(:,(period/2+1):period)))/(sum(simnkids(:,(period)))-sum(simnkids(:,(period/2))));
pquitfy=sum(sum(quitbirthf(:,1:(period/2))))/sum(simnkids(:,(period/2)));
pquitfo=sum(sum(quitbirthf(:,(period/2+1):period)))/(sum(simnkids(:,(period)))-sum(simnkids(:,(period/2))));

moment2=[pquitmy,pquitmo,pquitfy,pquitfo]-[pquitmyobs,pquitmoobs,pquitfyobs,pquitfoobs];




% (3) Leaving durations, conditional on childbirth quitting
meanquitdurmy=mean(quitdurm(quitdurm>0&quitdurm<10&quitagem<period/2));
meanquitdurmo=mean(quitdurm(quitdurm>0&quitdurm<10&quitagem>=period/2));
meanquitdurfy=mean(quitdurf(quitdurf>0&quitdurf<10&quitagef<period/2));
meanquitdurfo=mean(quitdurf(quitdurf>0&quitdurf<10&quitagef>=period/2));

moment3=[meanquitdurmy,meanquitdurmo,meanquitdurfy,meanquitdurfo]-[meanquitdurmyobs,meanquitdurmoobs,meanquitdurfyobs,meanquitdurfoobs];





%(5) Employment rate conditional on child's age

empmkidage=zeros(numel(at),1);
empfkidage=zeros(numel(at),1);

for n = 1:numel(at)
aux=zeros(nsim,period);
aux((simwork==1|simwork==2)&simat==n&simnkids>0)=1;
aux2=zeros(nsim,period);
aux2((simwork==1|simwork==3)&simat==n&simnkids>0)=1;
aux3=zeros(nsim,period);
aux3(simat==n&simnkids>0)=1;
empmkidage(n)=sum(sum(aux))./sum(sum(aux3));
empfkidage(n)=sum(sum(aux2))./sum(sum(aux3));
end
clear aux aux2 aux3

moment5=[empmkidage.', empfkidage.']-[empmkidageobs,empfkidageobs];


%(6) Wage distribution (mean and var) conditional on child's age and
%parent's age
%FOR THIS MOMENT CONDITION, DEFINE YOUNG<35, OLD>=35
meanwagemykid=zeros(numel(at),1);
meanwagefykid=zeros(numel(at),1);
varwagemykid=zeros(numel(nt),1);
varwagefykid=zeros(numel(nt),1);
meanwagemokid=zeros(numel(nt),1);
meanwagefokid=zeros(numel(nt),1);
varwagemokid=zeros(numel(nt),1);
varwagefokid=zeros(numel(nt),1);

for n = 1:numel(nt)
aux4=zeros(nsim,period);
aux4((simwork==1|simwork==2)&simnkids==(n-1))=simwm((simwork==1|simwork==2)&simnkids==(n-1));
aux5=aux4>0;
aux6=zeros(nsim,period);
aux6((simwork==1|simwork==3)&simnkids==(n-1))=simwf((simwork==1|simwork==3)&simnkids==(n-1));
aux7=aux6>0;
aux8=aux4(:,15:period);
aux9=aux6(:,15:period);
meanwagemykid(n)=mean(wm(aux4(aux5(:,1:14))));
varwagemykid(n)=var(wm(aux4(aux5(:,1:14))));
meanwagefykid(n)=mean(wf(aux6(aux7(:,1:14))));
varwagefykid(n)=var(wf(aux6(aux7(:,1:14))));
meanwagemokid(n)=mean(wm(aux8(aux5(:,15:period))));
varwagemokid(n)=var(wm(aux8(aux5(:,15:period))));
meanwagefokid(n)=mean(wf(aux9(aux7(:,15:period))));
varwagefokid(n)=var(wf(aux9(aux7(:,15:period))));
end

clear aux4 aux5 aux6 aux7 aux8 aux9

moment6=[meanwagemykid(1:3).',meanwagefykid(1:3).',meanwagemokid(1:3).',meanwagefokid(1:3).',...
    varwagemykid(1:3).',varwagefykid(1:3).',varwagemokid(1:3).',varwagefokid(1:3).']...
    -[meanwagemykidobs, meanwagefykidobs, meanwagemokidobs, meanwagefokidobs,...
    varwagemykidobs, varwagefykidobs, varwagemokidobs, varwagefokidobs];


% (7) Number of children at the end of career 
kids0=numel(simnkids(:,20)==0)/nsim;
kids1=numel(simnkids(:,20)==1)/nsim;
kids2=numel(simnkids(:,20)==2)/nsim;
kids3=numel(simnkids(:,20)==3)/nsim;
kids4=numel(simnkids(:,20)==4)/nsim;

moment7=[kids0,kids1,kids2,kids3,kids4]-[kids0obs,kids1obs,kids2obs,kids3obs,kids4obs];


%Define objective function

setmoment = [moment1,moment2,moment3,moment5,moment6,moment7];
smmobjective = setmoment*setmoment.';

end


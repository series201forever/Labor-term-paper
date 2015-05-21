function[q,works,wages]=simul(parahat,epsstarhat,epsstarstarhat,nind,nsim,data)
rng(278);
%This function construct matrices of simulated in-sample experience,
%work/not, and wages. The size of matrices is number of
%individual*period*number of simulations.
%The inputs are estimated parameters, two thresholds at estimated
%parameters, number of individuals, number of simulations and data.

%Draw sequence of wages and measurement errors.
sigmaepshat=parahat(7);
sigmaetahat=parahat(8);
simeps=sigmaepshat*randn([nind,15,nsim]);
simeta=sigmaetahat*randn([nind,15,nsim]);

%Construct matrices to store results. q has one additional dimension
%because I need to refer to q(16) in calculating period 15 value.
q=zeros(nind,16,nsim);
works=zeros(nind,15,nsim);
wages=zeros(nind,15,nsim);

%Forward simulation. Note that epsstar and epsstarstar(i,1,s) corresponds
%to period 0, while q=0 corresponds to period 0.
 for i=1:nind
    s=data(1+15*(i-1),7);
    r=data(1+15*(i-1),8);
    inite=data(1+15*(i-1),6);
    for k=1:nsim
        %foreach simulation, starting from period 1, simulate the decision
        %forward.
        for t=1:15
            if epsstarhat(i,q(i,t,k)+1,t)>epsstarstarhat(i,q(i,t,k)+1,t)&&simeps(i,t,k)<epsstarhat(i,q(i,t,k)+1,t)
            works(i,t,k)=0;
            q(i,t+1,k)=q(i,t,k);
            wages(i,t,k)=0;
            elseif epsstarhat(i,q(i,t,k)+1,t)>epsstarstarhat(i,q(i,t,k)+1,t)&&simeps(i,t,k)>epsstarhat(i,q(i,t,k)+1,t)
            q(i,t+1,k)=q(i,t,k)+1;
            works(i,t,k)=1;
            wages(i,t,k)=exp(parahat(1)+parahat(2)*s+parahat(3)*(inite+q(i,t,k))+parahat(4)*(inite+q(i,t,k))*r+parahat(9)*(inite+q(i,t,k))^2+simeps(i,t,k)+simeta(i,t,k));
            elseif epsstarhat(i,q(i,t,k)+1,t)<epsstarstarhat(i,q(i,t,k)+1,t)&&simeps(i,t,k)<epsstarstarhat(i,q(i,t,k)+1,t)
            q(i,t+1,k)=q(i,t,k);
            wages(i,t,k)=0;
            works(i,t,k)=0;
            elseif epsstarhat(i,q(i,t,k)+1,t)<epsstarstarhat(i,q(i,t,k)+1,t)&&simeps(i,t,k)>epsstarstarhat(i,q(i,t,k)+1,t)
            q(i,t+1,k)=q(i,t,k)+1;
            works(i,t,k)=1;
            wages(i,t,k)=exp(parahat(1)+parahat(2)*s+parahat(3)*(inite+q(i,t,k))+parahat(4)*(inite+q(i,t,k))*r+parahat(9)*(inite+q(i,t,k))^2+simeps(i,t,k)+simeta(i,t,k));
            end
        end
    end
 end
end


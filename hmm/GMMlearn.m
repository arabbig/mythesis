function [m,sigma2,stran,phtgV1T]=GMMlearn(v,S,Tskip,opts)
%SARLEARN EM training of a Switching AR model
% [a,sigma2,stran,phtgV1T]=SARlearn(v,L,S,opts)
%
% Inputs:
% v :  a single timeseries is contained in the row vector v
% L :  order of tha AR model 
% S : number of AR models.
% Tskip forces the switches to make a transition only at times t for mod(t,Tskip)==0
% opts.maxit
% opts.plotprogress
% 
% Outputs:
% a : learned AR coefficients
% sigma2 : learned innovation noise
% stran : learned transition distribution
% phtgV1T : smoothed posterior p(h(t)|v(1:T))
% See also demoSARinference.m
import brml.*
m=zeros(S,1); % set the AR coefficients
stran=condp(ones(S,S)); % switch transition
sprior=condp(ones(S,1)); % switch prior
sigma2=var(v)*ones(1,S);
T=size(v,2);

loglikem = zeros(1,opts.maxit);
for emloop=1:opts.maxit
    % Inference using HMM structure:
    [logalpha,loglik]=HMMforwardGMM(v,stran,sprior,m,sigma2,Tskip);
    logbeta=HMMbackwardGMM(v,stran,m,sigma2,Tskip);
    [phtgV1T,phthtpgV1T]=HMMsmoothGMM(logalpha,logbeta,m,sigma2,stran,v,Tskip);

    loglikem(emloop)=loglik;
    if opts.plotprogress; plot(loglikem); title('log likelihood'); drawnow; end
    for s=1:S
        m_sum=0; m_num=0; sigma_sum=0; sigma_num=0;
        for t=2:T
            sigma_sum = sigma_sum+phtgV1T(s,t)*(v(t)-v(t-1)).^2;
            sigma_num = sigma_num + phtgV1T(s,t);
            m_sum = m_sum+phtgV1T(s,t)*v(t);
            m_num = m_num + phtgV1T(s,t);
        end
        sigma2(s)=sigma_sum/sigma_num;
        m(s)=m_sum/m_num; 
    end
    t=1:T-1; tt=t(find(mod(t+1,Tskip)==0));
    stran=condp(sum(phthtpgV1T(:,:,tt),3)');
end

end
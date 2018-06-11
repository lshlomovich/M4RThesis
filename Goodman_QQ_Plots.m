function theoretical_values = Goodman_QQ_Plots(mean,observed_values,kk)

L = length(observed_values);

pp=(1:L)./(L+1);

% find quantile of Goodman dist for coherence
% for the AR rho^2=0.36 for all frequencies!

theoretical_values = zeros(1,L);

for j=1:L
    prob=pp(j);
    theoretical_values(j) = fzero(@(x) estcoh(x,kk,2,prob,mean),[0.000001,0.999]);
end
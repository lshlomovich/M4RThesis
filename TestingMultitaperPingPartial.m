%here we test the multitapering function for partial coherence
%we now can generate random poisson processes

%ITERATIVE APPROACH

% arrival rate per second
lambdaAB=1/100;  %X1 rate
lambdaBC=1/100;  %X2 rate
lambdaCA=1/100; %X3 rate

%make an array of lambda values

%lambdas = [lambdaAB; lambdaBC
lambdas = [lambdaAB; lambdaBC; lambdaCA];
N = length(lambdas);

%simulation time in seconds (1 hour)
T = 10*3600; 

%we will bin the counts to 0:1/m:1, ie rate per unit time is lambda*T/m
m = 10000;

n=500;
f=1500;
R_partial = ones(n,3);
R_partialest = zeros(n,3);
ntapers = 10;
for i=1:n
    %times of events on (0,1) open interval
    AB = transpose(rand(1,poissrnd(lambdaAB*T)));
    BC = transpose(rand(1,poissrnd(lambdaBC*T)));
    CA = transpose(rand(1,poissrnd(lambdaCA*T)));
    
    %add a ping - regular pulse
    %f gives frequency so that delta = f/m
    f1 = 10;
    delta1 = f1/m;
    %set a random start time t_0 between 0 and delta
    t0a = rand(1)*delta1;
    
    times1 = transpose(t0a:(1/m * 10):1);
    
    %select which of these show by bernoulli trial
    %probability of a 0 given by prob0
    prob01 = 0.1; %100*(1-p)% of spikes included from original impulse train
    ping1 = rand((length(times1)),1);
    ping1(ping1 < prob01) = 0;
    ping1(ping1 >= prob01) = 1;
    
    prob02 = 0.1; %100*(1-p)% of spikes included from original impulse train
    ping2 = rand((length(times1)),1);
    ping2(ping2 < prob01) = 0;
    ping2(ping2 >= prob01) = 1;
    
    
    %set up new ping times array
    indexes1 = find(ping1); %this finds indices of non-zero entries
    bernpingtimes1 = zeros(length(indexes1),1);
    indexes2 = find(ping2); %this finds indices of non-zero entries
    bernpingtimes2 = zeros(length(indexes2),1);
    
    
    for j = 1:(length(indexes1))
        bernpingtimes1(j) = times1(indexes1(j));
    end

    for jj = 1:(length(indexes2))
        bernpingtimes2(jj) = times1(indexes2(jj));
    end

    
    %we can now bin these times
    
    %we will bin the counts to 0:1/m:1, ie rate per unit time is lambda*T/m
    m = 10000;
    lambdas = lambdas*T/m;
    
    %this generates independent poisson processes, X_1, X_2, X_3
    X1times = [AB;bernpingtimes1];
    X2times = [BC;times1];
    X3times = [CA;bernpingtimes2];
    
    X1freq = transpose(hist(X1times,0:(1/m):1));
    X2freq = transpose(hist(X2times,0:(1/m):1));
    X3freq = transpose(hist(X3times,0:(1/m):1));
    
    %put these frequencies in an array
    freqs = [X1freq, X2freq, X3freq];
    
    %initialise counts
    counts = zeros(m+1,N);
    
    %counts contains the zero mean frequencies
    for k=1:N
        counts(:,k) = freqs(:,k) - mean(freqs(:,k));
    end
    
    [X_ests, X_cross_ests, X_matrix] = spectrum_est(counts(:,1),counts(:,2),counts(:,3),ntapers);
    %upweighting
    a1 = max(X_ests(:,1));
    a2 = max(X_ests(:,2));
    a3 = max(X_ests(:,3));

    A = [a1 0 0; 0 a2 0; 0 0 a3];
    xi = 0.000010;
    B = xi * A;
    X_matrix = X_matrix + B;
    
    inverted = inv(X_matrix(:,:,f));
    %iwishf = iwishrnd(abs(inverted),10);
    
    %use values to compute partial coherence
    
    R_partial(i,1) = (abs(inverted(1,2)).^2)./(inverted(1,1).*inverted(2,2));
    R_partial(i,2) = (abs(inverted(1,3)).^2)./(inverted(1,1).*inverted(3,3));
    R_partial(i,3) = (abs(inverted(2,3)).^2)./(inverted(2,2).*inverted(3,3));
    %R_partialest(i,1) = (abs(iwishf(1,2).^2))./(iwishf(1,1).*iwishf(2,2));
    %R_partialest(i,2) = (abs(iwishf(1,3).^2))./(iwishf(1,1).*iwishf(3,3));
    %R_partialest(i,3) = (abs(iwishf(2,3).^2))./(iwishf(2,2).*iwishf(3,3));
end

lvalues = 0:n;
Cs = 1-( 1- (( 0.95 ).^(1./(n-lvalues+1))) ).^(1./(ntapers-3+1));
sortedR = sort(R_partial(:,1));
for i=0:n-1
    if sortedR(n-i) < Cs(i+1)
        break
    end
end
critR = sortedR(n-i);
critC = Cs(i+1);

%R_partial = abs(R_partial);

%unbiased estimator - yields negative valued estimates 
%R_partial = (R_partial - (1/(ntapers - 3 + 2)))/(1- (1/(ntapers - 3 + 2)));
%set any negatives to zero - changes distribution fit 
%R_partial(R_partial<0) = 0;


%figure
%qqplot(R_partial(:,1),R_partialest(:,1))
%title('\rho_{12} estimate versus Wishart estimate')
%figure
%qqplot(R_partial(:,2),R_partialest(:,2))
%title('\rho_{13} estimate versus Wishart estimate')
%figure
%qqplot(R_partial(:,2),R_partialest(:,3))
%title('\rho_{23} estimate versus Wishart estimate')


partial_stat = ((ntapers-1).*R_partial(:,1))./(1-R_partial(:,1));
finvcdf = finv(linspace(0,1,n),2,2*(ntapers-1));
betainvcdf = betainv([1:n]./(n+1),1,ntapers-3+1);
figure
qqplot(partial_stat,finvcdf)
title('Estimate versus F distribution CDF')
figure
qqplot(R_partial(:,1),betainvcdf)
title('Estimate versus Beta distribution CDF')

partial12 = R_partial(:,1);
partial13 = R_partial(:,2);
partial23 = R_partial(:,3);
stat12 = ((ntapers-1).*partial12)./(1-partial12);
stat13 = ((ntapers-1).*partial13)./(1-partial13);
stat23 = ((ntapers-1).*partial23)./(1-partial23);
figure
plot(partial12,'b')
hold on
plot(partial13, 'g')
hold on
plot(partial23, 'k')
hold on
plot(repmat(critR,n,1),'r')
%plot(repmat(finv(0.95,2,2*(ntapers-1)),n,1),'r')
hold on
plot(repmat(critC,n,1))
hold off
title('Estimates verus Significance Line')
legend('\rho_{12}','\rho_{13}','\rho_{23}','Rsigline', 'Csigline')


%multitapering function for p=3 (3 times series)
function [spectrum_ests, cross_ests, spectral_mat] = spectrum_est(X_1,X_2,X_3,number_seq)
%N needs to be equal to length of sequence - sequences should be of same
%size when input 
N = length(X_1);
%Select a bandwidth, usually NW is 2, 5/2, 3, 7/2, or 4
time_halfbandwidth = number_seq/2;
%Obtain DPSSs
H = dpss(N, time_halfbandwidth, number_seq);
%Compute fourier transform
FS11 = zeros(N, number_seq);
FS22 = zeros(N, number_seq);
FS33 = zeros(N, number_seq);
FS12 = zeros(N, number_seq);
FS13 = zeros(N, number_seq);
FS23 = zeros(N, number_seq);
for i=1:number_seq
    %spectrum estimates
    FS11(:,i) = abs(fftshift(fft(H(:,i).*X_1))).^2;
    FS22(:,i) = abs(fftshift(fft(H(:,i).*X_2))).^2;
    FS33(:,i) = abs(fftshift(fft(H(:,i).*X_3))).^2;
    %cross-spectrum estimates
    FS12(:,i) = fftshift(fft(H(:,i).*X_1)).*conj(fftshift(fft(H(:,i).*X_2)));
    FS13(:,i) = fftshift(fft(H(:,i).*X_1)).*conj(fftshift(fft(H(:,i).*X_3)));
    FS23(:,i) = fftshift(fft(H(:,i).*X_2)).*conj(fftshift(fft(H(:,i).*X_3)));
end
FS11 = mean(FS11, 2);
FS22 = mean(FS22, 2);
FS33 = mean(FS33, 2);
FS12 = mean(FS12, 2);
FS13 = mean(FS13, 2);
FS23 = mean(FS23, 2);
%output array containing spectral estimates
spectrum_ests = [FS11, FS22, FS33];
cross_ests = [FS12, FS13, FS23];

spectral_mat = zeros(3,3,N);

spectral_mat(1,1,:) = FS11;
spectral_mat(2,2,:) = FS22;
spectral_mat(3,3,:) = FS33;

spectral_mat(1,2,:) = FS12;
spectral_mat(2,1,:) = conj(FS12);

spectral_mat(1,3,:) = FS13;
spectral_mat(3,1,:) = conj(FS13);

spectral_mat(2,3,:) = FS23;
spectral_mat(3,2,:) = conj(FS23);
    
end

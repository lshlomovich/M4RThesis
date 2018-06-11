%we now can generate random poisson processes
m = 10000;

% arrival rate per second
lambdaAB=1/35;  %X1 rate
lambdaBC=1/70;  %X2 rate
lambdaCA=1/50; %X3 rate

%make an array of lambda values
lambdas = [lambdaAB; lambdaBC; lambdaCA];
N = length(lambdas);

%simulation time in seconds (1 hour)
T = 10*3600; 

%times of events on (0,1) open interval 
AB = transpose(rand(1,poissrnd(lambdaAB*T)));
BC = transpose(rand(1,poissrnd(lambdaBC*T)));
CA = transpose(rand(1,poissrnd(lambdaCA*T)));

%add a ping - regular pulse
%f gives frequency so that delta = f1/100 gives periodicity 
f1 = 10;
delta1 = f1/m; 
%set a random start time t_0 between 0 and delta
t0a = rand(1)*delta1;
    
times1 = transpose(t0a:(1/m * f1):1);

%select which of these show by bernoulli trial 
%probability of a 0 given by prob0
prob01 = 0.1; %100*(1-p)% of spikes included from original impulse train
prob02 = 0.1; %second probability for independent bernoulli trial, use as equal to prob01 if only trialling once
ping1 = rand((length(times1)),1);
ping1(ping1 < prob01) = 0;
ping1(ping1 >= prob01) = 1;
ping2 = rand((length(times1)),1);
ping2(ping2 < prob02) = 0;
ping2(ping2 >= prob02) = 1;


%set up new ping times array
indexes1 = find(ping1); %this finds indices of non-zero entries
indexes2 = find(ping2);
bernpingtimes1 = zeros(length(indexes1),1);
bernpingtimes2 = zeros(length(indexes2),1);


for i = 1:(length(indexes1))
        bernpingtimes1(i) = times1(indexes1(i));
end

for i = 1:(length(indexes2))
        bernpingtimes2(i) = times1(indexes2(i));
end

%we can now bin these times

%we will bin the counts to 0:1/m:1, ie rate per unit time is lambda*T/m
lambdas = lambdas*T/m;

%this generates independent poisson processes, X_1, X_2, X_3
%and adds counts at the pulse times
X1times = [AB;bernpingtimes1];
X2times = [BC;times1];
X3times = [CA;bernpingtimes1];

X1freq = transpose(hist(X1times,0:(1/m):1));
X2freq = transpose(hist(X2times,0:(1/m):1));
X3freq = transpose(hist(X3times,0:(1/m):1));

%put these frequencies in an array
freqs = [X1freq, X2freq, X3freq];

%initialise counts
counts = zeros(m+1,N);

%counts contains the zero mean frequencies 
for i=1:N
    counts(:,i) = freqs(:,i) - mean(freqs(:,i));
end

%array of estimated spectral matrices in X_matrix
%X_ests contains the spectrum estimates, and X_cross_ests contains the
%cross spectral estimates
ntapers = 40;
[X_ests, X_cross_ests, X_matrix] = spectrum_est(counts(:,1),counts(:,2),counts(:,3),ntapers);

%for up-weighting access a_i = max_f(S_ii)
a1 = max(X_ests(:,1));
a2 = max(X_ests(:,2));
a3 = max(X_ests(:,3));

A = [a1 0 0; 0 a2 0; 0 0 a3];
xi = 0.000010;
B = xi * A;

%diagonal up-weighting
X_matrix = X_matrix + B;

S11 = X_ests(:,1) + B(1,1);
S22 = X_ests(:,2) + B(2,2);
S33 = X_ests(:,3) + B(3,3);
S12 = X_cross_ests(:,1);
S13 = X_cross_ests(:,2);
S23 = X_cross_ests(:,3);

%R_partial returns complex if this line commented
%X_matrix = real(X_matrix);

%calculate coherences
R12_est = (abs(S12.^2))./(S11.*S22);
R13_est = (abs(S13.^2))./(S11.*S33);
R23_est = (abs(S23.^2))./(S22.*S33);

%calculate partial coherences using inversion method

%col1 ~ partial coherence between 1 and 2
%col2 ~ partial coherence between 1 and 3
%col3 ~ partial coherence between 2 and 3
R_partial = zeros(m+1,3);

%use iterative method for p=3 case as matrix inversion is simply 1x1
for i=1:m+1
    %access spectral matrix for given f
    %invert and store
    inverted = inv(X_matrix(:,:,i));
    %use values to compute partial coherence
    R_partial(i,1) = (abs(inverted(1,2)).^2)./(inverted(1,1).*inverted(2,2));
    R_partial(i,2) = (abs(inverted(1,3)).^2)./(inverted(1,1).*inverted(3,3));
    R_partial(i,3) = (abs(inverted(2,3)).^2)./(inverted(2,2).*inverted(3,3));
end

%multiple hypothesis testing
lvalues = 0:m+1;
Cs = 1-( 1- (( 0.95 ).^(1./((m+1)-lvalues+1))) ).^(1./(ntapers-3+1));
%12
sortedR12 = sort(R_partial(:,1));
for i=0:m
    if sortedR12(m+1-i) < Cs(i+1)
        break
    end
end
critR12 = sortedR12((m+1)-i);

sortedR13 = sort(R_partial(:,2));
for i=0:m
    if sortedR13(m+1-i) < Cs(i+1)
        break
    end
end
critR13 = sortedR13((m+1)-i);
%critC = Cs(i+1);

sortedR23 = sort(R_partial(:,3));
for i=0:m
    if sortedR23(m+1-i) < Cs(i+1)
        break
    end
end
critR23 = sortedR23((m+1)-i);


% code to compute the partial coherence via the direct method
% R_partialdirect = zeros(m+1,3);
% % 
% S12c3 = S12 - S13.*(1./(S33)).*conj(S23);
% S13c2 = S13 - S12.*(1./(S22)).*S23;
% S23c1 = S23 - conj(S12).*(1./(S11)).*S13;
% % 
% S11c2 = S11 - S12.*conj(S12).*(1./S22);
% S11c3 = S11 - S13.*conj(S13).*(1./S33);
% S22c1 = S22 - conj(S12).*S12.*(1./S11);
% S22c3 = S22 - S23.*conj(S23).*(1./S33);
% S33c1 = S33 - conj(S13).*S13.*(1./S11);
% S33c2 = S33 - conj(S23).*S23.*(1./S22);
% 
% R_partialdirect(:,1) = ((abs(S12c3)).^2)./(S11c3.*S22c3);
% R_partialdirect(:,2) = ((abs(S13c2)).^2)./(S11c2.*S33c2);
% R_partialdirect(:,3) = ((abs(S23c1)).^2)./(S22c1.*S33c1);

S11theory = lambdas(1)*(1-lambdas(1)) + ((1-prob01)*delta1)*(1 - (1-prob02)*delta1);
S22theory = lambdas(2)*(1-lambdas(2)) + ((1-prob01)*delta1)*(1 - (1-prob02)*delta1);
S33theory = lambdas(3)*(1-lambdas(3)) + ((1-prob01)*delta1)*(1 - (1-prob02)*delta1);
S12theory = ((1-prob01)*delta1)*(1 - ((1-prob02)*delta1));
S13theory = ((1-prob01)*delta1)*(1 - ((1-prob02)*delta1));
S23theory = ((1-prob01)*delta1)*(1 - ((1-prob02)*delta1));

R12theory = S12theory^2/(S11theory*S22theory);
R13theory = S13theory^2/(S11theory*S33theory);
R23theory = S23theory^2/(S22theory*S33theory);

alpha = 0.95; %significance level
fcutoff = finv(alpha,2,2*(ntapers-1)); %significance cutoff
partial_significance = sqrt((fcutoff/(ntapers-1))/(1 + (fcutoff/(ntapers-1))));
bcutoff = betainv(alpha,1,ntapers-3+1);

%plot the result
figure
plot(linspace(-1/2,1/2,m+1),repmat(S11theory,m+1,1),'b')
hold on
plot(linspace(-1/2,1/2,m+1),S11,'r')
hold on
plot(linspace(-1/2,1/2,m+1),repmat(S22theory,m+1,1),'g')
hold on
plot(linspace(-1/2,1/2,m+1),S22,'k')
hold off
legend('S11theory', 'S11','S22theory','S22')
%saveas(gcf, 'bern.1ping.spectra.jpg')



%plot the result
figure
plot(linspace(-1/2,1/2,m+1),repmat(S12theory,m+1,1),'b')
hold on
plot(linspace(-1/2,1/2,m+1),real(S12),'r')
hold on
plot(linspace(-1/2,1/2,m+1),repmat(S23theory,m+1,1),'g')
hold on
plot(linspace(-1/2,1/2,m+1),real(S23),'k')
hold off
legend('S12theory', 'S12','S23theory','S23')
%saveas(gcf, 'bern.2pings.crossspectra.jpg')


%plot the result
figure
plot(linspace(-1/2,1/2,m+1),repmat(R12theory,m+1,1),'b')
hold on
plot(linspace(-1/2,1/2,m+1),R12_est,'r')
hold on
plot(linspace(-1/2,1/2,m+1),repmat(R13theory,m+1,1),'g')
hold on
plot(linspace(-1/2,1/2,m+1),R13_est,'k')
hold off
legend('R12theory', 'R12','R13theory','R13')
%ylim([0 1])
%saveas(gcf, 'bern.2pings.coherence.jpg')

%half frequency plot, as symmetric
% figure
% plot(linspace(0,1/2,m/2 + 1),repmat(S11theory2,m/2+1,1),'b')
% hold on
% plot(linspace(0,1/2,m/2 + 1),S11((m/2 + 1):end),'r')
% hold on
% plot(linspace(0,1/2,m/2 + 1),repmat(S22theory2,m/2+1,1),'g')
% hold on
% plot(linspace(0,1/2,m/2 + 1),S22((m/2 + 1):end),'k')
% hold off
% legend('S11theory', 'S11','S22theory','S22')



figure
plot(linspace(0,1/2,m+1),R_partial(:,1),'b')
hold on
plot(linspace(0,1/2,m+1),R_partial(:,2),'r')
hold on
plot(linspace(0,1/2,m+1),R_partial(:,3),'g')
hold on
plot(linspace(0,1/2,m+1),repmat(bcutoff,1,m+1),'k')
hold on
plot(linspace(0,1/2,m+1),repmat(critR12,1,m+1),'--b')
hold on
plot(linspace(0,1/2,m+1),repmat(critR13,1,m+1),'--r')
hold on
plot(linspace(0,1/2,m+1),repmat(critR23,1,m+1),'--g')
hold on
%plot(linspace(0,1/2,m+1),repmat(critC,1,m+1),'--k')
hold off
legend('partialR12','partialR13', 'partialR23', 'sig line', 'crit12', 'crit13', 'crit23' )
%saveas(gcf, 'bern.2pings.hyptesting.jpg')

% figure
% plot(linspace(0,1/2,m+1),R_partialdirect(:,1),'b')
% hold on
% plot(linspace(0,1/2,m+1),R_partialdirect(:,2),'r')
% hold on
% plot(linspace(0,1/2,m+1),R_partialdirect(:,3),'g')
% hold on
% plot(linspace(0,1/2,m+1),repmat(bcutoff,1,m+1),'k')
% hold on
% plot(linspace(0,1/2,m+1),repmat(critR,1,m+1),'--k')
% hold off
% legend('directpartialR12','directpartialR13', 'directpartialR23', 'sig line', 'multiple hyp line')


%multitapering function for p=3 (3 times series)
function [spectrum_ests, cross_ests, spectral_mat] = spectrum_est(X_1,X_2,X_3,number_seq)
%N needs to be equal to length of sequence - sequences should be of same
%size when input 
N = length(X_1);
%Select a bandwidth, usually NW is 2, 5/2, 3, 7/2, or 4
time_halfbandwidth = number_seq/2 + 1;
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

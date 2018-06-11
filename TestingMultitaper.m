%here we test the multitapering function by doing a goodness of fit with
%the goodman distribution
%we now can generate random poisson processes

% arrival rate per second
lambdaAB=1/10;  %X1 rate
lambdaBC=1/40;  %X2 rate
lambdaCA=1/30; %X3 rate

%make an array of lambda values

%lambdas = [lambdaAB; lambdaBC];
lambdas = [lambdaAB; lambdaBC; lambdaCA];
N = length(lambdas);

%simulation time in seconds (1 hour)
T = 10*3600; 

%we will bin the counts to 0:1/m:1, ie rate per unit time is lambda*T/m
m = 10000 - 1;
%n gives the number of iid estimates
n=1000;
f=200;
results = zeros(n,6);
ntapers = 10;
for i=1:n
    AB = rand(1,poissrnd(lambdaAB*T));
    BC = rand(1,poissrnd(lambdaBC*T));
    CA = rand(1,poissrnd(lambdaCA*T));
    ABfreq = transpose(hist(AB,0:(1/m):1));
    BCfreq = transpose(hist(BC,0:(1/m):1));
    CAfreq = transpose(hist(CA,0:(1/m):1));
    freqs = [ABfreq, BCfreq, CAfreq];
    newfreqs = freqs;
    %initialise counts
    counts = zeros(m+1,N);

    %counts contains the zero mean frequencies 
    for j=1:N
        counts(:,j) = newfreqs(:,j) - mean(newfreqs(:,j));
    end
    
    [X_ests, X_cross_ests, X_matrix] = spectrum_est(counts(:,1),counts(:,2),counts(:,3),ntapers);
    results(i,1) = X_ests(f,1);
    results(i,2) = X_ests(f,2);
    results(i,3) = X_ests(f,3);
    results(i,4) = X_cross_ests(f,1);
    results(i,5) = X_cross_ests(f,2);
    results(i,6) = X_cross_ests(f,3);
end


%for this plot, best to use 40 tapers, and we only need n=1
% figure
% plot(-(1/2) + (1/(m+1))*(0:m), X_ests(:,1), 'r')
% hold on
% plot(-(1/2) + (1/(m+1))*(0:m), repmat(lambdaAB*T/m, m+1,1), '-k')
% hold on
% plot(-(1/2) + (1/(m+1))*(0:m), X_ests(:,2), 'b')
% hold on
% plot(-(1/2) + (1/(m+1))*(0:m), repmat(lambdaBC*T/m, m+1,1),'--k')
% hold on
% plot(-(1/2) + (1/(m+1))*(0:m), X_ests(:,3), 'g')
% hold on
% plot(-(1/2) + (1/(m+1))*(0:m), repmat(lambdaCA*T/m, m+1,1),':k')
% legend('S11', 'S11theory', 'S22', 'S22theory', 'S33', 'S33theory')



R12_est = sort((abs(results(:,4).^2))./(results(:,1).*results(:,2)));
R13_est = sort((abs(results(:,5).^2))./(results(:,1).*results(:,3)));
R23_est = sort((abs(results(:,6).^2))./(results(:,2).*results(:,3)));

goodman12 = Goodman_QQ_Plots(0,R12_est,ntapers);
goodman13 = Goodman_QQ_Plots(0,R13_est,ntapers);
goodman23 = Goodman_QQ_Plots(0,R23_est,ntapers);

figure
plot(R12_est,goodman12)
hold on
plot(R13_est,goodman13)
hold on
plot(R23_est,goodman23)
hold on
line([0 1],[0 1])
hold off
legend('\gamma^2_{12}','\gamma^2_{13}','\gamma^2_{23}')
title('A QQ Plot: Coherence Estimates vs Goodman Distribution')

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

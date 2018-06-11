%we now can generate random poisson processes
m = 100000;

% arrival rate per second
lambdaAB=1/50;  %X1 rate
lambdaBC=1/30;  %X2 rate
lambdaCA=1/60; %X3 rate

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
%f gives frequency so that delta = f/m
f1 = 10;
delta1 = f1/m;
%set a random start time t_0 between 0 and delta
t0a = rand(1)*delta1;
    
times1 = transpose(t0a:(1/m * f1):1);

%we can now bin these times

%we will bin the counts to 0:1/m:1, ie rate per unit time is lambda*T/m
lambdas = lambdas*T/m;

%this generates independent poisson processes, X_1, X_2, X_3
X1times = [AB;times1];
X2times = [BC;times1];
X3times = [CA;times1];

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

%calculate coherences
R12_est = (abs(S12.^2))./(S11.*S22);
R13_est = (abs(S13.^2))./(S11.*S33);
R23_est = (abs(S23.^2))./(S22.*S33);

%generate theoretical values in line with those derived
lambdas2 = lambdas .* (1-lambdas) + (delta1 * (1-delta1));
S11baseline = repmat(lambdas2(1),m+1,1);
S11theory = S11baseline;
S22baseline = repmat(lambdas2(2),m+1,1);
S22theory = S22baseline;
S33baseline = repmat(lambdas2(3),m+1,1);
S33theory = S33baseline;
for i = 1:(m*delta1)
    S11theory(i/delta1) = pi*f1;
    S22theory(i/delta1) = pi*f1;
    S33theory(i/delta1) = pi*f1;
end
S11theory(m/2) = lambdas2(1);
S22theory(m/2) = lambdas2(2);
S33theory(m/2) = lambdas2(3);

%generate theoretical cross spectra
S12theory = repmat((delta1 * (1-delta1)), m+1,1);
S13theory = repmat((delta1 * (1-delta1)), m+1,1);
S23theory = repmat((delta1 * (1-delta1)), m+1,1);
for i = 1:(m*delta1)
    S12theory(i/delta1) = pi*f1;
    S13theory(i/delta1) = pi*f1;
    S23theory(i/delta1) = pi*f1;
end
S12theory(m/2) = (delta1 * (1-delta1)); 
S13theory(m/2) = (delta1 * (1-delta1));
S23theory(m/2) = (delta1 * (1-delta1));

%compute the coherence
R12theory = (S12theory.^2)./(S11theory.*S22theory);
R13theory = (S13theory.^2)./(S11theory.*S33theory);
R23theory = (S23theory.^2)./(S22theory.*S33theory);

%plot the spectra
figure
plot(linspace(-1/2,1/2,m+1),S11,'r')
hold on
plot(linspace(-1/2,1/2,m+1),S22,'g')
hold on
plot(linspace(-1/2,1/2,m+1),S11theory,'--b')
hold on
plot(linspace(-1/2,1/2,m+1),S22theory,'--k')
hold off
legend('S11', 'S22','S11theory','S22theory')
%saveas(gcf, 'Regpulse_spectra.jpg')

%plot the cross-spectra
figure
plot(linspace(-1/2,1/2,m+1),S12theory,'-k')
hold on
plot(linspace(-1/2,1/2,m+1),real(S12),'g')
hold off
legend('S12theory', 'S12')
%saveas(gcf, 'Regpulse_crossspectra.jpg')

% we could plot these with half frequency, as the functions are symmetric
% figure
% plot(linspace(0,1/2,m/2 + 1),S11theory((m/2 + 1):end),'b')
% hold on
% plot(linspace(0,1/2,m/2 + 1),S11((m/2 + 1):end),'r')
% hold on
% plot(linspace(0,1/2,m/2 + 1),S22theory((m/2 + 1):end),'g')
% hold on
% plot(linspace(0,1/2,m/2 + 1),S22((m/2 + 1):end),'k')
% hold off
% legend('S11theory', 'S11','S22theory','S22')

%plots to produce comparison of coherence estimate
figure
plot(linspace(-1/2,1/2,m+1),R12_est,'b')
hold on
plot(linspace(-1/2,1/2,m+1),R13_est,'r')
hold on
plot(linspace(-1/2,1/2,m+1),R12theory,'k')
hold on
plot(linspace(-1/2,1/2,m+1),R13theory,'--k')
hold off
legend('R12','R13', 'R12theory', 'R13theory')
%saveas(gcf, 'Regpulsecoherence.jpg')


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
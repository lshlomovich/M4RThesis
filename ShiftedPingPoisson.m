%we now can generate random poisson processes

% arrival rate per second
lambdaAB=1/75;  %X1 rate
lambdaBC=1/75;  %X2 rate
lambdaCA=1/75; %X3 rate

%make an array of lambda values

%lambdas = [lambdaAB; lambdaBC];
lambdas = [lambdaAB; lambdaBC; lambdaCA];
N = length(lambdas);

%simulation time in seconds (1 hour)
T = 10*3600; 

%times of events on (0,1) open interval 
AB = rand(1,poissrnd(lambdaAB*T));
BC = rand(1,poissrnd(lambdaBC*T));
CA = rand(1,poissrnd(lambdaCA*T));

%we will bin the counts to 0:1/m:1, ie rate per unit time is lambda*T/m
m = 100000;

%this generates independent poisson processes, X_1, X_2, X_3
ABfreq = transpose(hist(AB,0:(1/m):1));
BCfreq = transpose(hist(BC,0:(1/m):1));
CAfreq = transpose(hist(CA,0:(1/m):1));

%freqs = [ABfreq, BCfreq];
freqs = [ABfreq, BCfreq, CAfreq];

%add a ping - regular pulse
%f gives frequency so that delta = f/m
f1 = 10;
delta1 = f1/m;
%set a random start time t_0 between 0 and delta
t0a = rand(1)*delta1;

times1 = transpose(t0a:(1/m * 10):1);

%jitter the ping times 
epsilon = 0.0005;
shifttimes1 = times1 + epsilon;
shifttimes1 = shifttimes1(shifttimes1<=1);

%bin these times
ping1 = transpose(hist(times1,0:(1/m):1));
sping1 = transpose(hist(shifttimes1,0:(1/m):1));

%create matrix of pings to be added to each process
newfreqs = freqs + [ping1, sping1, ping1];

%initialise counts
counts = zeros(m+1,N);

%counts contains the zero mean frequencies 
for i=1:N
    counts(:,i) = newfreqs(:,i) - mean(newfreqs(:,i));
end

%array of estimated spectral matrices in X_matrix
%X_ests contains the spectrum estimates, and X_cross_ests contains the
%cross spectral estimates
[X_ests, X_cross_ests, X_matrix] = spectrum_est(counts(:,1),counts(:,2),counts(:,3),40);

S11 = X_ests(:,1);
S22 = X_ests(:,2);
S33 = X_ests(:,3);
S12 = X_cross_ests(:,1);
S13 = X_cross_ests(:,2);
S23 = X_cross_ests(:,3);

%S12 = real(X_cross_ests(:,1));
%S13 = real(X_cross_ests(:,2));
%S23 = real(X_cross_ests(:,3));

%R_partial returns complex if this line commented
X_matrix = real(X_matrix);

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
    %spectral_m = X_matrix((3*i+1):3*(i+1),:);
    
    %invert and store
    inverted = inv(X_matrix(:,:,i));
    %use values to compute partial coherence
    
    R_partial(i,1) = (abs(inverted(1,2).^2))./(inverted(1,1).*inverted(2,2));
    R_partial(i,2) = (abs(inverted(1,3).^2))./(inverted(1,1).*inverted(3,3));
    R_partial(i,3) = (abs(inverted(2,3).^2))./(inverted(2,2).*inverted(3,3));
end

%R_partial_12 = S12 - (S13.*(1./S33).*S23);

%plot the results
figure
plot(linspace(0,1/2,m+1),real(S12),'r')
hold on
plot(linspace(0,1/2,m+1),real(S13),'--b')
hold off
legend('S12', 'S13')

figure
plot(linspace(0,1/2,m+1),S11,'r')
hold on
plot(linspace(0,1/2,m+1),S22,'--b')
hold off
legend('S11', 'S22')

%by plotting this we show that the coherence between 1 and 2 is the same as
%between 1 and 3
figure
plot(linspace(0,1/2,m+1),R12_est,'k')
hold on
plot(linspace(0,1/2,m+1),R13_est,'--r')
hold off
legend('R12', 'R13')
%saveas(gcf, 'shiftpingcoherence.jpeg')

figure
plot(linspace(0,1/2,m+1),R_partial(:,1),'k')
hold on
plot(linspace(0,1/2,m+1),R_partial(:,2),'--r')
hold off
legend('partialR12', 'partialR13')
%saveas(gcf, 'shiftpingpartial.jpeg')


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
spectral_mat(2,1,:) = FS12;

spectral_mat(1,3,:) = FS13;
spectral_mat(3,1,:) = FS13;

spectral_mat(2,3,:) = FS23;
spectral_mat(3,2,:) = FS23;
    
end
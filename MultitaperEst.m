%Generate known processes

%Simulate MA(q)
q = 2;
%simulation time/steps
S = 1024*10;
x = ones(S,1);
x(1) = 1;
x(2) = 0.5;
theta1 = 0.2;
theta2 = 0.2;
sigma = 1;
mu = 0;
epsx = normrnd(mu, sigma, S, 1);

%Start the loop running from obs. q+1 to S
for t=(q+1):S    
    %The MA(q) model
    x(t) = epsx(t) - theta1*epsx(t-1) - theta2*epsx(t-2);    
end

%Calculate true sdf of MA(q) from literature
fvaluesx = transpose(linspace(-0.5, 0.5, S+1));
sdfx = sigma*(abs((1 - theta1*exp(-1*sqrt(-1)*2*pi*fvaluesx) - theta1*exp(-1*sqrt(-1)*4*pi*fvaluesx)).^2));


%Simulate AR(p)
p = 2;
T = 1024*10;                        %Number of observations
y = ones(T,1);                      %Vector stores simulations
y(1) = 1;                           %Set the first obs. to 1    
y(2) = 0.5;                         %Set the second obs. to 0.5
phi1 = 0.2;                         %Set the value of phi1 (coefficient on y(t-1))
phi2 = 0.2;                         %Set the value of phi2 (coefficient on y(t-2))
sigma_e = 1;                        %Set s.d. of the error term
mu_e = 0;                           %Set mean of the error term
eps = normrnd(mu_e, sigma_e, T, 1); %Creat a vector of normal random numbers with mean, mu_e and s.d. sigma. Dimension Tx1 

%Start the loop running from obs. p+1 to T
for t=(p+1):T
    %The AR(p) model
    y(t) = phi1*y(t-1) + phi2*y(t-2) + eps(t);    
end

%Calculate true sdf of AR(p) from literature 
fvaluesy = transpose(linspace(-0.5, 0.5, T+1));
sdfy = sigma./(abs(1 - phi1*exp(-1*sqrt(-1)*2*pi*fvaluesy) - phi1*exp(-1*sqrt(-1)*4*pi*fvaluesy)).^2);



%Multitapering estimates using fft and fftshift

SAB = x; %line only when used with 'Generate.m'
SBC = y;

%N needs to be equal to length of sequence
NSAB = length(SAB);
NSBC = length(SBC);
%number of sequences 
num_seq = 40;
%Select a bandwidth, usually NW is 2, 5/2, 3, 7/2, or 4
time_halfbandwidth = num_seq/2 + 1;
%Obtain DPSSs
H = dpss(NSAB, time_halfbandwidth, num_seq);
%Compute fourier transform
FSAB = zeros(NSAB, num_seq);
FSBC = zeros(NSBC, num_seq);
%FSAC = zeros(NSAB, num_seq);
for i=1:num_seq
    FSAB(:,i) = abs(fftshift(fft(H(:,i).*SAB))).^2;
    FSBC(:,i) = abs(fftshift(fft(H(:,i).*SBC))).^2;
    %FSAC(:,i) = fftshift(fft(H(:,i).*SAB)).*conj(fftshift(fft(H(:,i).*SBC)));
end
FSAB = mean(FSAB, 2);
FSBC = mean(FSBC, 2);
%FSAC = mean(FSAC, 2);
%Rhat = abs((FSAC).^2)./(FSAB.*FSBC);

figure
plot(-(1/2) + (1/NSAB)*(0:NSAB-1), FSAB, 'r')
hold on
plot(fvaluesx,sdfx,'b')
xlim([0 0.5])
title('Spectral Estimation for MA(2)')
legend('Multitaper Estimate','True SDF')
saveas(gcf, 'MAestimate.jpg')

figure
plot(-(1/2) + (1/NSBC)*(0:NSBC-1), FSBC, 'r')
hold on
plot(fvaluesy,sdfy,'b')
xlim([0 0.5])
title('Spectral Estimation for AR(2)')
legend('Multitaper Estimate','True SDF')
saveas(gcf, 'ARestimate.jpg')

% figure
% title('Spectral Estimation for L=40')
% hold on
% plot(-(1/2) + (1/NSBC)*(0:NSBC-1), FSBC, 'b')
% hold on
% plot(-(1/2) + (1/NSAB)*(0:NSAB-1), Rhat, 'k')
% legend('SAB', 'SBC', 'R')
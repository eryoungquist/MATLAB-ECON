function GBM = GeometricBrownianMotion(S0,mu,sigma,T,N)
%EM Euler-Maruyama method on linear SDE
%
% SDE is dX = lambda*X dt + mu*X dW, X(0) = Xzero,
% where lambda = 2, mu = 1 and Xzero = 1.
%
% Discretized Brownian path over [0,1] has dt = 2^(-8).
% Euler-Maruyama uses timestep R*dt.
randn('state',100)
lambda = 2; mu = 1; Xzero = 1; % problem parameters
T = 1; N = 2^8; dt = 1/N;
dW = sqrt(dt)*randn(1,N); % Brownian increments
%W = cumsum(dW); % discretized Brownian path
%Xtrue = Xzero*exp((lambda-0.5*mu^2)*([dt:dt:T])+mu*W);
%plot([0:dt:T],[Xzero,Xtrue],'m-'), hold on
R = 4; Dt = R*dt; L = N/R; % L EM steps of size Dt = R*dt
Xem = zeros(1,L); % preallocate for efficiency
Xtemp = Xzero;

r=0;
pz = r*Dt;
Z = rand(1,L) > pz;
for j = 1:L
    Winc = sum(dW(R*(j-1)+1:R*j));
    Xtemp = Xtemp + Dt*lambda*Xtemp + mu*Xtemp*Winc;
    %For some probability
    % if condition is met Xtemp = Xzero;
    if Z(j) == 0
        Xtemp = Xzero;
    end
    Xem(j) = Xtemp;

end

plot([0:Dt:T],[Xzero,Xem],'r--*'), hold off
xlabel('Time','FontSize',12)
ylabel('Position','FontSize',16,'Rotation',0,'HorizontalAlignment','right')
title('Simulated GBM')
%emerr = abs(Xem(end)-Xtrue(end))

clear;
randn('state',100)
 mu = 1; sigma = 1; Xzero = 1; % problem parameters
dt = 0.001;
Xthresh = 2; %threshold for function to stop
Xr = Xzero;
r = 1;
pz = r*dt;

n = 1000;
FPT = zeros(1,n);
for j = 1:n
    clear g;
    Xtemp = Xzero;
    t = 0;
    i = 1;
    %g = zeros(2,16); %preallocate array
    while Xtemp < Xthresh
        Z = pz >= rand;
        %z(i) = Z;
        dW = sqrt(dt)*randn;
        Xtemp = Xtemp + (1-Z)*Xtemp*(mu*dt+sigma*dW)+Z*(Xthresh-Xtemp);

        i = i+1;
        %g(:,i) = [i*dt Xtemp]';
        %if size(g(:,1)) == i
        %   g = [g,zeros(2,i)]; %resize array
        %end
    end
     FPT(j) = i*dt;
 end
%z
%plot(g(1,:),g(2,:),'--*')

MFPT = sum(FPT)/n
%The equations based on the paper
q1= (sqrt((sigma.^2-2.*mu).^2 +(8.*r.*sigma.^2))+(sigma.^2-2.*mu))./(2.^sigma.^2);
Tx0=(Xzero./Xthresh).^q1;
TxR=(Xr./Xthresh).^q1;
Tr = (1-Tx0)./(r.*TxR)



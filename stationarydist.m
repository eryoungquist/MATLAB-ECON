syms x t x0 S m xr ru
S=1;
m=1;
P0= exp((-log(x/x0)-(m-(S^2/2)*t)^2)/(t*2*S^2))/(x*sqrt(2*t*pi*S^2));

Pr = exp(-r*t)* P0 + r*integral(exp(-ru)*P0,0,t);
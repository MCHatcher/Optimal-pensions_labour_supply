//OLG pensions: Dynare replication of the first-best. Last updated: April 16, 2024. 
%Written by M. Hatcher (m.c.hatcher@soton.ac.uk)

var y, cy, co, k, l, P, prod, y1, cy1, co1, k1, l1, MPK, gap, gap0, gap1, gap2, gap3;
//predetermined_variables tau_c;
varexo e; //u;

parameters alpha, betta, chi, omega_bar, theta, n,  rho, rho1, cystar, costar, kstar, lstar, ystar, tau_p, tau_c;

alpha = 0.30;
betta = 0.85;
chi = 2;
theta = 4;
omega_bar = 0.995;
n = 0.05;
rho = 0.95;
rho1 = 0.95;

lstar = ( (1+betta/omega_bar)*(1-alpha)/ (theta*(1-alpha*omega_bar)) )^(1/(1+chi)); 

kstar = (alpha*omega_bar*(1+n)^(-alpha))^(1/(1-alpha))*lstar;

ystar = (kstar/(1+n))^alpha*lstar^(1-alpha);

cystar = omega_bar/(betta + omega_bar)*(1-alpha*omega_bar)*ystar;

costar = betta*(1+n)/omega_bar*cystar;

tau_p = ( betta*(1-alpha) - alpha*omega_bar*(1+betta) )/ ( betta*(1-alpha) - alpha*omega_bar*(1-omega_bar) );

tau_c = -tau_p;

model;

//Shocks
prod = prod(-1)^rho*exp(e);

//Planner
co = (betta/omega_bar)*(1+n)*cy;

y = prod*(k(-1)/(1+n))^alpha*l^(1-alpha);

1/cy = betta*( alpha*y(+1)*(1+n)/k*1/co(+1) );

cy =  y - k - co/(1+n);

theta*l^(1+chi) = (1-alpha)*y/cy;

//Decentralized economy
y1 = prod*(k1(-1)/(1+n))^alpha*l1^(1-alpha);

1/( (1+tau_c)*cy1 ) = betta*MPK(+1) / ( (1+tau_c(+1))*co(+1) );

P = (1+n)*(tau_p*(1-alpha)*y1 + tau_c*(1-alpha*omega_bar)*y1);

theta*l1^(1+chi) = (1-alpha)*(1-tau_p)*y1/((1+tau_c)*cy1);

co1 = (MPK*k1(-1) + P)/(1+tau_c);

cy1 = ( (1-alpha)*(1-tau_p)*y1 - k1 )/(1+tau_c);

MPK = alpha*(1+n)*y1/k1(-1);

//Checks
gap = y - y1;

gap0 = cy - cy1;

gap1 = co - co1;

gap2 = l - l1;

gap3 = 1/( (1+tau_c)*cy1 ) - betta*MPK(+1) / ( (1+tau_c(+1))*co(+1) ) ; 

end;

initval;
y = ystar;
cy = cystar;
co = costar;
l = lstar; 
k = kstar;
y1 = ystar;
cy1 = cystar;
co1 = costar;
l1 = lstar; 
k1 = kstar;
prod = 1;
MPK = (1+n)/omega_bar;
gap = 0;
gap0 = 0;
gap1 = 0;
gap2 = 0;
gap3 = 0;
end;

shocks;
var e; stderr 0.01;
//var u; stderr 0.01;
end;

steady; 
check;

stoch_simul(order=2, periods = 1000, irf=0);




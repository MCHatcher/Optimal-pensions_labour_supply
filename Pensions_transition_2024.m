%Pensions transition 2023. Last updated: April 16, 2024. 
%Written by M. Hatcher (m.c.hatcher@soton.ac.uk)

clear, clc, %close all

%Calibration
alpha = 0.30;
betta = 0.85;
chi = 2;
theta = 4;
omega = 0.995;
n = 0.05;

%Other parameters
T0 = 20; T_pre = 10;
T = 1e5;
T_plot = 30;

%Steady state (planner)
lss = ( (1+betta/omega)*(1-alpha)/ (theta*(1-alpha*omega)) )^(1/(1+chi)); 
kss = (alpha*omega*(1+n)^(-alpha))^(1/(1-alpha))*lss;
yss = (kss/(1+n))^alpha*lss^(1-alpha);
cyss = omega/(betta + omega)*(1-alpha*omega)*yss;
coss = betta*(1+n)/omega*cyss;
betta_crit = alpha*omega/(1-alpha*(1+omega));

%Initial values
k_init = 1.2*kss;
lstar = lss;
A = ones(T+1,1);
tau_c = 0;  %see below
tau_p = 0; %see below
betta_hat = betta*alpha / ( alpha + (1-alpha)*tau_p + (1-alpha*omega)*tau_c );
l = ( (1+betta_hat)/theta )^(1/(1+chi));
l1 = ( (1+betta) / theta )^(1/(1+chi));

%Post transtion
tau_p1 = ( betta*(1-alpha) - alpha*omega*(1+betta) )/ ( betta*(1-alpha) - alpha*omega*(1-omega) );
tau_c1 = -tau_p1;
betta_hat1 = betta*alpha / ( alpha + (1-alpha)*tau_p1 + (1-alpha*omega)*tau_c1 );
l_new = ( (1+betta_hat1)/theta )^(1/(1+chi));

%Initial matrices
ystar = NaN(T,1); labstar = ystar; cystar = ystar; costar = ystar; kstar = ystar; ystar_prime = ystar; 
costar_prime = ystar; Ustar = ystar; Vstar = Ustar; y = ystar; cy = ystar; co = ystar; k = ystar; 
y_prime = ystar; co_prime = ystar; U = ystar; V = U; lambda = U; lab = U;
y1 = U; cy1 = U; co1 = U; k1 = U; y1_prime = U; co1_prime = U; U1 = U; V1 = U; lab1 = U;

%-------------------
%Planner solution
%-------------------

ystar(1) = A(1)*(k_init/(1+n))^alpha*lstar^(1-alpha);

labstar(1) = lstar;

cystar(1) = omega/(betta+omega)*(1-alpha*omega)*ystar(1);

costar(1) = betta/omega*(1+n)*cystar(1);

U_init_star = betta*log(costar(1));

kstar(1) = alpha*omega*ystar(1);

ystar_prime(1) = A(2)*(kstar(1)/(1+n))^alpha*lstar^(1-alpha);

costar_prime(1) = betta*(1+n)/(betta+omega)*(1-alpha*omega)*ystar_prime(1);

Ustar(1) = log(cystar(1)) + betta*log(costar_prime(1)) - (theta/(1+chi))*lstar^(1+chi);

Vstar(1) = Ustar(1); 

%-----------------------
%Decentralized economy
%-----------------------

y(1) = A(1)*(k_init/(1+n))^alpha*l^(1-alpha);

lab(1) = l;

k(1) = betta_hat*(1-alpha)*(1-tau_p)/(1+betta_hat)*y(1);

cy(1) = k(1)/( betta_hat*(1+tau_c) );

co(1) = (alpha + (1-alpha)*tau_p + (1-alpha*omega)*tau_c )*(1+n)/(1+tau_c)*y(1);

U_init = betta*log(co(1));

y_prime(1) = A(2)*(k(1)/(1+n))^alpha*l^(1-alpha);

co_prime(1) = (alpha + (1-alpha)*tau_p + (1-alpha*omega)*tau_c )*(1+n)/(1+tau_c)*y_prime(1);

U(1) = log(cy(1)) + betta*log(co_prime(1)) - (theta/(1+chi))*l^(1+chi);

V(1) = U(1);

%-----------------------
%Laissez faire
%-----------------------

y1(1) = A(1)*(k_init/(1+n))^alpha*l1^(1-alpha);

lab1(1) = l1;

k1(1) = betta/(1+betta)*(1-alpha)*y1(1);

cy1(1) = k1(1)/betta;

co1(1) = alpha*(1+n)*y1(1);

U1_init = betta*log(co1(1));

y1_prime(1) = A(2)*(k1(1)/(1+n))^alpha*l1^(1-alpha);

co1_prime(1) = alpha*(1+n)*y1_prime(1);

U1(1) = log(cy1(1)) + betta*log(co1_prime(1)) - (theta/(1+chi))*l1^(1+chi);

V1(1) = U1(1);

%-------------------
%Simulation loop
%-------------------

    for t=2:T

        %----------------------
        %Dcentralized economy
        %----------------------

        lab(t) = l;

        y(t) = A(t)*(k(t-1)/(1+n))^alpha*l^(1-alpha); 

        k(t) = betta_hat*(1-alpha)*(1-tau_p)/(1+betta_hat)*y(t);

        cy(t) = k(t)/( betta_hat*(1+tau_c) );

        co(t) = (alpha + (1-alpha)*tau_p + (1-alpha*omega)*tau_c)*(1+n)/(1+tau_c)*y(t);

        y_prime(t) = A(t+1)*(k(t)/(1+n))^alpha*l^(1-alpha);

        co_prime(t) = (alpha + (1-alpha)*tau_p + (1-alpha*omega)*tau_c)*(1+n)/(1+tau_c)*y_prime(t);

        U(t) = log(cy(t)) + betta*log(co_prime(t)) - (theta/(1+chi))*l^(1+chi);

        V(t) = omega^(t-1)*U(t);

        if t >= T0

            lab(t) = l_new;

            y(t) = A(t)*(k(t-1)/(1+n))^alpha*l_new^(1-alpha);

            k(t) = betta_hat1*(1-alpha)*(1-tau_p1)/(1+betta_hat1)*y(t);

            cy(t) = k(t)/( betta_hat1*(1+tau_c1) );

            co(t) = (alpha + (1-alpha)*tau_p1 + (1-alpha*omega)*tau_c1)*(1+n)/(1+tau_c1)*y(t);

            y_prime(t) = A(t+1)*(k(t)/(1+n))^alpha*l_new^(1-alpha);

            co_prime(t) = (alpha + (1-alpha)*tau_p1 + (1-alpha*omega)*tau_c1)*(1+n)/(1+tau_c1)*y_prime(t);

            U(t) = log(cy(t)) + betta*log(co_prime(t)) - (theta/(1+chi))*l_new^(1+chi);

            V(t) = omega^(t-1)*U(t);

            if t==T0
                U_init0 = betta*log( (alpha + (1-alpha)*tau_p1 + (1-alpha*omega)*tau_c1)*(1+n)/(1+tau_c1)*y(T0) ); 
            end

        end

        %-----------------------
        %Laissez faire
        %-----------------------

        y1(t) = A(t)*(k1(t-1)/(1+n))^alpha*l1^(1-alpha);

        lab1(t) = l1;

        k1(t) = betta/(1+betta)*(1-alpha)*y1(t);

        cy1(t) = k1(t)/betta;

        co1(t) = alpha*(1+n)*y1(t);

        y1_prime(t) = A(t+1)*(k1(t)/(1+n))^alpha*l1^(1-alpha);

        co1_prime(t) = alpha*(1+n)*y1_prime(t);
        
        U1(t) = log(cy1(t)) + betta*log(co1_prime(t)) - (theta/(1+chi))*l1^(1+chi);

        V1(t) = omega^(t-1)*U1(t);

        if t==T0
            U1_init0 = betta*log( (alpha)*(1+n)*y1(T0) );
        end

        %---------
        %Planner
        %---------

        ystar(t) = A(t)*(kstar(t-1)/(1+n))^alpha*lstar^(1-alpha);

        labstar(t) = lstar;

        cystar(t) = omega/(betta+omega)*(1-alpha*omega)*ystar(t);

        costar(t) = betta/omega*(1+n)*cystar(t);

        kstar(t) = alpha*omega*ystar(t);

        if t==T0-1
            kstar(t) = k(t);
        end

        ystar_prime(t) = A(t+1)*(kstar(t)/(1+n))^alpha*lstar^(1-alpha);

        costar_prime(t) = betta*(1+n)/(betta+omega)*(1-alpha*omega)*ystar_prime(t);

        Ustar(t) = log(cystar(t)) + betta*log(costar_prime(t)) - (theta/(1+chi))*lstar^(1+chi);

        Vstar(t) = omega^(t-1)*Ustar(t);

        if t==T0
            Ustar_init0 = betta*log( betta*(1+n)/(betta+omega)*(1-alpha*omega)*ystar(T0) );
        end

    end
        
%--------------------
%Welfare analysis
%--------------------
lambda_init0 = 100*( exp( (U_init0-U1_init0)/(1+betta) ) -1);
        
    for t=1:T
        
        lambda(t) = 100*( exp( (U(t)-U1(t))/(1+betta) ) -1);
    end

%---------
%Checks
%---------
Wstar = Ustar_init0/omega + sum(Vstar(T0:end)); 
W = U_init0/omega + sum(V(T0:end)); 

gap = Wstar - W
kstar(end) - kss
max(abs(kstar(T0-1:end) - k(T0-1:end)))

%---------
%Plots
%---------
Welfare = [lambda(T0-(T_pre-1):T0-1); lambda_init0; lambda(T0:T_plot)];
Time = [-T_pre:T_plot-T0];
set(0,'DefaultLineLineWidth',1)
lab1_plot = lab1; labstar_plot = labstar; kstar_plot = kstar; cy1_plot = cy1; k1_plot = k1;
lab1_plot(T0-T_pre:T0-1) = NaN; labstar_plot(T0-T_pre:T0-1) = NaN;  kstar_plot(T0-T_pre:T0-1) = NaN; 
%cy1_plot(T0-T_pre:T0-1) = NaN; k1_plot(T0-T_pre:T0-1) = NaN;

figure(1)
subplot(2,3,1), plot(Time,Welfare,'k'), title('Welfare effects'), xlabel('Generations (d.o.b.)'), ylabel('% c.e.') 
subplot(2,3,4), hold on, plot(Time,kstar_plot(T0-T_pre:T_plot)/k1(T0),'k'), plot(Time,k(T0-T_pre:T_plot)/k1(T0),'k'), plot(Time,k1_plot(T0-T_pre:T_plot)/k1(T0),'--k'), 
title('Capital: k_{t+1}'), xlabel('Time, t') 
subplot(2,3,5), hold on, plot(Time,cy(T0-T_pre:T_plot)/cy1(T0),'k'), plot(Time,co(T0-T_pre:T_plot)/co1(T0),'--k'), plot(Time,cy1_plot(T0-T_pre:T_plot)/cy1(T0),'--k'), 
title('Consumptions'), xlabel('Time, t') 
subplot(2,3,6), hold on, plot(Time,labstar_plot(T0-T_pre:T_plot)/l1,'k'), plot(Time,lab(T0-T_pre:T_plot)/l1,'k'), plot(Time,lab1_plot(T0-T_pre:T_plot)/l1,'--k'), 
title('Labour: l_t'), xlabel('Time, t'), axis([-inf,inf,-inf,inf]) 

figure(2)
subplot(2,3,1), plot(Time,Welfare,'k'), title('Welfare effects'), xlabel('Generations (d.o.b.)'), ylabel('% c.e.'), hold on,
subplot(2,3,2), hold on, plot(Time,k(T0-T_pre:T_plot)/kstar(end),'k'), %plot(Time,kstar_plot(T0-T_pre:T_plot)/kstar(end),'k')
title('Capital: k_{t+1}'), xlabel('Time, t'), hold on
subplot(2,3,4), hold on, plot(Time,cy(T0-T_pre:T_plot)/cystar(end),'k'), title('Consumption (young)'), xlabel('Time, t') 
subplot(2,3,5), hold on, plot(Time,co(T0-T_pre:T_plot)/costar(end),'k'), title('Consumption (old)'), xlabel('Time, t') 
subplot(2,3,6), hold on, plot(Time,labstar_plot(T0-T_pre:T_plot)/lstar,'k'), plot(Time,lab(T0-T_pre:T_plot)/lstar,'k'), 
title('Labour: l_t'), xlabel('Time, t'), %axis([-inf,inf,-inf,inf]) 


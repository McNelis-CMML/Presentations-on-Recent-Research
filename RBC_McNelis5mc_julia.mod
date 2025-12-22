

//==============================================================
// Small open economy with CES production, Q-theory
// investment, and current account dynamics, oil and TVP shocks
// Calvo pricing, exchange rate dynamics and Taylor rule
//
// This version does NOT automatically run perfect_foresight
// Use with Julia script that loads shocks manually
//==============================================================


var
    c       // consumption  (Z)
    k       // capital
    y       // output
    oil     // imported oil input
    i       // investment
    q       // Tobin's q
    yhat    // CES composite (pre-TFP)
    mpk     // marginal product of capital
    mpoil   // marginal product of oil
    l       // labor
    w       // real wage
    Ck      // derivative of adj. cost wrt K_{t-1}
    r       // net interest rate on foreign assets
    nfa     // net foreign assets (level)
    a       // TFP shock (log deviation, around zero)
    a_oil
    lambda
    ymarg
    mc
    pi
    rdom
    s
    p
    rk;


    


varexo eps_a eps_oil;  // TFP shock

parameters delta Phi  alphaK alphaL alphaO kbar lbar obar rhoc ymargbar yhatbar Abar mpkbar mplbar mpobar ibar tauc tauw tauk cbar sigma_c phiL r_ss PhiB rho_a
lambdabar gamma Poilbar sigma_mc mc_ss theta pow_uc rho_oil y_ss rho_r phi_pi phi_y beta sigma_a sigma_oil r_lb;
    

delta   = 0.04;          
Phi     = 10;             
alphaK  = 0.35;
alphaL  = 0.60;
alphaO  = 0.05;
kbar = 3.71;
lbar = .33;
obar = 0.05;
rhoc    = -0.1;         
ymargbar =(1/rhoc)* (alphaK*kbar^rhoc+alphaL*lbar^rhoc+alphaO *obar ^rhoc)^(1/rhoc-1);
yhatbar =(alphaK *kbar^rhoc+alphaL*lbar^rhoc+alphaO *obar^rhoc)^(1/rhoc);
Abar = 1/yhatbar;
ymargbar = alphaK*kbar^rhoc + alphaL*lbar^rhoc + alphaO*obar^rhoc;
yhatbar = ymargbar^(1/rhoc);
y_ss = 1;
mpkbar = Abar *ymargbar^(1/rhoc - 1) * alphaK * kbar^(rhoc - 1);
mplbar = Abar * ymargbar^(1/rhoc - 1) * alphaL * lbar^(rhoc - 1);
mpobar = Abar * ymargbar^(1/rhoc - 1) * alphaO * obar^(rhoc - 1);

Poilbar = mpobar;
ibar = delta * kbar;
tauc	 = 0.1;
tauw	= 0.25;
tauk	 = 0.15;
cbar  = (1 -obar *mpobar - delta*kbar -tauw * mplbar *lbar - tauk * kbar * mpkbar)/(1+tauc);
sigma_c = 2;            
phiL = .5; 
    
PhiB    = 0.001;        
rho_a   = 0.9;    
rho_oil = .9;       
lambdabar =cbar^(-sigma_c)/(1+tauc);
gamma  =mplbar *(1-tauw) *lambdabar/(lbar^(1/phiL));
beta = inv(mpkbar + (1-delta));
r_ss = 1/beta - 1;
sigma_mc = 1/(1-rhoc);
pow_uc = 1/(1-sigma_mc);    


mc_ss = (1/Abar) * ( alphaL^sigma_mc * mplbar^(1-sigma_mc) + alphaK^sigma_mc * mpkbar^(1-sigma_mc)+ alphaO^sigma_mc * mpobar^(1-sigma_mc) )^pow_uc;
theta = .75;
phi_pi = 1.5;
phi_y = .5;
rho_r = .9;
sigma_a = .005;
sigma_oil = .005;
r_lb = 0;


model;
a = rho_a * a(-1) + sigma_a *eps_a;
a_oil = rho_oil * a_oil(-1) + sigma_oil * eps_oil;
yhat = ymarg^(1/rhoc);
ymarg =   alphaK*k(-1)^rhoc + alphaL*l^rhoc + alphaO*oil^rhoc;
y = Abar * exp(a) * yhat;
mpk   = Abar * exp(a) * (ymarg^(1/rhoc - 1)) * alphaK * k(-1)^(rhoc - 1);
mpoil = Abar * exp(a) * (ymarg^(1/rhoc - 1)) * alphaO * oil^(rhoc - 1);
w     = Abar * exp(a) * (ymarg^(1/rhoc - 1)) * alphaL * l^(rhoc - 1);
s/p * Poilbar * exp(a_oil) = mpoil;
k = (1 - delta)*k(-1) + i;
i = delta * k(-1) + (q - 1) * (k(-1)/Phi);
Ck = (Phi/2) * ((i - delta*k(-1))/k(-1))^2 - Phi * ((i - delta*k(-1))/k(-1)) * (i/k(-1));
r = r_ss - PhiB * nfa(-1);
// [mcp = 'rdom > 0']
rdom = (1-rho_r) * r_ss + rho_r * rdom(-1) + (1-rho_r) * phi_pi * pi(+1) + (1-rho_r) * phi_y * log(y(+1)/y_ss);
lambda = c^(-sigma_c)/(1+tauc); 
lambda = beta * (1 + rdom) * lambda(+1);
q = beta * ( mpk + (1 - delta)*q(+1) - Ck);
rk = mc * mpk;
l = (w(-1) *(1-tauw)*lambda/gamma)^phiL; 
nfa = (1 + r) * nfa(-1) + (y - (1+tauc) * c - i - (s/p)*Poilbar*exp(a_oil)* oil - tauw * w * l - tauk * (rk) * k(-1));
mc = (1/Abar * exp(a)) * ( alphaL^sigma_mc * w^(1-sigma_mc) + alphaK^sigma_mc * mpk^(1-sigma_mc)+ alphaO^sigma_mc *(exp(a_oil)*Poilbar)^(1-sigma_mc))^pow_uc;
pi = beta * pi(+1) + ((1-theta)*(1-beta*theta)/theta) * log(mc(-1)/mc_ss);
rdom - r = s(+1) - s;
p = p(-1) + pi;
end;



steady_state_model;
    a    = 0;
    k    = kbar;
    l    = lbar;
    oil  = obar;
    yhat = yhatbar;
    y = 1; 
    mpk   = mpkbar;
    mpoil = mpobar;
    w     = mplbar;
    nfa = 0;
    r   = r_ss;
    q   = 1;
    i   = delta * k;
    Ck  = 0;
    c   = cbar;
    lambda = lambdabar;
    ymarg = ymargbar;
    mc = mc_ss;
    rdom = r_ss;
    pi = 0;
    s = 1; 
    p = 1;
    rk = mpk * mc;
    a_oil = 0;
end;

initval;
 a    = 0;
    k    = kbar;
    l    = lbar;
    oil  = obar;
    yhat = yhatbar;
    y = 1; 
    mpk   = mpkbar;
    mpoil = mpobar;
    w     = mplbar;
    nfa = 0;
    r   = r_ss;
    q   = 1;
    i   = delta * k;
    Ck  = 0;
    c   = cbar;
    lambda = lambdabar;
    ymarg = ymargbar;
    mc = mc_ss;
    rdom = r_ss;
    pi = 0;
    s = 1; 
    p = 1;
    rk = mpk * mc;
    a_oil = 0;

end;

endval;
 a    = 0;
    k    = kbar;
    l    = lbar;
    oil  = obar;
    yhat = yhatbar;
    y = 1; 
    mpk   = mpkbar;
    mpoil = mpobar;
    w     = mplbar;
    nfa = 0;
    r   = r_ss;
    q   = 1;
    i   = delta * k;
    Ck  = 0;
    c   = cbar;
    lambda = lambdabar;
    ymarg = ymargbar;
    mc = mc_ss;
    rdom = r_ss;
    pi = 0;
    s = 1; 
    p = 1;
    rk = mpk * mc;
    a_oil = 0;

end;

steady;
check;

shocks;
  var eps_a;
    stderr 1;

  var eps_oil;
    stderr 1;
end;

stoch_simul(order=1, irf=40, periods = 1000);





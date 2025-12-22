//==============================================================
// Real small open economy with CES production, Q-theory
// investment, and current account dynamics
// (based on RBC1_OpenEconomy_TVPShock_Closure_sigma1.xlsx)
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
    a      // TFP shock (log deviation, around zero)
    lambda
    ymarg;
    


varexo
    eps_a;  // TFP shock

parameters delta Phi beta alphaK alphaL alphaO kbar lbar obar rhoc ymargbar yhatbar Abar mpkbar mplbar mpobar ibar tauc tauw tauk cbar sigma_c phiL r_ss PhiB rho_a
lambdabar gamma Poilbar;
    

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
lambdabar =cbar^(-sigma_c)/(1+tauc);
gamma  =mplbar *(1-tauw) *lambdabar/(lbar^(1/phiL));
beta = inv(mpkbar + (1-delta));
r_ss = 1/beta - 1;

model;
a = rho_a * a(-1) + eps_a;
yhat = ymarg^(1/rhoc);
ymarg =   alphaK*k(-1)^rhoc + alphaL*l^rhoc + alphaO*oil^rhoc;
y = Abar * exp(a) * yhat;
mpk   = Abar * exp(a) * (ymarg^(1/rhoc - 1)) * alphaK * k(-1)^(rhoc - 1);
mpoil = Abar * exp(a) * (ymarg^(1/rhoc - 1)) * alphaO * oil^(rhoc - 1);
w     = Abar * exp(a) * (ymarg^(1/rhoc - 1)) * alphaL * l^(rhoc - 1);
Poilbar = mpoil;
k = (1 - delta)*k(-1) + i;
i = delta * k(-1) + (q - 1) * (k(-1)/Phi);
Ck = (Phi/2) * ((i - delta*k(-1))/k(-1))^2 - Phi * ((i - delta*k(-1))/k(-1)) * (i/k(-1));
r = r_ss - PhiB * nfa(-1);
lambda = c^(-sigma_c)/(1+tauc); 
lambda = beta * (1 + r) * lambda(1);
q = beta * ((1 - delta)*q(1) + mpk - Ck);
l = (w(-1) *(1-tauw)*lambda/gamma)^phiL; 
nfa = (1 + r) * nfa(-1) + y - (1+tauc) * c - i - Poilbar*oil - tauw * w * l - tauk * (mpk) * k(-1);
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

end;

//--------------------
// Shocks
//--------------------

shocks;
    var eps_a; stderr 0.01;
end;

//--------------------
// Solve & simulate
//--------------------

steady;
check;

stoch_simul(order=1, pruning, irf=40, periods = 1000)
    y c i k q oil nfa r;
//send_endogenous_variables_to_workspace;
//send_exogenous_variables_to_workspace;
//send_irfs_to_workspace;

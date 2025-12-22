//==============================================================
// SOE with GK + SIMPLIFIED CCyB 
// CCyB via time-varying capital requirement: phi1_t = phi1_ss + phi2 * gap
// Bank balance sheet: q*k = N + D + KB, where KB = phi1 * q*k
//==============================================================

var
    c k y oil i q yhat mpk mpoil l w Ck r nfa a lambda ymarg
    D N Rk spr nu leverage phi1 KB credit_gap;

varexo
    eps_a;

parameters 
    delta Phi beta alphaK alphaL alphaO kbar lbar obar rhoc 
    ymargbar yhatbar Abar mpkbar mplbar mpobar ibar 
    tauc tauw tauk cbar sigma_c phiL r_ss PhiB rho_a
    lambdabar gamma Poilbar
    thetaN theta nu_ss Nbar Nss Dss KBss Rk_ss
    phi1_ss phi2 credit_ss;

//==============================================================
// CALIBRATION
//==============================================================

delta   = 0.04;          
Phi     = 10;             
alphaK  = 0.35;
alphaL  = 0.60;
alphaO  = 0.05;
kbar    = 3.71;
lbar    = 0.33;
obar    = 0.05;
rhoc    = -0.1;         

ymargbar = alphaK*kbar^rhoc + alphaL*lbar^rhoc + alphaO*obar^rhoc;
yhatbar  = ymargbar^(1/rhoc);
Abar     = 1/yhatbar;

mpkbar  = Abar * ymargbar^(1/rhoc - 1) * alphaK * kbar^(rhoc - 1);
mplbar  = Abar * ymargbar^(1/rhoc - 1) * alphaL * lbar^(rhoc - 1);
mpobar  = Abar * ymargbar^(1/rhoc - 1) * alphaO * obar^(rhoc - 1);

Poilbar   = mpobar;
ibar      = delta * kbar;
tauc      = 0.1;
tauw      = 0.25;
tauk      = 0.15;
cbar      = (1 - obar*mpobar - delta*kbar - tauw*mplbar*lbar - tauk*kbar*mpkbar)/(1+tauc);
sigma_c   = 2;            
phiL      = 0.5; 
PhiB      = 0.001;        
rho_a     = 0.9;           
lambdabar = cbar^(-sigma_c)/(1+tauc);
gamma     = mplbar*(1-tauw)*lambdabar/(lbar^(1/phiL));

thetaN = 0.972;
theta  = 0.381;  

// CCyB PARAMETERS
phi1_ss = 0.02;          // Steady state capital requirement (TEMPORARILY 0 FOR TESTING)
phi2    = 0.01;           // CCyB responsiveness (set to 0 to disable CCyB)

beta  = 1/(mpkbar + (1-delta));
r_ss  = 1/beta - 1;

// STEADY STATE WITH CAPITAL REQUIREMENTS
credit_ss = kbar;
KBss      = phi1_ss * credit_ss;                    // KB_ss = phi1_ss * L_ss
nu_ss     = (mpkbar + 1 - delta)/(1 + r_ss);

// STANDARD GK leverage constraint (unchanged)
Nss       = (theta/nu_ss) * kbar;

// Deposits from balance sheet: q*k = N + D + KB
Dss       = kbar - Nss - KBss;

// Nbar must satisfy net worth evolution in steady state
// N_ss = thetaN * (Rk_ss*k - (1+r_ss)*D_ss) + (1-thetaN)*Nbar
Rk_ss = nu_ss * (1 + r_ss);
Nbar = (Nss - thetaN * (Rk_ss*kbar - (1 + r_ss)*Dss)) / (1 - thetaN);

//==============================================================
// MODEL
//==============================================================

model;

// Production
a = rho_a * a(-1) + eps_a;
ymarg = alphaK*k(-1)^rhoc + alphaL*l^rhoc + alphaO*oil^rhoc;
yhat  = ymarg^(1/rhoc);
y     = Abar * exp(a) * yhat;
mpk   = Abar * exp(a) * (ymarg^(1/rhoc - 1)) * alphaK * k(-1)^(rhoc - 1);
mpoil = Abar * exp(a) * (ymarg^(1/rhoc - 1)) * alphaO * oil^(rhoc - 1);
w     = Abar * exp(a) * (ymarg^(1/rhoc - 1)) * alphaL * l^(rhoc - 1);
Poilbar = mpoil;

// Capital
k = (1 - delta)*k(-1) + i;
i = delta * k(-1) + (q - 1) * (k(-1)/Phi);
Ck = (Phi/2) * ((i - delta*k(-1))/k(-1))^2 
     - Phi * ((i - delta*k(-1))/k(-1)) * (i/k(-1));

// Household
lambda = c^(-sigma_c)/(1+tauc); 
lambda = beta * (1 + r) * lambda(+1);
r = r_ss - PhiB * nfa(-1);
l = (w*(1-tauw)*lambda/gamma)^phiL; 

// FINANCIAL SECTOR - WITH TIME-VARYING CAPITAL REQUIREMENTS

// Credit gap (log deviation from steady state)
credit_gap = log(q(-1)*k(-1)) - log(credit_ss);

// CCyB Rule: capital requirement responds to credit gap
phi1 = phi1_ss + phi2 * credit_gap;

// Regulatory capital required
KB = phi1 * q * k;

// STANDARD GK LEVERAGE CONSTRAINT (unchanged)
N = (theta/nu) * q * k;

// Balance sheet identity: deposits are residual
D = q*k - N - KB;

// Net worth evolution - KEEP SIMPLE FORM, D lagged is just a variable
N = thetaN * (Rk*q(-1)*k(-1) - (1 + r(-1))*D(-1)) + (1 - thetaN)*Nbar;

// Return on capital
Rk = ((1 - delta)*q + mpk - Ck) / q(-1);

// No-arbitrage condition for nu
((1 - delta)*q(+1) + mpk(+1) - Ck(+1)) / q = nu * (1 + r);

// Spread (GK friction only - CCyB works through phi1)
spr = log(nu);

// Leverage ratio
leverage = q*k / N;

// Current account
nfa = (1 + r) * nfa(-1) + y - (1+tauc)*c - i - Poilbar*oil 
      - tauw*w*l - tauk*mpk*k(-1);

end;

//==============================================================
// STEADY STATE
//==============================================================
steady_state_model;
    a     = 0;
    k     = kbar;
    l     = lbar;
    oil   = obar;
    yhat  = yhatbar;
    ymarg = ymargbar;
    y     = 1; 
    mpk   = mpkbar;
    mpoil = mpobar;
    w     = mplbar;
    nfa   = 0;
    r     = r_ss;
    q     = 1;
    i     = delta * k;
    Ck    = 0;
    c     = cbar;
    lambda = lambdabar;
    
    // Financial sector steady state
    phi1  = phi1_ss;
    KB    = KBss;
    N     = Nss;
    D     = Dss;
    nu    = nu_ss;
    Rk    = Rk_ss;
    spr   = log(nu_ss);
    leverage = kbar / Nss;
    credit_gap = 0;
end;

//==============================================================
// SHOCKS
//==============================================================
shocks;
    var eps_a; stderr 0.01;
end;

//==============================================================
// SOLUTION
//==============================================================
steady;
check;

stoch_simul(order=1,irf=20, periods=1000)
    y c i k q N nu spr Rk leverage phi1 KB credit_gap D;

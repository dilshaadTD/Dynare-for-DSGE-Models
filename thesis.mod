//////////////////////////////////////////////////////////////////////////////
//This model is a DSGE model with Financial Accelerator Mechanism. It is
//developed in order to examine the effectiveness of reserve requirements on
//current account deficit in small open economies. 
//Dilþat Tugba Dalkýran
//Bilkent University 2011-2012
//Master's Thesis
//////////////////////////////////////////////////////////////////////////////
var c, i, g, ex, y, rer, i_D, pi, r_K, i_IB, q, k, n, a, h, pi_d, L, i_R, rrMP, rr, i_f, s, B, i_L, CA, d, br, pp, pph, D_IB; 
varexo e_a, e_ex, e_g, e_i_f, e_rr, e_mon;
parameters beta, psi1, psi2,delta, alpha, epsilon,tilda, phi, v, gamma, mu,psiB, X, eta, theta;

psi1     = -0.015;                                                  //parameter for the deposit banks cost of deviating from required reserve
psi2     = 0.01;                                                    //parameter for the deposit banks cost of deviating from required reserve
delta    = 0.025;                                                  // depreciation rate
beta     = 0.985;                                                   //discount factor
alpha    = 0.33;                                                   //capital share in production
phi      = 1;                                                        //inverse of Frish Labor elasticity
theta    = 0.75;                                                   // calvo
v        = 0.97;                                                       //Entrepreneur's Survival rate
X        = 0.25;                                                       // sensitivity of changes in capital to fluctuations in I/K ratio
eta      = 0.05;                                                     //elasticity comes from contract problem
psiB     = 1;                                                    // sensitivity parameter for cost of holding bond
gamma    = 0.75;                                                   //share of domestic consumption
epsilon  = 6.5;                                                  //markup=e/(e-1) number is from g,g,natalucci (2007)
tilda    = 1;                                                    //sensitivity of labor supply in utility function 
mu       = 0.12;                                                 //fraction of auditing cost


///////////////////////////////MODEL//////////////////////////////////////////////////////////////////////////////////////////////////
model (linear);

# i_Dss   = 1/beta;                                              //interest rate on deposits
# rrss    = 0.1;                                                //reserve ratio that deposit banks decides to holds
# rrMPss  = 0.1;                                                //reserve requirements that CB dictates
# Grr     = psi1*(rrss-rrMPss) + (psi2/2)*(rrss-rrMPss)^2;       //cost of deviating from ss RR
# i_IBss  = i_Dss + rrss*(0.015)+ Grr;                            //interbank rate
# i_Rss   = i_IBss-0.015;                                        //spread between i_IB and i_R
# by      = 0.4;
# qss     = 1;                                                   // price of capital
# S       = 1;                                                   //exchange rate
# gy      = 0.22;                                                // goverment share in production   G/Y
# i_Lss   = (i_IBss-rrMPss*i_Rss)/(1-rrMPss);        //interest rate on loans
# r_Kss   = 1.005*i_Lss;                                        //external finance premium
# yk      = (r_Kss+delta-1)/alpha;                               //Y/K
# iy      = delta / yk;                                          //investment share in GDP  I/Y
# leverage= 0.88;                                                //leverage ratio nss/qss*kss
# ly      = (1-leverage)*qss/yk;                                 //loans to GDP ratio  L/Y
# i_fss   = i_IBss;  //escude 2007 (argem)
# omega   = 0.01;
# adjcosttoy = (yk^(-1))*((X/2)*((iy/(yk^-1))-delta)^2 + omega*mu*r_Kss*qss)+Grr; 
# dy      = ly/((1-rrss)^2) - S*by/(1-rrss);                      // deposits to GDP ratio D/Y
//# exy     = (i_fss-1)*by + adjcosttoy;                 // X/Y
//#  cy = 1-alpha+1/epsilon-gy+((1-rrss)*i_IBss+2*rrss*i_Rss-1-rrss);

# nexy = (i_fss-1)*S*by;
# cy      = (1 - (iy+gy + nexy + adjcosttoy));
# exy = 1-gamma*(cy+iy+gy)-gamma*adjcosttoy;
# hss     = (cy^(1/(-phi-1)))*((1-alpha)/tilda)^(1/(1+phi));      //labor supply
# ass     = 1;                                                    //technology in ss
# yss     = ass^(1/(1-alpha)) * yk^(alpha/(alpha-1)) *hss;        //GDP at ss
# zss     = alpha* yk;                                            //rental rate
# kss     = yss / yk;
# gbar    = 0;                                                    //transfers from dying entrepreneurs
# Vss     = qss*kss*leverage/v - (1-v)*gbar/v;                    //worth of surviving entrepreneurs
# p_h     = 1;
# p_f     = 1;
# p       = 1;
# Lss     = ly*yss;
# nss     = leverage*kss*qss;
# Bss     = by*yss;
# Dss     = dy*yss;
# iss     = iy*yss;
# adjcostss = adjcosttoy*yss;  
# my      = exy-nexy; 

                           //Aggregate Demand Equations//
y = gamma*(cy+gy+iy)*pp + gamma*(cy*c + iy*i + gy*g) + exy*(s+ex) - pph;
//y=gamma*(cy+gy+iy)*pp+gamma*(cy*c + iy*i + gy*g)+exy*(s+ex)+gamma*adjcosttoy*(pp+(1/adjcostss)*(kss*mu*r_Kss*omega*(k(-1)+r_K+q(-1)) + psi1*(rrss*rr-rrMPss*rrMP)+psiB*Bss*(s + B - )*pph))-pph;
//y = gamma*(cy*c + iy*i + gy*g + (1-gamma)*rer)+(1-gamma)*(rer + ex); // goods market clearing
c = c(+1) - i_D + pi(+1);                                     //Consumer's Euler Equation
r_K(+1) - i_L + pi(+1) = eta*(q+k-n);                        //Leverage and External Finance Premium
r_K = (y-k(-1))*zss/r_Kss + q*(1-delta)/r_Kss + q(-1);               //Price of Capital
q = X*(i - k(-1));                                                   //Investment Demand
                       
                           //Aggregate Supply Equations//
y = a + alpha*k(-1) + (1-alpha)*h;                                   //Production Function
y - h - c = phi*h;                                                   //Labor Supply
pi_d = beta*pi_d(+1) + ((1-theta)*(1-theta*beta)/theta)*(alpha*(y-k(-1)) + (1-alpha)*(y-h) - a); //Phillip's Curve

                         //Evolution of State Variables//
k = delta*i + (1-delta)*k(-1);                                       //Evolution of Capital

nss*n = r_Kss*v*(1-mu)*kss*(r_K+k(-1)+q(-1)) - i_Lss*v*Lss*(i_L+pp(-1)-pp + L(-1));
//n = v*n(-1) + r_K +L*v*Lss/nss + (1-v*(Lss+nss)/nss)*(q(-1)+k(-1));  //Entrepsreneurs: Net Worth

                         //The others//

i_IBss*i_IB = i_Rss*i_R - psi2*(rr*rrss-rrMP*rrMPss);              //reserve requirements
i_Dss*i_D = (1-rrss)*i_IBss*i_IB - rrss*i_IBss*rr + rrss*i_Rss*(i_R + rr) - psi1*(rrss*rr - rrMPss*rrMP); //Deposit Rate
q + k = (Lss/(Lss+nss))*L + (nss/(Lss+nss))*n;                          //Entrepreneurs' Balance Sheet
rer - rer(-1) = s - s(-1) - pi_d;                                               //Real Exchange Rate
pi = pi_d*gamma + (1-gamma)*(s-s(-1));                                      //CPI-Inflation Rate
//(yss/(Bss*i_fss))*(y + pph -(cy+iy+gy)*pp-(cy*c+gy*g+iy*i)+by*(s+B)-adjcosttoy*(pp+(1/adjcostss)*(kss*mu*r_Kss*omega*(k(-1)+r_K+q(-1)) + psi1*(rrss*rr-rrMPss*rrMP)+psiB*Bss*(s + B - pp)*pph)))=s+B(-1)+i_f(-1);     //Balance of payment identity
(yss/(Bss*i_fss)) * (y + pph - (cy+iy+gy)*pp - (cy*c+gy*g+iy*i) + by*(s+B)) = s + B(-1) + i_f(-1); //Balance of payment identity

// if Banks borrow from abroad
i_IBss*i_IB - psiB*Bss*(br) = i_fss*(s(+1)- s + i_f); //uncovered interest parity condition          
Lss*L = S*Bss*((1-rrMPss)*(B + s) - rrMPss*rrMP)+Dss*((1-rrMPss)*(-rrss*rr+(1-rrss)*d)-(1-rrss)*rrMPss*rrMP); //lending bank's balance sheet
//i_L = i_IB*i_IBss-i_Rss*rrMPss*(i_R+rrMP)/(i_IBss-i_Rss*rrMPss)+rrMPss*rrMP/(1-rrMPss);


i_R = i_IB; //spread
CA  = -B + B(-1);  //liability olarak olarak tanýmlandýðý için iþaret deðiþtirir
br  = s - s(-1) + B - B(-1) - pi + br(-1);
d = D_IB;

pp  = pi + pp(-1);
pph = pi_d + pph(-1);


//Policies
i_IB = 3.3*pi + 0.3*y + e_mon;
rrMP = 1.9*L+e_rr; //financial and price stability objective

  //i_IB = 2.5*pi + 0.2*y;
  //rrMP = 19.7*L+e_rr; //only price stability objective

   //Shocks//
a = 0.89*a(-1)+e_a;
ex = 0.80*ex(-1)+e_ex;
g = .86*g(-1)+e_g;
i_f = .88*i_f(-1)+e_i_f;

end;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
check;

steady;

//initval;
//end; 

shocks;
//var e_a = 1.13;   //technology shock
//var e_i_f = 0.43;     //foreign interest rate shock
//var e_ex = 5.01;  //export shock
//var e_g = 4.63;  //government shock
var e_rr = 1.63;   // shock to reserve requirements
//var e_mon = 1.5;  // monetary policy shock
end;

stoch_simul(order=1,irf=30) i_L, i_R, i_IB, pi, y, CA, B, k, s, L, D_IB, i_D,n;
//stoch_simul(order=1,irf=40,periods=100000) CA, i_D, i_L, i_R, i_IB, rrMP, br, pi;
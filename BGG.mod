

// Fabio Verona, Ines Drummond, Manuel Mota Freitas - FEP - University of Porto (Portugal)
// May 08

// BGG model, code for Dynare

var cU , iU , yU , rnU , kU , xU , qU , hU , piU , rU , efpU , nU , ceU , rkU , aS,gS, mpS;
varexo e_mpS,e_gS,e_aS;

parameters  cy,cey,iy,gy,kn,nk1,yn,markup,Rk,epsilon,bbeta,eta,alphaK,delta,vsigma,rho,vphi,kappa,omega,v,gamma,sigma, R;

cy	=	0.605770191	; // C/Y
cey	=	0.01	; // Ce/Y
iy	=	0.184229809	; //% I/Y
gy	=	0.2	; //% G/Y
kn	=	2.081759973	; //%K/N
nk1	=	0.480362776	; //%N/K
yn	=	0.282494996	; //%Y/N
markup	=	1.1	; //%X
Rk	=	1.018177298	;
epsilon	=	0.957593537	;
bbeta	=	0.99	;
eta	=	3	;
//beta1 = 1 ; %coeficiente associado a H na funçao utilizadade (se=1 => ln (1-H))
alphaK	=	0.35	;
delta	=	0.025	;
vsigma=0.11;
rho	=	0.9	;
vphi	=	0.25	;
kappa	=	0.08583	;
omega	=	0.984615	;
v	=	0.052092347	;
//v=0;
gamma	=	0.9728	;
sigma	=	1	;

R = 1/bbeta ;


model;
efpU = rkU(+1) - rU;         //1
gS = 0.95*gS(-1) + e_gS;     //2
aS = aS(-1)  + e_aS;         //3                                           
mpS = e_mpS;                 //4
yU = cy*cU + iy*iU + gy*gS + cey*ceU;  //5
sigma*cU = -rU + sigma*cU(+1);         //6 
ceU = nU;                              //7
rkU(+1) - rU = -v*(nU-qU-kU);          //8
rkU = (1-epsilon)*(yU-kU(-1)-xU) + epsilon*qU - qU(-1);  //9
qU = vphi*(iU-kU(-1));             //10
yU = aS + alphaK*kU(-1)+(1-alphaK)*omega*hU;   //11
yU - hU - xU - cU = eta^(-1)*hU;               //12
piU = -kappa*xU + bbeta*piU(+1);               //13
kU = delta*iU + (1-delta)*kU(-1);              //14
nU = gamma/bbeta*nU(-1) + ((gamma/bbeta)-((gamma*kn)/bbeta))*rU(-1) + ((gamma*kn)/bbeta+((gamma*kn)*(Rk-(1/bbeta))))*rkU +((gamma*kn)*(Rk-(1/bbeta)))*qU(-1) + ((gamma*kn)*(Rk-(1/bbeta)))*kU(-1) + ((1-alphaK)*(1-omega)*yn/markup)*(yU-xU);
rnU = rho*rnU(-1) + vsigma*piU(-1) - mpS;   //16
rnU = rU + piU(+1);  //17
end;


//initval;

//end;


shocks; 
var e_gS; stderr 0.0155;
//var e_gS; stderr 0;
var e_aS; stderr 0.0043;
//var e_aS; stderr 0;
var e_mpS;stderr 0.06;
end;

check;
steady(solve_algo=2);

stoch_simul(order=1,irf=500) yU iU efpU nU rnU piU rU rkU;





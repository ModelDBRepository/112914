%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% "Capturing the bursting dynamics of a two-cell inhibitory network 
%                   using a one-dimensional map"
%
%       Victor Matveev (1), Amitabha Bose (1), Farzan Nadim(1,2)             
%      (1) Dept Math Sci, NJIT (2) Dept Bio Sci, Rutgers-Newark
%  
%  Dynamic variable order: Y = [V1 V2 w1 w2 h1 h2 s1 s2]
%  Parameter array: [ gbarsyn gtbar Vh tauSyn tgamma tlo thi vthresh ]
%
%                         February 3, 2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dydt = burstODE(t, Y, Params)
global Vh gbarsyn gtbar;

% See Eqs. 2-6 and Appendix A1 for parameter list

gbarsyn = Params(1); gtbar = Params(2); Vh = Params(3);
tauSyn = Params(4); tgamma = Params(5);
tlo = Params(6); thi = Params(7);
vthresh = Params(8);

VK=-84;     VL=-60;     VCa=120;    VNa=180;
C=2;        phi=0.6667; V1=-12;
V2=18;      V3=-8;      V4=6.0;
gCa=4.0;    gK=8.0;     gL=2.0;
Cinv=1.0/C; einh=-80;
vax=Vh;     vix=Vh;   
slope=4;    zap=14;

minf=0.5*(1+tanh((Y(1)-V1)/V2));
winf=0.5*(1+tanh((Y(1)-V3)/V4));
kW=cosh((Y(1)-V3)/(2*V4));
minf1=0.5*(1+tanh((Y(2)-V1)/V2));
winf1=0.5*(1+tanh((Y(2)-V3)/V4));
kW1=cosh( (Y(2) - V3) / (2 * V4) );
ax=0.5*(1+tanh(slope*(Y(1)-vax)));
ax1=0.5*(1+tanh(slope*(Y(2)-vax)));

dydt=[Cinv*(-gCa*minf*(Y(1)-VCa) - gK*Y(3)*(Y(1)-VK) - gL*(Y(1)-VL) - gtbar*ax*Y(5)*(Y(1)-VCa) - gbarsyn*Y(8)*(Y(1)-einh) + zap);
Cinv * (-gCa*minf1*(Y(2)-VCa) - gK*Y(4)*(Y(2)-VK) - gL*(Y(2)-VL) - gtbar*ax1*Y(6)*(Y(2)-VCa) - gbarsyn*Y(7)*(Y(2)-einh) + zap);
phi * (winf - Y(3)) * kW; phi * (winf1 - Y(4)) * kW1;
(1-Y(5))*0.5*(1 + tanh(slope*(vix-Y(1))))/tlo - Y(5)*0.5*(1 + tanh(slope*(Y(1)-vix)))/thi;
(1-Y(6))*0.5*(1 + tanh(slope*(vix-Y(2))))/tlo - Y(6)*0.5*(1 + tanh(slope*(Y(2)-vix)))/thi;
(1-Y(7))*0.5*(1 + tanh(slope*(Y(1)-vthresh)))/tgamma - Y(7)*0.5*(1 + tanh(slope*(vthresh-Y(1))))/tauSyn;
(1-Y(8))*0.5*(1 + tanh(slope*(Y(2)-vthresh)))/tgamma - Y(8)*0.5*(1 + tanh(slope*(vthresh-Y(2))))/tauSyn; ];


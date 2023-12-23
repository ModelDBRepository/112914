%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% "Capturing the bursting dynamics of a two-cell inhibitory network
%                   using a one-dimensional map"
%      Victor Matveev (1), Amitabha Bose (1), Farzan Nadim(1,2)
%      (1) Dept Math Sci, NJIT (2) Dept Bio Sci, Rutgers-Newark
%
%                         February 3, 2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function wlist = Vnullcline(interval, Vex, h, s)
global gbarsyn gtbar Vh;

VK=-84;     VL=-60;     VCa=120;
C=2;        phi=0.6667; V1=-12;
V2=18;      V3=-8;      V4=6.0;
gCa=4.0;    gK=8.0;     gL=2.0;
Cinv=1.0/C; einh=-80;  
slope=4;    zap=14;

ax=0.5*(1+tanh(slope*(Vex-Vh)));
wsol = 0; wlist = [];

for V = interval
    minf=0.5*(1+tanh((V-V1)/V2));
    Vprime = @(x) Cinv*(-gCa*minf*(V-VCa) - gK*x*(V-VK) - gL*(V-VL) - gtbar*ax*h*(V-VCa) - gbarsyn*s*(V-einh) + zap);
    wsol = fzero(Vprime, wsol);
    wlist = [wlist, wsol];
end;




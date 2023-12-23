%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% "Capturing the bursting dynamics of a two-cell inhibitory network 
%                   using a one-dimensional map"
%      Victor Matveev (1), Amitabha Bose (1), Farzan Nadim(1,2)             
%      (1) Dept Math Sci, NJIT (2) Dept Bio Sci, Rutgers-Newark
%  
%             Multistable bursting (Figs. 1, 10 and 11)
%                         February 3, 2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global gbarsyn gtbar Vh;

% List of initial conditions for different bursting states
ICS(1,1:8) = [-23.1 -54.1 0.00104 0      0.0855 0.111  0.107 0 ];
ICS(2,1:8) = [-38.7 -56.9 0.00848 0      0.122  0.0774 0.303 0 ];
ICS(3,1:8) = [-58.3 -34.1 0       0.425  0.0951 0.126  0     0.647 ]; 
ICS(4,1:8) = [-56    5.9  0       0.58   0.161  0.0958 0     0.916 ]; 
ICS(5,1:8) = [-58.1 -18.6 0       0.591  0.137  0.121  0     0.78 ]; 

SN = 0;
while SN < 7 || SN > 11
   SN = input('*** How many spikes per burst? (7, 8, 9, 10 or 11):_');
end;
IC = ICS(SN-6,:);

Period = 50; 
twindow = Period * 1.2;
T = 2 *(Period + twindow);
figure(1); set(1, 'position', [250, 500, 640, 320]);

options = [];
options = odeset(options,'RelTol',1e-4);

Params = [1.1 1.38 -52 1 0.2 100 20 -3];
[t,y] = ode45(@burstODE, [0 T], IC, options, Params);

Vinterval = [-65:2:30];
Yw = 0.5*(1+tanh((Vinterval+8)/6));

i = find(t>twindow/2, 1);

while t(i) < T - twindow/2;
    i = i + 6;
    subplot('position',[0.1 0.15 0.45 0.78]);
    hold off;
    plot(Vinterval, Yw, 'b-', 'linewidth', 2);
    xlabel('V','fontsize',14); ylabel('w','fontsize',14);
    hold on;
    plot([Vh Vh], [-0.1 0.8], 'r--', 'linewidth', 2);
    
    V1 = y(i,1); w = y(i,3);
    plot(V1, w,'om', 'linewidth', 2);
    V2 = y(i,2); w = y(i,4);
    plot(V2, w,'ok', 'linewidth', 2);

    ix = y(i,5); s  = y(i,8); 
    Y = Vnullcline(Vinterval, V1, ix, s);
    plot(Vinterval, Y, '-m', 'linewidth', 2);
    ix = y(i,6); s  = y(i,7);
    Y = Vnullcline(Vinterval, V2, ix,s);
    plot(Vinterval, Y, '-k', 'linewidth', 2);

    axis([-65 30 -0.1 0.7]);
    title(['\bf Inhibitory network with T-currents: ', num2str(SN), '-spike burst']);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subplot('position', [0.6 0.6 0.37 0.25]);
    hold off;
    tt = t - t(i) + twindow/2;
    vmin = min(y(:,1))*1.1; vmax = max(y(:,1))*1.1;
    plot([twindow/2 twindow/2], [vmin vmax], 'y-', 'linewidth', 2); hold on;
    plot([0 twindow], [Vh Vh], 'b-', 'linewidth', 2);
    plot(tt,y(:,1),'m-','linewidth',1);
    plot(tt,y(:,2),'k-','linewidth',1);
    axis([0 twindow vmin vmax]);
    title('Voltage: V_1(t), V_2(t)','fontsize',12);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    hmax = max(y(:,5))*1.1;
    subplot('position', [0.6 0.15 0.37 0.25]);
    hold off;
    tt = t - t(i) + twindow/2;
    plot([twindow/2 twindow/2], [0 hmax], 'y-', 'linewidth', 2);
    hold on;
    plot(tt,y(:,5),'m-','linewidth',1);
    plot(tt,y(:,6),'k-','linewidth',1);
    axis([0 twindow 0 hmax]);
    title('Inactivation: h_1(t), h_2(t)','fontsize',12);
    xlabel('time (ms)', 'fontsize', 12);
    drawnow;

end;

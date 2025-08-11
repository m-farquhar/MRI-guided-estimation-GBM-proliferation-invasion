%% Creates Figure 5, all data combinations
% Figure5.eps
%
% Files required:
% mouseData_radius.mat
% mouse1dataopt_parameters_gauss.mat
% mouse2dataopt_parameters_gauss.mat
% mouse3dataopt_parameters_gauss.mat

%% Close all figures and clear workspace
close all
clear all

%% Plot properties
linewidth = 1.25;
fontsize = 10;
interpreter = 'Latex';
markersize = 7;
colour1 = [0, 0.447 0.7410];
colour2 = [0.85,0.325,0.098];
colour3 = [0.9290,0.6940,0.1250];
colour4 = [0.4940,0.1840,0.5560];
colour5 = [0.4660,0.6740,0.1880];
colour6 = [0.3010,0.7450,0.9330];
colour7 = [0.6350,0.0780,0.1840];
marker1 = '*';
marker2 = 'o';
marker3 = 'v';
marker4 = '^';
marker5 = '<';
marker6 = '>';
marker7 = 'x';
marker8 = '+';
marker9 = 's';
marker10 = 'd';
marker11 = 'p';
marker12 = 'h';
ax = [0,25,0,4];
legToken = [10,18];

L = 5;
N = 10001;
dr = L/(N-1);
r = 0:dr:L;

V0 = 2;
N0 = 100000;

c0 = N0/V0;
cstar = 8000;
b = (V0/pi^(3/2))^(1/3);

figure
pos = get(gcf, 'position');
pos(3) = pos(3)*1.05;
pos(4) = pos(4)*1.05;
set(gcf, 'position',pos)
t1 = tiledlayout(3,3);
t1.TileSpacing = 'compact';
t1.Padding = 'none';
%% Load Data
load('mouseData_radius.mat')
ms1opt = load('mouse1dataopt_parameters_gauss');

tspan1 = [7,14,19,21];
raddat1 = rall(1,:);

ms2opt = load('mouse2dataopt_parameters_gauss');

tspan2 = [7,14,19,21];
raddat2 = rall(2,:);

ms3opt = load('mouse3dataopt_parameters_gauss');

tspan3 = [7,14,19];
raddat3 = rall(3,1:3);


%% Generate best fits using Fisher's formula
x0_1 = max(ms1opt.rhoOptV);
v1 = (raddat1(end) - raddat1(end-1))/(tspan1(end) - tspan1(end-1));
tol = 1e-3;
options = optimoptions(@lsqcurvefit,'DiffMinChange', tol, 'display','iter');
[rhoBest1, resnorm1] = lsqcurvefit(@(x,xdata) b*analytic_sol_gauss_r_nondim(v1.^2/4/x/x/b^2, x*xdata, cstar/c0), x0_1, tspan1, raddat1, 0, inf, options);
Dbest1 = v1.^2/4/rhoBest1;
alpha1 = Dbest1/rhoBest1/b^2;
tau1 = rhoBest1*tspan1;
beta = cstar/c0;
rad1 = analytic_sol_gauss_r_nondim(alpha1,tau1,beta);

x0_2 = max(ms2opt.rhoOptV);
v2 = (raddat2(end) - raddat2(end-1))/(tspan2(end) - tspan2(end-1));
tol = 1e-3;
options = optimoptions(@lsqcurvefit,'DiffMinChange', tol, 'display','iter');
[rhoBest2, resnorm2] = lsqcurvefit(@(x,xdata) b*analytic_sol_gauss_r_nondim(v2.^2/4/x/x/b^2, x*xdata, cstar/c0), x0_2, tspan2, raddat2, 0, inf, options);
Dbest2 = v2.^2/4/rhoBest2;
alpha2 = Dbest2/rhoBest2/b^2;
tau2 = rhoBest2*tspan2;
rad2 = analytic_sol_gauss_r_nondim(alpha2,tau2,beta);

x0_3 = max(ms3opt.rhoOptV);
v3 = (raddat3(end) - raddat3(end-1))/(tspan3(end) - tspan3(end-1));
tol = 1e-3;
options = optimoptions(@lsqcurvefit,'DiffMinChange', tol, 'display','iter');
[rhoBest3, resnorm3] = lsqcurvefit(@(x,xdata) b*analytic_sol_gauss_r_nondim(v3.^2/4/x/x/b^2, x*xdata, cstar/c0), x0_3, tspan3, raddat3, 0, inf, options);
Dbest3 = v3.^2/4/rhoBest3;
alpha3 = Dbest3/rhoBest3/b^2;
tau3 = rhoBest3*tspan3;
rad3 = analytic_sol_gauss_r_nondim(alpha3,tau3,beta);


%% Plot Curves using all data points
nexttile

hold on
fplot(@(x) b*analytic_sol_gauss_r_nondim(ms1opt.DoptV(1)/ms1opt.rhoOptV(1)/b^2, ms1opt.rhoOptV(1)*x, cstar/c0),[0,25], 'linewidth', linewidth, 'color', colour1)
plot(tspan1, raddat1, ['k', marker1], 'markersize', markersize, 'linewidth', linewidth)

fplot(@(x) b*analytic_sol_gauss_r_nondim(Dbest1/rhoBest1/b^2, rhoBest1*x, cstar/c0),[0,25], '--', 'linewidth', linewidth, 'color', [0,0,0]) %optimal using Fisher formula
axis(ax)
ytickformat('%.1f')
set(gca, 'fontsize',fontsize, 'ticklabelinterpreter', 'latex')
fig1 = get(gca, 'position');
title('\textbf{Mouse 1}', 'fontsize', fontsize, 'interpreter', interpreter)
text(0,4.2, '\textbf{A}', 'fontsize', fontsize,'interpreter', interpreter)
box on

nexttile
hold on
fplot(@(x) b*analytic_sol_gauss_r_nondim(ms2opt.DoptV(1)/ms2opt.rhoOptV(1)/b^2, ms2opt.rhoOptV(1)*x, cstar/c0),[0,25], 'linewidth', linewidth, 'color', colour1)
plot(tspan2, raddat2, ['k', marker1], 'markersize', markersize, 'linewidth', linewidth)

fplot(@(x) b*analytic_sol_gauss_r_nondim(Dbest2/rhoBest2/b^2, rhoBest2*x, cstar/c0),[0,25], '--', 'linewidth', linewidth, 'color', [0,0,0]) % Optimal using Fisher formula 
axis(ax)

ytickformat('%.1f')
set(gca, 'fontsize',fontsize, 'ticklabelinterpreter', 'latex')
fig1 = get(gca, 'position');
box on
title('\textbf{Mouse 2}', 'fontsize', fontsize, 'interpreter', interpreter)
text(0,4.2, '\textbf{D}', 'fontsize', fontsize,'interpreter', interpreter)

nexttile

plot(-10,-10,['k',marker1],'markersize',markersize,'linewidth',linewidth)
hold on
plot(-10,-10, 'linewidth',linewidth, 'color', colour1)
axis(ax)
axis off
lg = legend('Data', '1,2,3,4', 'fontsize',fontsize, 'interpreter', interpreter, 'location', 'eastoutside');
title('\textbf{Mouse 3}', 'fontsize', fontsize, 'interpreter', interpreter)
lg.ItemTokenSize = legToken;


%% Plot Curves using 3 data points

nexttile
hold on
colmat = [colour1;colour2;colour3;colour4];
for i = 1:4
    if ms1opt.DoptV(i+1)~= 0
        fplot(@(x) b*analytic_sol_gauss_r_nondim(ms1opt.DoptV(i+1)/ms1opt.rhoOptV(i+1)/b^2, ms1opt.rhoOptV(i+1)*x, cstar/c0),[0,25], 'linewidth', linewidth, 'color', colmat(i,:))
    end
end
plot(tspan1, raddat1, ['k', marker1], 'markersize', markersize, 'linewidth', linewidth)
fplot(@(x) b*analytic_sol_gauss_r_nondim(Dbest1/rhoBest1/b^2,rhoBest1*x, cstar/c0),[0,25], '--', 'linewidth', linewidth, 'color', [0,0,0])
axis(ax)
ytickformat('%.1f')
box on
set(gca, 'fontsize',fontsize, 'ticklabelinterpreter', 'latex')
text(0,4.2, '\textbf{B}', 'fontsize', fontsize,'interpreter', interpreter)

nexttile
hold on
colmat = [colour1;colour2;colour3;colour4];
for i = 1:4
    if ms2opt.DoptV(i+1)~= 0
        fplot(@(x) b*analytic_sol_gauss_r_nondim(ms2opt.DoptV(i+1)/ms2opt.rhoOptV(i+1)/b^2, ms2opt.rhoOptV(i+1)*x, cstar/c0),[0,25], 'linewidth', linewidth, 'color', colmat(i,:))
    end
end
plot(tspan2, raddat2, ['k', marker1], 'markersize', markersize, 'linewidth', linewidth)
fplot(@(x) b*analytic_sol_gauss_r_nondim(Dbest2/rhoBest2/b^2,rhoBest2*x, cstar/c0),[0,25], '--', 'linewidth', linewidth, 'color', [0,0,0])
axis(ax)
ytickformat('%.1f')
box on
set(gca, 'fontsize',fontsize, 'ticklabelinterpreter', 'latex')
text(0,4.2, '\textbf{E}', 'fontsize', fontsize,'interpreter', interpreter)

nexttile
hold on
colmat = [colour1;colour2;colour3;colour4];
plot(-10,-10,'k*', 'markersize',markersize,'linewidth',linewidth)
plot(-10,-10, 'linewidth',linewidth, 'color', colour1)
plot(-10,-10, 'linewidth',linewidth, 'color', colour2)
plot(-10,-10, 'linewidth',linewidth, 'color', colour3)
plot(-10,-10, 'linewidth',linewidth, 'color', colour4)

for i = 1:4
    if ms3opt.DoptV(i+1)~= 0
        fplot(@(x) b*analytic_sol_gauss_r_nondim(ms3opt.DoptV(i+1)/ms3opt.rhoOptV(i+1)/b^2, ms3opt.rhoOptV(i+1)*x, cstar/c0),[0,25], 'linewidth', linewidth, 'color', colmat(i,:))
    end
end
plot(tspan3, raddat3, ['k', marker1], 'markersize', markersize, 'linewidth', linewidth)
fplot(@(x) b*analytic_sol_gauss_r_nondim(Dbest3/rhoBest3/b^2,rhoBest3*x, cstar/c0),[0,25], '--', 'linewidth', linewidth, 'color', [0,0,0])
axis(ax)
ytickformat('%.1f')

lg = legend('Data', '1,2,3', '1,2,4', '1,3,4', '2,3,4', 'fontsize',fontsize, 'interpreter', interpreter, 'location', 'eastoutside');
lg.ItemTokenSize = legToken;
box on
set(gca, 'fontsize',fontsize, 'ticklabelinterpreter', 'latex')
text(0,4.2, '\textbf{G}', 'fontsize', fontsize,'interpreter', interpreter)

%% Plot curves using 2 data points
nexttile
plot(tspan1, raddat1, ['k', marker1], 'markersize', markersize, 'linewidth', linewidth)
hold on
colmat = [colour1;colour2;colour3;colour4;colour5;colour6];
for i = 1:6
    if ms1opt.DoptV(i+5)~= 0
        fplot(@(x) b*analytic_sol_gauss_r_nondim(ms1opt.DoptV(i+5)/ms1opt.rhoOptV(i+5)/b^2, ms1opt.rhoOptV(i+5)*x, cstar/c0),[0,25], 'linewidth', linewidth, 'color', colmat(i,:))
    end
end
plot(tspan1, raddat1, ['k', marker1], 'markersize', markersize, 'linewidth', linewidth)
ylim([0,3])
fplot(@(x) b*analytic_sol_gauss_r_nondim(Dbest1/rhoBest1/b^2, rhoBest1*x, cstar/c0),[0,25], '--', 'linewidth', linewidth, 'color', [0,0,0])
axis(ax)
ytickformat('%.1f')
set(gca, 'fontsize',fontsize, 'ticklabelinterpreter', 'latex')
box on
text(0,4.2, '\textbf{C}', 'fontsize', fontsize,'interpreter', interpreter)


nexttile
plot(tspan2, raddat2, ['k', marker1], 'markersize', markersize, 'linewidth', linewidth)
hold on
colmat = [colour1;colour2;colour3;colour4;colour5;colour6];
for i = 1:6
    if ms2opt.DoptV(i+5)~= 0
        fplot(@(x) b*analytic_sol_gauss_r_nondim(ms2opt.DoptV(i+5)/ms2opt.rhoOptV(i+5)/b^2, ms2opt.rhoOptV(i+5)*x, cstar/c0),[0,25], 'linewidth', linewidth, 'color', colmat(i,:))
    end
end
plot(tspan2, raddat2, ['k', marker1], 'markersize', markersize, 'linewidth', linewidth)
ylim([0,3])
%     ax = axis;
%     ax(3) = 0;
fplot(@(x) b*analytic_sol_gauss_r_nondim(Dbest2/rhoBest2/b^2, rhoBest2*x, cstar/c0),[0,25], '--', 'linewidth', linewidth, 'color', [0,0,0])
axis(ax)
ytickformat('%.1f')
set(gca, 'fontsize',fontsize, 'ticklabelinterpreter', 'latex')
box on
text(0,4.2, '\textbf{F}', 'fontsize', fontsize,'interpreter', interpreter)

nexttile

plot(-10,-10,'k*', 'markersize',markersize,'linewidth',linewidth)
hold on
plot(-10,-10, 'linewidth',linewidth, 'color', colour1)
plot(-10,-10, 'linewidth',linewidth, 'color', colour2)
plot(-10,-10, 'linewidth',linewidth, 'color', colour3)
plot(-10,-10, 'linewidth',linewidth, 'color', colour4)
plot(-10,-10, 'linewidth',linewidth, 'color', colour5)
plot(-10,-10, 'linewidth',linewidth, 'color', colour6)
markermat = {marker1, marker8, marker9, marker10, marker11, marker12};
colmat = [colour1;colour2;colour3;colour4;colour5;colour6];
for i = 1:6
    if ms3opt.DoptV(i+5)~= 0
        fplot(@(x) b*analytic_sol_gauss_r_nondim(ms3opt.DoptV(i+5)/ms3opt.rhoOptV(i+5)/b^2, ms3opt.rhoOptV(i+5)*x, cstar/c0),[0,25], 'linewidth', linewidth, 'color', colmat(i,:))
    end
end
plot(tspan3, raddat3, ['k', marker1], 'markersize', markersize, 'linewidth', linewidth)
ylim([0,3])
fplot(@(x) b*analytic_sol_gauss_r_nondim(Dbest3/rhoBest3/b^2, rhoBest3*x, cstar/c0),[0,25], '--', 'linewidth', linewidth, 'color', [0,0,0])
axis(ax)
ytickformat('%.1f')
set(gca, 'fontsize',fontsize, 'ticklabelinterpreter', 'latex')
box on
text(0,4.2, '\textbf{H}', 'fontsize', fontsize,'interpreter', interpreter)

lg = legend('Data', '1,2', '1,3','1,4','2,3','2,4','3,4', 'fontsize',fontsize, 'interpreter', interpreter, 'location', 'eastoutside');
lg.ItemTokenSize = legToken;

xlabel(t1,'$t$', 'fontsize', fontsize, 'interpreter', interpreter)
ylabel(t1,'$r^*$', 'fontsize', fontsize, 'interpreter', interpreter)

exportgraphics(gcf, 'Figure5.eps','contenttype','vector')

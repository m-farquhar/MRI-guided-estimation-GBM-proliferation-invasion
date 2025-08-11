%% Clear workspace
clear all

%% Point combinations
Points = {[1,2,3,4], [1,2,3], [1,2,4], [1,3,4], [2,3,4], [1,2], [1,3], [1,4], [2,3], [2,4], [3,4],[1,2,3,4], [1,2,3], [1,2,4], [1,3,4], [2,3,4], [1,2], [1,3], [1,4], [2,3], [2,4], [3,4],[2,3]};

%% Initial Parameters
N0 = 100000;
cstar = 8000;

%% data files
filenames = {'mouse1data', 'mouse2data','mouse3data'};
disp('$\mathcal{I}_k$ & $D$ & $\rho$ &  Max Error& $D$ & $\rho$ &  Max Error& $D$ & $\rho$ &  Max Error\\')

disp('\midrule')
load('mouseData_radius.mat')
tspan1 = [7,14,19,21];
tspan2 = [7,14,19,21];
tspan3 = [7,14,19];
rad1 = rall(1,:);
rad2 = rall(2,:);
rad3 = rall(3,1:3);
load([filenames{1}, 'opt_parameters_gauss.mat'], 'rhoOptV', 'DoptV', 'radOptV')
rhoOptV1 = rhoOptV;
DoptV1 = DoptV;
radOptV1 = radOptV;
MEV1 = max( abs( repmat(rad1,11,1) - radOptV1 ),[],2);

load([filenames{2}, 'opt_parameters_gauss.mat'], 'rhoOptV', 'DoptV', 'radOptV')
rhoOptV2 = rhoOptV;
DoptV2 = DoptV;
radOptV2 = radOptV;
MEV2 = max( abs( repmat(rad2,11,1) - radOptV2 ),[],2);

load([filenames{3}, 'opt_parameters_gauss.mat'], 'rhoOptV', 'DoptV', 'radOptV')
rhoOptV3 = rhoOptV;
DoptV3 = DoptV;
radOptV3 = radOptV;
MEV3 = max( abs( repmat(rad3,11,1) - radOptV3 ),[],2);

results = [DoptV1, rhoOptV1, MEV1, DoptV2, rhoOptV2, MEV2, DoptV3, rhoOptV3, MEV3];



for i = 1:size(results,1)
    Pts = Points{i};
          ptsstr = num2str(Pts,'%1.g,');
            ptsstr = ptsstr(1:end-1);           
    if any(Pts==4)
        disp(['$', ptsstr, '$ & '])%, num2str(Dopt, '%.3f'), ' & ', num2str(rhoOpt, '%.3f'),  ' & ', num2str(maxerr, '%.3f'), ' \\'])
        fprintf('%.3f & %.3f & %.3f & %.3f & %.3f & %.3f & - & - & - \\\\ \n', results(i,1:6))
    else
        disp(['$', ptsstr, '$ & '])%, num2str(Dopt, '%.3f'), ' & ', num2str(rhoOpt, '%.3f'),  ' & ', num2str(maxerr, '%.3f'), ' \\'])
        fprintf('%.3f & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f \\\\ \n', results(i,:))
    end
end
    disp('\bottomrule')

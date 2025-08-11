%% Plot properties
linewidth = 1;
fontsize = 10;
interpreter = 'Latex';
markersize = 5;
colour1 = [0, 0.447 0.7410];
colour2 = [0.85,0.325,0.098];
colour3 = [0.9290,0.6940,0.1250];
colour4 = [0.4940,0.1840,0.5560];
colour5 = [0.4660,0.6740,0.1880];
colour6 = [0.3010,0.7450,0.9330];
colourb = [0.15,0.15,0.15];
funstr = '$D/(\rho{b}^2)$';
funCol = @(D, rho,b) D./(rho)./b.^2;
deltavals1 = 0.35; % Mouse 1 delta plot values
deltavals2 = 0.35; % Mouse 2 delta plot values

%% Reduce saturation in colourmap for low values
cmapsize = 100;
% convert to HSV so we can damp down the saturation
colmat = rgb2hsv(reshape(cool(cmapsize), cmapsize,1,3));
colmat(:,1,2) = colmat(:,1,2) .*[linspace(0,1,cmapsize/2).^1' ; ones(cmapsize/2,1)];
colmat = squeeze(hsv2rgb(colmat));



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
legToken = [10,18];

%% boxplot labels
boxlabel1 = 'BI';
boxlabel2 = 'LS';

%% Model Parameters
V0 = 2;
N0 = 100000;
c0 = N0/V0;
cstar = 8000;
b = (V0/pi^(3/2))^(1/3);
tvec = linspace(0,25,100);


%% Plots
figure
pos=get(gcf,'position');
pos(3) = 1.0*pos(3);
pos(4) = 1.65*pos(4);
pos(2) = pos(2) - 0.4*pos(4);
set(gcf, 'position',pos)
t1 = tiledlayout(3,9);
t1.TileSpacing = 'tight';
t1.Padding = 'none';

axMs = [
    0.065, 0.12, 0.11, 0.15     % mouse 1
    0.065, 0.12, 0.13, 0.17    % mouse 2
    0.065, 0.12, 0.16, 0.20    % mouse 3
    ];

load mouseData_radius.mat
days = [7,14,19,21];

for mouse = 1:3
    %% Add mouse label to row
    nexttile
    axis([-1,1,-1,1])
    h=text(0,0, ['Mouse ', num2str(mouse)], 'fontsize',fontsize, 'interpreter', interpreter,'horizontalalignment','center', 'verticalalignment', 'bottom');
    set(h,'Rotation',90);
    axis off

    %% Plot parameter density maps from Tim's code
    if mouse == 1
        load mousefile C1 % Load Bayesian Inference results
        load mouse1dataopt_parameters_gauss DoptV rhoOptV % Load least squares results
        C = C1;
    elseif mouse == 2
        load mousefile C2 % Load Bayesian Inference results
        load mouse2dataopt_parameters_gauss DoptV rhoOptV % Load least squares results
        C = C2;
    elseif mouse == 3
        load mousefile C3 % Load Bayesian Inference results
        load mouse3dataopt_parameters_gauss DoptV rhoOptV % Load least squares results
        C = C3;
    else
        error('Bad mouse')
    end


    nbins = 55;
    nlevels = 50;
    nexttile([1,3]);
    [handle,N,Xedges,Yedges] = doplot(C, axMs(mouse,:), nbins, nlevels);

    colormap(colmat)
    clim([0,inf])
    hold on
    markervec = {marker2, marker3, marker4, marker5, marker6, marker7, marker8, marker9, marker10, marker11, marker12};
    colourvec = {colour1, colour1, colour2, colour3, colour4, colour1, colour2, colour3, colour4, colour5, colour6};
    for i = 1:11
        if DoptV(i) ~=0
            plot(DoptV(i), rhoOptV(i), markervec{i}, 'linewidth', linewidth, 'markersize', markersize, 'color',colourvec{i}) % Plot least squares solution on figure
        end
    end

    ylabel('$\rho$', 'interpreter', interpreter, 'fontsize', fontsize)
    xlabel('$D$', 'interpreter', interpreter, 'fontsize', fontsize)
    xlim(axMs(mouse,1:2))
    ytickformat('%.2f')
    set(gca, 'fontsize',fontsize, 'ticklabelinterpreter', 'latex')
    box on

    % Add subfigure labels and legend (legend is on mouse two so that it
    % includes all combination, could be included on mouse one instead
    if mouse==3
        title('\textbf{C}', 'interpreter', interpreter, 'fontsize', fontsize,'fontweight','bold')
        ax = gca;
        ax.TitleHorizontalAlignment = 'left';
    elseif mouse == 2
        title('\textbf{B}', 'interpreter', interpreter, 'fontsize', fontsize,'fontweight','bold')
        ax = gca;
        ax.TitleHorizontalAlignment = 'left';
        lg = legend({boxlabel1,'1,2,3,4', '1,2,3', '1,2,4', '1,3,4', '2,3,4', '1,2', '1,3', '1,4', '2,3', '2,4','3,4'} , 'fontsize',fontsize, 'interpreter', interpreter,'numColumns',6,'orientation','horizontal');
        lg.ItemTokenSize = legToken;
        lg.Layout.Tile = 'South';
    else
        title('\textbf{A}', 'interpreter', interpreter, 'fontsize', fontsize,'fontweight','bold')
        ax = gca;
        ax.TitleHorizontalAlignment = 'left';
    end
    drawnow

    %% Plot solution curves from parameter maps
    nexttile([1,3]);
    Xmids = (Xedges(1:end-1) + Xedges(2:end))/2;
    Ymids = (Yedges(1:end-1) + Yedges(2:end))/2;
    [YY,XX] = meshgrid(Ymids, Xmids);
    Nfine = N(:);
    [Nfine, i] = sort(Nfine(:));
    Nfine = Nfine(:,1);
    Nfine(Nfine<0.05*max(Nfine)) = nan;
    XX = XX(i);
    YY = YY(i);
    XX = XX(~isnan(Nfine));
    YY = YY(~isnan(Nfine));
    Nfine = Nfine(~isnan(Nfine));

    [T,D] = meshgrid(tvec, XX);
    [~,Rho] = meshgrid(tvec, YY);
    R = analytic_sol_gauss_r_dim(T, D, Rho, b, cstar,c0);
    CData = repmat(Nfine,1,length(tvec));
    contourf(T,R,CData, nlevels,'LineStyle','none')
    grid off
    axis([0,25.05,-0.02,4])
    clim([0,inf])
    set(gca,'clipping', 'on')
    box on
    xlabel('$t$', 'interpreter', interpreter, 'fontsize', fontsize)
    ylabel('$r^*$', 'interpreter', interpreter, 'fontsize', fontsize)
    set(gca, 'fontsize',fontsize, 'ticklabelinterpreter', 'latex')

    hold on
    plot(days, rall(mouse,:), 'k*', 'markersize', 8, 'linewidth', 1.)

    %% Add in best v^2/4/rho approximation
    x0 = max(rhoOptV);
    v = (rall(mouse,end) - rall(mouse,end-1))/(days(end) - days(end-1));
    if mouse==3
    v = (rall(mouse,end-1) - rall(mouse,end-2))/(days(end-1) - days(end-2));
    end
    tol = 1e-3;
    options = optimoptions(@lsqcurvefit,'DiffMinChange', tol, 'display','none');
    [rhoBest, resnorm2] = lsqcurvefit(@(x,xdata) analytic_sol_gauss_r_dim(xdata,v.^2/4/x, x, b, cstar,c0), x0, days(~isnan(rall(mouse,:))), rall(mouse,~isnan(rall(mouse,:))), 0, inf, options);
    Dbest = v.^2/4/rhoBest;
    disp([Dbest, rhoBest])
    fplot(@(x) analytic_sol_gauss_r_dim(x, Dbest,rhoBest,b,cstar, c0),[0,25], '--', 'linewidth', linewidth, 'color', [0,0,0])

    % Add subfigure labels
    if mouse ==1
        title('\textbf{D}', 'interpreter', interpreter, 'fontsize', fontsize,'fontweight','bold')
        ax = gca;
        ax.TitleHorizontalAlignment = 'left';
    elseif mouse==2

        title('\textbf{E}', 'interpreter', interpreter, 'fontsize', fontsize,'fontweight','bold')
        ax = gca;
        ax.TitleHorizontalAlignment = 'left';
    else
        title('\textbf{F}', 'interpreter', interpreter, 'fontsize', fontsize,'fontweight','bold')
        ax = gca;
        ax.TitleHorizontalAlignment = 'left';
    end
    set(gca, 'clipping', 'on', 'ClippingStyle','rectangle')
    drawnow


    %% Plot boxplots
    load("mouse" + mouse + "dataopt_parameters_gauss")
    reps = 100;
    alpha = C(:,2) ./ C(:,3) / b^2;
    alpha_pcrtiles = prctile(alpha, [25,50,75]);
    alpha_bayes = [repmat(alpha_pcrtiles(1), 1, reps), alpha_pcrtiles(2), repmat(alpha_pcrtiles(3), 1, reps)]';  % to get the IQR right for the visualisation
    alpha_LS = DoptV./rhoOptV/b^2;
    alphas = [alpha_bayes ; alpha_LS];
    cats = [0*alpha_bayes; 0*alpha_LS+1];
    nexttile([1,2]);
    h = boxchart(cats, alphas);
    xticks([0,1])
    xticklabels({boxlabel1, boxlabel2})

    set(h, 'linewidth', linewidth)
    set(gca, 'ColorOrder', [0,0,0;0,0,0])
    yyaxis left
    ylim([(min(alphas)+max(alphas))/2-0.18, (min(alphas)+max(alphas))/2+0.18])
    yticks(sort(uniquetol([min(alphas(cats==0)), median(alphas(cats==0)), max(alphas(cats==0))],1e-2)))
    ytickformat('%.2f')

    yyaxis right
    ylim([(min(alphas)+max(alphas))/2-0.18, (min(alphas)+max(alphas))/2+0.18])
    yticks(sort(uniquetol([min(alphas(cats==1)), median(alphas(cats==1)), max(alphas(cats==1))],1e-2)))
    ytickformat('%.2f')

    ylabel('$\alpha$', 'fontsize', fontsize, 'interpreter', interpreter)
    xlabel(' ', 'fontsize', fontsize, 'interpreter', interpreter)
    title(' ', 'fontsize',fontsize, 'interpreter', interpreter)
    set(gca, 'fontsize',fontsize, 'ticklabelinterpreter', 'latex')
    box on
    % Add subfigure labels
    if mouse ==1
        title('\textbf{G}', 'interpreter', interpreter, 'fontsize', fontsize,'fontweight','bold')
        ax = gca;
        ax.TitleHorizontalAlignment = 'left';
    elseif mouse==2
        title('\textbf{H}', 'interpreter', interpreter, 'fontsize', fontsize,'fontweight','bold')
        ax = gca;
        ax.TitleHorizontalAlignment = 'left';
    else
        title('\textbf{I}', 'interpreter', interpreter, 'fontsize', fontsize,'fontweight','bold')
        ax = gca;
        ax.TitleHorizontalAlignment = 'left';
    end
    drawnow
end

exportgraphics(gcf, 'Figure6.eps', ContentType='vector')


function [handle,N,Xedges,Yedges] = doplot(C, ax, nbins, nlevels)
mp = @(t) (t(1:end-1) + t(2:end))/2;
dx = (ax(2)-ax(1))/nbins;
dy = (ax(4)-ax(3))/nbins;
Xedges = ax(1):dx:ax(2);
Yedges = ax(3):dy:ax(4);

[X,Y] = meshgrid(mp(Xedges), mp(Yedges));

xi = [X(:),Y(:)];

[N, xi] = ksdensity(C(:,2:3),xi);
N = reshape(N,nbins,nbins)';
[~,handle] = contourf(reshape(xi(:,1),nbins,nbins), reshape(xi(:,2),nbins,nbins), N', nlevels, 'LineStyle', 'none');
axis(ax)
end
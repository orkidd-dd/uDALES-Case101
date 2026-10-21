function figs = plot_tke_publication(S,E,outdir,doExport)
% S(k): label, xH, horizontal, zH, vertical (dimensionless).
% E(k): xH, horizontal, zH, vertical (digitized experimental curves).
% No smoothing, fitted curves or artificial uncertainty bands.
if nargin<4, doExport=true; end
if ~exist(outdir,'dir'), mkdir(outdir); end
% One temperature palette shared by every figure. Source is encoded by marks.
cc=[0.1 0.1 0.1; 0 0.35 0.65; 0.80 0.35 0; 0.45 0.20 0.60];
tempStyles={':','-','--','-.'};
kh=[]; kv=[]; eh=[]; ev=[];
for k=1:4
    kh=[kh; S(k).horizontal(:); E(k).horizontal(:)]; %#ok<AGROW>
    in=S(k).zH>=0 & S(k).zH<=1;
    ie=E(k).zH>=0 & E(k).zH<=1;
    kv=[kv; S(k).vertical(in); E(k).vertical(ie)]; %#ok<AGROW>
    eh=[eh; E(k).horizontal(:)]; %#ok<AGROW>
    v=E(k).vertical(ie); ev=[ev; v(:)]; %#ok<AGROW>
end
hmax=niceMax(kh); vmax=niceMax(kv);
% Experimental overview panels have their own limits. Expand if needed so
% future datasets are not silently clipped. Comparison panels remain shared.
expHmax=max(2,niceMax(eh)); expVmax=max(4,niceMax(ev));
figs=gobjects(1,3);
for direction=1:2
    figs(direction)=figure('Color','w','Units','centimeters', ...
        'Position',[2 2 18 16]);
    tl=tiledlayout(figs(direction),2,2,'TileSpacing','compact','Padding','compact');
    for k=1:4
        a=nexttile(tl); hold(a,'on');
        if direction==1
            is= S(k).zH>=0 & S(k).zH<=1;
            ie= E(k).zH>=0 & E(k).zH<=1;
            he=plot(a,E(k).vertical(ie),E(k).zH(ie),'o','Color',cc(k,:), ...
                'MarkerSize',3.3,'MarkerFaceColor','w','LineWidth',0.65);
            hs=plot(a,S(k).vertical(is),S(k).zH(is),'LineStyle',tempStyles{k}, ...
                'Color',cc(k,:),'LineWidth',1.25);
            xlim(a,[0 vmax]); ylim(a,[0 1]); yticks(a,0:.2:1);
            xlabel(a,'$TKE_{2D}/u_*^2$','Interpreter','latex');
            ylabel(a,'$z/H$','Interpreter','latex');
        else
            he=plot(a,E(k).xH,E(k).horizontal,'o','Color',cc(k,:), ...
                'MarkerSize',3.3,'MarkerFaceColor','w','LineWidth',0.65);
            hs=plot(a,S(k).xH,S(k).horizontal,'LineStyle',tempStyles{k}, ...
                'Color',cc(k,:),'LineWidth',1.25);
            xlim(a,[-.5 .5]); ylim(a,[0 hmax]); xticks(a,-.5:.25:.5);
            xlabel(a,'$x/H$','Interpreter','latex');
            ylabel(a,'$TKE_{2D}/u_*^2$','Interpreter','latex');
        end
        title(a,sprintf('(%c)  \\DeltaT = %s','a'+k-1,S(k).label), ...
            'FontWeight','normal','Interpreter','tex');
        styleAxes(a);
        if k==1
            lg=legend(a,[he hs],{'Experiment (digitized)','Simulation'}, ...
                'Location','best','Box','off'); lg.FontSize=8;
        end
    end
    if direction==1, name='TKE_vertical_comparison'; else, name='TKE_horizontal_comparison'; end
    finish(figs(direction),outdir,name,doExport);
end
% Overview matching the lower sketch: columns=direction, rows=data source.
figs(3)=figure('Color','w','Units','centimeters','Position',[3 3 18 16]);
tl=tiledlayout(figs(3),2,2,'TileSpacing','compact','Padding','compact');
for panel=1:4
    a=nexttile(tl); hold(a,'on');
    for k=1:4
        if panel<=2, d=S(k); else, d=E(k); end
        if mod(panel,2)==1
            ii=d.zH>=0 & d.zH<=1;
            plot(a,d.vertical(ii),d.zH(ii),'Color',cc(k,:), ...
                'LineStyle',tempStyles{k},'LineWidth',1.15,'DisplayName',S(k).label);
        else
            plot(a,d.xH,d.horizontal,'Color',cc(k,:), ...
                'LineStyle',tempStyles{k},'LineWidth',1.15,'DisplayName',S(k).label);
        end
    end
    if mod(panel,2)==1
        if panel<=2
            xlim(a,[0 vmax]);
        else
            xlim(a,[0 expVmax]); xticks(a,linspace(0,expVmax,5));
        end
        ylim(a,[0 1]); yticks(a,0:.2:1);
        xlabel(a,'$TKE_{2D}/u_*^2$','Interpreter','latex'); ylabel(a,'$z/H$','Interpreter','latex');
        direction='Vertical';
    else
        xlim(a,[-.5 .5]); xticks(a,-.5:.25:.5);
        if panel<=2
            ylim(a,[0 hmax]);
        else
            ylim(a,[0 expHmax]); yticks(a,linspace(0,expHmax,5));
        end
        xlabel(a,'$x/H$','Interpreter','latex'); ylabel(a,'$TKE_{2D}/u_*^2$','Interpreter','latex');
        direction='Horizontal';
    end
    if panel<=2, source='Simulation'; else, source='Experiment'; end
    title(a,sprintf('(%c)  %s: %s','a'+panel-1,source,direction),'FontWeight','normal');
    styleAxes(a);
    if panel==1, legend(a,'Location','best','Box','off','FontSize',8); end
end
xlabel(tl,'Experimental panels use different TKE axis limits.', ...
    'FontName','Times New Roman','FontSize',8,'FontWeight','normal');
finish(figs(3),outdir,'TKE_overview',doExport);
end

function styleAxes(a)
set(a,'FontName','Times New Roman','FontSize',10,'LineWidth',0.75, ...
    'TickDir','in','TickLength',[.018 .018],'Box','on','XGrid','off','YGrid','off');
a.XLabel.FontSize=11; a.YLabel.FontSize=11; a.Title.FontSize=10;
end

function lim=niceMax(v)
v=v(isfinite(v)); assert(~isempty(v),'No finite TKE data.');
m=max(v)*1.06;
if m<=0, lim=1; return; end
step=10^floor(log10(m))/5; lim=ceil(m/step)*step;
end

function finish(f,outdir,name,doExport)
drawnow;
if doExport
    exportgraphics(f,fullfile(outdir,[name '.pdf']),'ContentType','vector');
    exportgraphics(f,fullfile(outdir,[name '.png']),'Resolution',600);
    savefig(f,fullfile(outdir,[name '.fig']));
end
end

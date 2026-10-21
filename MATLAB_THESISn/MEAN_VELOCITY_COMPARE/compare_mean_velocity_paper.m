%% Four-temperature mean-velocity comparison (single-file workflow)
% Experiment V is vertical velocity, called W here to match uDALES w.
% All local functions are included below. udbase is not required.
clear; clc;
here=fileparts(mfilename('fullpath'));
outdir=fullfile(here,'mean_velocity_figures');
if ~exist(outdir,'dir'), mkdir(outdir); end

%% Settings
root='H:\udales_data\u-dales-output';
folders={'101_128','101_70k','101_170k','101_240k'};
temps=[0 70 170 240];
H=16; xc=32; yc=32; zmid=H/2;
ustar=[.14 .14 .14 .14];
timeWindow=[600 3500];
recomputeSimulation=true; % false: replot from the saved cache
doExport=true;
cache=fullfile(outdir,'mean_velocity_profiles.mat');
settings=struct('root',root,'folders',{folders},'temps',temps,'H',H, ...
    'xc',xc,'yc',yc,'zmid',zmid,'ustar',ustar,'timeWindow',timeWindow);
% Exactly the same palette and line styles as the TKE figures.
colors=[.1 .1 .1; 0 .35 .65; .80 .35 0; .45 .20 .60];
styles={':','-','--','-.'};

%% Experimental curves (already normalized by experimental u*)
E=struct([]);
for k=1:4
    u=readtable(fullfile(here,sprintf('U_vertical_%dK_nodes.csv',temps(k))));
    w=readtable(fullfile(here,sprintf('V_horizontal_%dK_nodes.csv',temps(k))));
    E(k).z=u.z_over_H; E(k).u=u.U_over_ustar;
    E(k).x=w.x_over_H; E(k).w=w.V_over_ustar;
end

%% Simulation: use each component's native staggered grid
if recomputeSimulation
    assert(all(isfinite(ustar)&ustar>0),'Friction velocity must be positive.');
    S=struct([]);
    logid=fopen(fullfile(outdir,'sampling_log.txt'),'w');
    assert(logid>=0,'Cannot create sampling log.');
    logCleanup=onCleanup(@() fclose(logid));
    for k=1:4
        file=fullfile(root,folders{k},'tdump.101.nc');
        t=readCoordinate(file,'time');
        ti=find(t>=timeWindow(1)&t<=timeWindow(2));
        assert(~isempty(ti)&&all(diff(ti)==1),'Empty/noncontiguous time selection.');
        % ut is on (xm,yt,zt); wt is on (xt,yt,zm).
        [U,xu,zu,uy]=meanYPlane(file,'ut','xm','zt',yc,ti);
        [W,xw,zw,wy]=meanYPlane(file,'wt','xt','zm',yc,ti);
        uv=interp1(xu,U,xc,'linear',NaN).';
        wh=interp1(zw,W.',zmid,'linear',NaN).';
        assert(any(isfinite(uv))&&any(isfinite(wh)),'Profile outside simulation domain.');
        iz=zu>=0&zu<=H; ix=xw>xc-H/2&xw<xc+H/2;
        S(k).label=sprintf('%d K',temps(k));
        S(k).z=zu(iz)/H; S(k).u=uv(iz)/ustar(k);
        S(k).x=(xw(ix)-xc)/H; S(k).w=wh(ix)/ustar(k);
        S(k).records=numel(ti); S(k).timeRange=t(ti([1 end]));
        % No artificial no-slip points or spatial extrapolation are added.
        writetable(table(S(k).z,S(k).u,'VariableNames',{'z_over_H','U_over_ustar'}), ...
            fullfile(outdir,sprintf('simulation_U_%dK.csv',temps(k))));
        writetable(table(S(k).x,S(k).w,'VariableNames',{'x_over_H','W_over_ustar'}), ...
            fullfile(outdir,sprintf('simulation_W_%dK.csv',temps(k))));
        msg=sprintf(['%d K: %d records, %.6f..%.6f s; ' ...
            'U: x=%.2f m, y=[%.2f %.2f] -> %.2f m; ' ...
            'W: z=%.2f m, y=[%.2f %.2f] -> %.2f m.\n'], ...
            temps(k),numel(ti),t(ti(1)),t(ti(end)),xc,uy(1),uy(2),yc,zmid,wy(1),wy(2),yc);
        fprintf('%s',msg); fprintf(logid,'%s',msg);
    end
    clear logCleanup;
    save(cache,'S','settings');
else
    C=load(cache,'S','settings');
    assert(isequaln(C.settings,settings),'Cache settings differ; recompute first.');
    S=C.S;
end

%% Matched experimental nodes: interpolation, no extrapolation
rows=cell(8,5); r=0;
for k=1:4
    for dir=1:2
        if dir==1
            q=E(k).z; e=E(k).u; s=interp1(S(k).z,S(k).u,q,'linear',NaN);
            ok=q>=0&q<=1&isfinite(e)&isfinite(s);
            physical=H*q; coord='z_over_H'; kind='U_vertical';
        else
            q=E(k).x; e=E(k).w; s=interp1(S(k).x,S(k).w,q,'linear',NaN);
            ok=abs(q)<=.5&isfinite(e)&isfinite(s);
            physical=xc+H*q; coord='x_over_H'; kind='W_horizontal';
        end
        writetable(table(q(ok),physical(ok),e(ok),s(ok), ...
            'VariableNames',{coord,'coordinate_m','Experiment','Simulation'}), ...
            fullfile(outdir,sprintf('matched_%s_%dK.csv',kind,temps(k))));
        err=s(ok)-e(ok); r=r+1;
        if isempty(err), rmse=NaN; mae=NaN;
        else, rmse=sqrt(mean(err.^2)); mae=mean(abs(err)); end
        rows(r,:)={temps(k),kind,sum(ok),rmse,mae};
    end
end
metrics=cell2table(rows,'VariableNames',{'DeltaT_K','Profile','N','RMSE','MAE'});
writetable(metrics,fullfile(outdir,'comparison_metrics.csv'));
disp(metrics);

%% Plot three figures
uSim=[]; uExp=[]; wSim=[]; wExp=[];
for k=1:4
    uSim=[uSim;S(k).u(:)]; wSim=[wSim;S(k).w(:)]; %#ok<AGROW>
    iz=E(k).z>=0&E(k).z<=1;
    uExp=[uExp;E(k).u(iz)]; wExp=[wExp;E(k).w(:)]; %#ok<AGROW>
end
uCommon=velocityLimit([uSim;uExp]); wCommon=velocityLimit([wSim;wExp]);
figs=gobjects(1,3);
for dir=1:2
    figs(dir)=figure('Color','w','Units','centimeters','Position',[2 2 18 16]);
    tl=tiledlayout(figs(dir),2,2,'TileSpacing','compact','Padding','compact');
    for k=1:4
        a=nexttile(tl); hold(a,'on');
        if dir==1
            ie=E(k).z>=0&E(k).z<=1;
            he=plot(a,E(k).u(ie),E(k).z(ie),'o','Color',colors(k,:), ...
                'MarkerFaceColor','w','MarkerSize',3.3,'LineWidth',.7);
            hs=plot(a,S(k).u,S(k).z,'Color',colors(k,:), ...
                'LineStyle',styles{k},'LineWidth',1.3);
            setupProfile(a,1,uCommon);
        else
            he=plot(a,E(k).x,E(k).w,'o','Color',colors(k,:), ...
                'MarkerFaceColor','w','MarkerSize',3.3,'LineWidth',.7);
            hs=plot(a,S(k).x,S(k).w,'Color',colors(k,:), ...
                'LineStyle',styles{k},'LineWidth',1.3);
            setupProfile(a,2,wCommon);
        end
        title(a,sprintf('(%c)  \\DeltaT = %d K','a'+k-1,temps(k)), ...
            'FontWeight','normal','Interpreter','tex');
        if k==1
            legend(a,[he hs],{'Experiment (digitized)','Simulation'}, ...
                'Location','best','Box','off','FontSize',8);
        end
    end
    if dir==1, name='U_vertical_comparison'; else, name='W_horizontal_comparison'; end
    exportFigure(figs(dir),outdir,name,doExport);
end

figs(3)=figure('Color','w','Units','centimeters','Position',[3 3 18 16]);
tl=tiledlayout(figs(3),2,2,'TileSpacing','compact','Padding','compact');
limits=[velocityLimit(uSim),velocityLimit(wSim), ...
    max(5.5,velocityLimit(uExp)),max(5,velocityLimit(wExp))];
for panel=1:4
    a=nexttile(tl); hold(a,'on');
    for k=1:4
        if panel<=2, d=S(k); else, d=E(k); end
        if mod(panel,2)==1
            iz=d.z>=0&d.z<=1; x=d.u(iz); y=d.z(iz);
        else
            x=d.x; y=d.w;
        end
        plot(a,x,y,'Color',colors(k,:),'LineStyle',styles{k}, ...
            'LineWidth',1.25,'DisplayName',sprintf('%d K',temps(k)));
    end
    setupProfile(a,1+mod(panel+1,2),limits(panel));
    if panel<=2, source='Simulation'; else, source='Experiment'; end
    if mod(panel,2)==1, name='Horizontal velocity'; else, name='Vertical velocity'; end
    title(a,sprintf('(%c)  %s: %s','a'+panel-1,source,name),'FontWeight','normal');
    if panel==1, legend(a,'Location','best','Box','off','FontSize',8); end
end
xlabel(tl,'Velocity axis limits may differ between simulation and experiment.', ...
    'FontName','Times New Roman','FontSize',8);
exportFigure(figs(3),outdir,'mean_velocity_overview',doExport);
fprintf('Saved figures and matched profiles to:\n%s\n',outdir);

%% Local functions
function q=readCoordinate(file,name)
q=double(ncread(file,name)); q=q(:);
assert(all(isfinite(q))&&all(diff(q)>0),'Non-increasing coordinate: %s.',name);
end

function [F,x,z,ybracket]=meanYPlane(file,var,xname,zname,yc,ti)
x=readCoordinate(file,xname); y=readCoordinate(file,'yt'); z=readCoordinate(file,zname);
j=find(y<=yc,1,'last');
assert(~isempty(j)&&j<numel(y),'Cannot bracket target y.');
w=(yc-y(j))/(y(j+1)-y(j)); ybracket=y([j j+1]);
info=ncinfo(file,var); names={info.Dimensions.Name};
target={xname,'yt',zname,'time'}; order=zeros(1,4);
assert(numel(names)==4,'Expected four-dimensional velocity data.');
start=ones(1,4); count=info.Size;
for d=1:4
    id=find(strcmp(names,target{d})); assert(isscalar(id),'Unexpected staggered grid for %s.',var);
    order(d)=id;
end
assert(ti(end)<=count(order(4)),'Time dimension mismatch.');
start(order(2))=j; count(order(2))=2;
start(order(4))=ti(1); count(order(4))=numel(ti);
A=double(ncread(file,var,start,count)); A=permute(A,order);
assert(size(A,1)==numel(x)&&size(A,2)==2&&size(A,3)==numel(z));
% Interpolate each record to the same y plane before time averaging.
P=(1-w)*A(:,1,:,:)+w*A(:,2,:,:);
F=reshape(mean(P,4,'omitnan'),numel(x),numel(z));
end

function lim=velocityLimit(v)
v=v(isfinite(v)); assert(~isempty(v),'No finite velocity data.');
lim=max(1,ceil(1.08*max(abs(v))*2)/2);
end

function setupProfile(a,dir,lim)
if dir==1
    xlim(a,[-lim lim]); ylim(a,[0 1]); yticks(a,0:.2:1);
    xlabel(a,'$U/u_*$','Interpreter','latex'); ylabel(a,'$z/H$','Interpreter','latex');
    xline(a,0,':','Color',[.75 .75 .75],'HandleVisibility','off','LineWidth',.5);
else
    xlim(a,[-.5 .5]); ylim(a,[-lim lim]); xticks(a,-.5:.25:.5);
    xlabel(a,'$x/H$','Interpreter','latex'); ylabel(a,'$W/u_*$','Interpreter','latex');
    yline(a,0,':','Color',[.75 .75 .75],'HandleVisibility','off','LineWidth',.5);
end
set(a,'FontName','Times New Roman','FontSize',10,'LineWidth',.75, ...
    'TickDir','in','Box','on','XGrid','off','YGrid','off');
a.XLabel.FontSize=11; a.YLabel.FontSize=11; a.Title.FontSize=10;
end

function exportFigure(f,outdir,name,enabled)
drawnow;
if enabled
    exportgraphics(f,fullfile(outdir,[name '.pdf']),'ContentType','vector');
    exportgraphics(f,fullfile(outdir,[name '.png']),'Resolution',600);
    savefig(f,fullfile(outdir,[name '.fig']));
end
end

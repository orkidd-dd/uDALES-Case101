function plot_temperature_publication(here,showFigure,mode)
% Compare experimental and LES mean temperature fields for DW heating, H/W = 1.
% mode = 'equivalent' (default): LES shown as equivalent experimental temperature
%        T_eq = T0_exp + theta* x dT_exp, using the same isotherm values as the experiment.
% mode = 'normalized': LES shown as theta* = (theta - theta_ref)/dT_LES.
% Isotherms are drawn in two tiers: thin, dense lines (~0.01 in theta*) to match the
% line density of the experimental ink traces, and thicker labelled lines at round values.
if nargin<1 || isempty(here), here=fileparts(mfilename('fullpath')); end
if nargin<2 || isempty(showFigure), showFigure=true; end
if nargin<3 || isempty(mode), mode='equivalent'; end
mode=validatestring(mode,{'equivalent','normalized'});

% ---------------- User settings ----------------
Href = 16;               % building height (m)
xc   = 32;               % canyon centre (m)
W    = 16;               % canyon width (m), H/W = 1
thetaRefC = 295-273.15;  % LES reference potential temperature (deg C)
T0expC    = 300-273.15;  % approx. experimental reference temperature, T0 ~ 300 K (Fellini et al., 2020)
cases = [70 170 240];    % experimental Delta T (K)
% Labelled (thick) isotherm values (deg C) for the equivalent-temperature mode,
% one cell per case: every 3 deg C for 70 K and 170 K, every 5 deg C for 240 K.
% Values outside the data range are dropped automatically.
eqLevels   = {30:3:60, 30:3:90, 30:5:100};
starLevels = 0:0.05:0.6; % labelled (thick) isotherm values for the normalized mode
starStep   = 0.01;       % target spacing of the thin isotherms, in theta* (both modes)
% Label placement: every labelled isotherm gets one label where the neighbouring
% thin lines are at least labelMinGap (in units of H) apart, so labels do not
% cover the packed lines next to the heated wall.
labelMinGap = 0.03;
labelMax    = 20;        % max labels per panel
% Extra dashed isotherms in the well-mixed vortex core (where the thin lines are
% far apart). They subdivide the thin-line interval and are drawn only there.
coreFine   = false;      % set true to switch them on
coreSub    = 5;          % fine (dashed) lines per thin-line interval
coreGapMin = 0.08;       % core = where thin lines are at least this far apart (units of H)
fineColor=[.45 .45 .45]; fineWidth=.3;
% Light grey shading under the isotherms. The grey scale is stretched over the
% temperature range of the core (same core definition as above), so small
% differences inside the vortex become visible; values outside saturate.
shadeCore  = true;            % set false for plain white panels
shadeQuant = [.02 .98];       % core quantiles mapped to the lightest/darkest grey
shadeGrey  = [1 .70];         % grey level for cooler -> warmer (1 = white)
shadeUp    = 4;               % upsampling of the LES grid for a smooth shade
% Line styles: thin dense isotherms and thick labelled isotherms
minorColor=[.35 .35 .35]; minorWidth=.3;
majorColor=[.1 .1 .1];    majorWidth=.6;
% Publication layout (mm and pt)
figureWidth=183; figureHeight=263; panelWidth=74;
axisFont=9; titleFont=10; contourFont=9;
xs=[14 103]; ys=[180 96 12];
% -----------------------------------------------

out=fullfile(here,'figures'); if ~exist(out,'dir'), mkdir(out); end
D=load(fullfile(here,'data','experimental_ink_traces.mat')); E=D.E;
sources=readtable(fullfile(fileparts(here),'Figure11_all_updated','Figure11_sources_and_differences.csv'));
sources=sources(strcmp(sources.Wall,'DW'),:);
assert(isequal(sources.ExperimentalLabel_K.',cases));

fields=readTemperatureFields(sources,xc,W,Href);

% Normalized and equivalent temperatures, plus canyon-averaged warming for the text
meanTstar=zeros(numel(cases),1); coreTstar=meanTstar;
for j=1:numel(fields)
    dTles=sources.LES_deltaT_K(j);
    fields(j).Tstar=(fields(j).T-thetaRefC)/dTles;
    fields(j).Teq=T0expC+fields(j).Tstar*cases(j);
    meanTstar(j)=mean(fields(j).Tstar,'all');
    inCore=(fields(j).z>.08 & fields(j).z<.92) & (abs(fields(j).x)<.38).';
    coreTstar(j)=mean(fields(j).Tstar(inCore));
    fprintf('%d K: canyon-mean theta* = %.3f, core-mean theta* = %.3f\n',cases(j),meanTstar(j),coreTstar(j));
end
sources.MinTheta_C=arrayfun(@(c)min(c.T,[],'all'),fields).';
sources.MaxTheta_C=arrayfun(@(c)max(c.T,[],'all'),fields).';
sources.CanyonMeanThetaStar=meanTstar;
sources.CoreMeanThetaStar=coreTstar;
writetable(sources(:,{'ExperimentalLabel_K','LES_deltaT_K','Records','Start_s','End_s', ...
    'MinTheta_C','MaxTheta_C','CanyonMeanThetaStar','CoreMeanThetaStar','Source'}), ...
    fullfile(out,'temperature_sources.csv'));

f=figure('Visible','off','Color','w','Units','centimeters', ...
    'Position',[2 2 figureWidth/10 figureHeight/10]);
panels={};
for row=1:3
    for col=1:2
        n=(row-1)*2+col;
        ax=axes(f,'Units','centimeters','Position',[xs(col)/10 ys(row)/10 panelWidth/10 panelWidth/10]);
        hold(ax,'on');
        if col==1
            e=E(row); rgb=repmat(1-e.ink,1,1,3);
            image(ax,[e.x(1) e.x(end)],[e.z(1) e.z(end)],rgb); set(ax,'YDir','normal');
            kind='Experiment';
        else
            c=fields(row);
            switch mode
                case 'equivalent'
                    Z=c.Teq; major=eqLevels{row}; fmt='%g';
                    % Subdivide the labelled step so that the thin lines are ~starStep apart
                    % in theta* and every labelled line lies on the same grid.
                    dMajor=median(diff(major));
                    nSub=max(1,round(dMajor/(starStep*cases(row))));
                    dMinor=dMajor/nSub; base=major(1);
                case 'normalized'
                    Z=c.Tstar; major=starLevels; fmt='%.2f';
                    dMinor=starStep; base=0;
            end
            zmin=min(Z(:)); zmax=max(Z(:));
            major=major(major>=zmin & major<=zmax);
            minor=base+(ceil((zmin-base)/dMinor):floor((zmax-base)/dMinor))*dMinor;
            tol=1e-6*dMinor;
            minor=minor(all(abs(minor(:)-major(:).')>tol,2));   % no doubled lines
            if numel(major)<2
                warning('%d K: fewer than two labelled isotherms fall in the data range; check the level settings.',cases(row));
            end
            fprintf('%d K: thin isotherm spacing = %.3g (theta* = %.4f), %d thin + %d labelled lines\n', ...
                cases(row),dMinor,dMinor/(strcmp(mode,'equivalent')*cases(row)+strcmp(mode,'normalized')), ...
                numel(minor),numel(major));
            % Local distance between neighbouring thin isotherms (units of H)
            [gx,gz]=gradient(Z,c.x,c.z);
            gap=dMinor./max(hypot(gx,gz),eps);
            if shadeCore
                v=Z(gap>=coreGapMin);
                if numel(v)<10, v=Z(:); end          % no clear core: use the whole field
                v=sort(v(:)); qi=@(p)v(max(1,min(numel(v),round(p*(numel(v)-1))+1)));
                cl=[qi(shadeQuant(1)) qi(shadeQuant(2))];
                if cl(2)<=cl(1), cl=cl(1)+[-1 1]*max(eps(cl(1)),1e-6); end
                xf=linspace(c.x(1),c.x(end),shadeUp*numel(c.x));
                zf=linspace(c.z(1),c.z(end),shadeUp*numel(c.z));
                Zf=interp2(c.x,c.z,Z,xf,zf.','linear');
                imagesc(ax,xf,zf,Zf); set(ax,'YDir','normal');
                g=linspace(shadeGrey(1),shadeGrey(2),256).';
                colormap(ax,[g g g]); set(ax,'CLim',cl);
                fprintf('   grey shade: %.3g (white) to %.3g (darkest)\n',cl(1),cl(2));
            end
            if coreFine
                dFine=dMinor/coreSub;
                k=ceil((zmin-base)/dFine):floor((zmax-base)/dFine);
                k=k(mod(k,coreSub)~=0);              % skip values already drawn as thin lines
                fine=base+k*dFine;
                Zc=Z; Zc(gap<coreGapMin)=NaN;       % restrict to the quiet core
                if ~isempty(fine) && any(isfinite(Zc(:)))
                    contour(ax,c.x,c.z,Zc,levelArg(fine),'LineStyle','--', ...
                        'Color',fineColor,'LineWidth',fineWidth);
                end
                fprintf('   dashed core isotherms every %.3g (core = %.0f%% of panel)\n', ...
                    dFine,100*mean(isfinite(Zc(:))));
            end
            if ~isempty(minor)
                contour(ax,c.x,c.z,Z,levelArg(minor),'Color',minorColor,'LineWidth',minorWidth);
            end
            if ~isempty(major)
                C=contour(ax,c.x,c.z,Z,levelArg(major),'Color',majorColor,'LineWidth',majorWidth);
                placeLabels(ax,C,major,fmt,contourFont,c.x,c.z,gap,labelMinGap,labelMax);
            end
            plot(ax,[-.5 -.5 .5 .5],[1 0 0 1],'Color',[.45 .45 .45],'LineWidth',1.3);
            kind='Simulation';
        end
        axis(ax,'equal'); xlim(ax,[-.5 .5]); ylim(ax,[0 1]);
        xticks(ax,[-.5 0 .5]); yticks(ax,[0 .25 .5 .75 1]);
        set(ax,'FontName','Times New Roman','FontSize',axisFont,'FontWeight','bold', ...
            'Box','on','TickDir','in','LineWidth',.6,'Layer','top');
        xlabel(ax,'x/H','FontSize',10,'FontWeight','bold');
        ylabel(ax,'z/H','FontSize',10,'FontWeight','bold');
        title(ax,sprintf('(%c) %s, \\DeltaT_{exp} = %d K','a'+n-1,kind,cases(row)), ...
            'FontSize',titleFont,'FontWeight','bold');
        ax.XTickLabel={}; ax.YTickLabel={};
        for v=[-.5 0 .5]
            text(ax,v,-.025,sprintf('%g',v),'FontName','Times New Roman','FontSize',axisFont, ...
                'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','top','Clipping','off');
        end
        for v=[0 .25 .5 .75 1]
            text(ax,-.52,v,sprintf('%g',v),'FontName','Times New Roman','FontSize',axisFont, ...
                'FontWeight','bold','HorizontalAlignment','right','VerticalAlignment','middle','Clipping','off');
        end
        ax.XLabel.Position=[0 -.075 0]; ax.YLabel.Position=[-.625 .5 0];
        if row<3, ax.XLabel.String=''; end
        drawnow; ax.Units='points'; b=ax.Position;
        panels{end+1}=struct('id',char('a'+n-1),'bbox_pt',[b(1:2) b(1:2)+b(3:4)], ...
            'grid_id','main','row_start',row-1,'row_stop',row,'col_start',col-1,'col_stop',col); %#ok<AGROW>
    end
end


drawnow; f.Units='points'; sz=f.Position(3:4);
L=struct('schema_version',1,'backend','matlab', ...
    'figure',struct('width_pt',sz(1),'height_pt',sz(2)),'panels',{panels});
fid=fopen(fullfile(out,'temperature.layout.json'),'w'); fprintf(fid,'%s',jsonencode(L)); fclose(fid);

base=fullfile(out,['HW1_DW_temperature_' mode]);
exportgraphics(f,[base '.png'],'Resolution',600,'Padding','figure');
exportgraphics(f,[base '.pdf'],'ContentType','vector','Padding','figure');
savefig(f,[base '.fig']);
if showFigure, set(f,'Visible','on'); else, close(f); end
end

% -------------------------------------------------------------------------
function v=levelArg(v)
% contour() treats a scalar as "number of levels"; repeat it to request one level.
if isscalar(v), v=[v v]; end
end

% -------------------------------------------------------------------------
function fields=readTemperatureFields(sources,xc,W,Href)
% Time- and spanwise-averaged potential temperature in the canyon (deg C).
for j=1:height(sources)
    file=fullfile(sources.Source{j},'tdump.101.nc');
    t=double(ncread(file,'time')); t=t(:);
    assert(numel(t)>=2,'Need at least two records to infer the averaging interval.');
    % Start time of each averaging record. The first record's start is inferred
    % from the dump interval, which also handles restarted runs.
    left=[t(1)-(t(2)-t(1)); t(1:end-1)];
    selected=find(left>=sources.Start_s(j)-1e-8 & t<=sources.End_s(j)+1e-8);
    assert(numel(selected)==sources.Records(j) && all(diff(selected)==1), ...
        'Selected records do not match the expected averaging window.');

    x=double(ncread(file,'xt')); y=double(ncread(file,'yt')); z=double(ncread(file,'zt'));
    ix=find(x>xc-W/2 & x<xc+W/2);
    iz=find(z>0 & z<Href);
    ny=numel(y);
    weight=(t(selected)-left(selected))/(t(selected(end))-left(selected(1)));

    info=ncinfo(file,'thlt');
    dims={info.Dimensions.Name};
    order=cellfun(@(name)find(strcmp(dims,name)),{'xt','yt','zt','time'});
    start=ones(1,4); count=start;
    start(order(1:3))=[ix(1) 1 iz(1)];
    count(order(1:3))=[numel(ix) ny numel(iz)];

    T=zeros(numel(ix),ny,numel(iz));
    for k=1:numel(selected)
        start(order(4))=selected(k);
        slab=permute(double(ncread(file,'thlt',start,count)),order);
        assert(all(isfinite(slab),'all') && all(slab>200 & slab<400,'all'));
        T=T+weight(k)*slab;
    end
    plane=squeeze(mean(T,2));            % spanwise average (2D canyon)
    fields(j).x=(x(ix)-xc)/Href; %#ok<AGROW>
    fields(j).z=z(iz)/Href;
    fields(j).T=plane.'-273.15;          % [nz x nx], deg C (potential temperature)
end
end

% -------------------------------------------------------------------------
function placeLabels(ax,C,labelVals,fmt,fontSize,gx,gz,gap,minGap,maxLabels)
% One label per isotherm, placed where neighbouring lines are not packed too
% tightly and away from labels already placed.
placed=[]; num=0;
for v=labelVals
    ptr=1; p=[];
    while ptr<size(C,2)
        np=C(2,ptr);
        if abs(C(1,ptr)-v)<1e-7*max(1,abs(v)), p=[p C(:,ptr+1:ptr+np)]; end %#ok<AGROW>
        ptr=ptr+np+1;
    end
    if isempty(p), continue; end
    inside=p(1,:)>-.45 & p(1,:)<.45 & p(2,:)>.04 & p(2,:)<.96;
    g=interp2(gx,gz,gap,p(1,:),p(2,:));
    valid=find(inside & g>=minGap);
    if isempty(valid)
        % Fall back to the least crowded point on this line, if it is not too tight
        cand=find(inside & g>=.5*minGap);
        if isempty(cand), continue; end
        [~,k]=max(g(cand)); valid=cand(k);
    end
    % Prefer roomy spots, spread along a diagonal so labels do not stack up
    target=[-.3+.2*mod(num,4); .85-.2*mod(num,4)];
    cost=sum((p(:,valid)-target).^2,1)-.02*min(g(valid)/minGap,3);
    [~,ord]=sort(cost);
    for ii=valid(ord)
        q=p(:,ii).';
        if ~isempty(placed) && any(sum(((placed-q)./[.13 .075]).^2,2)<1), continue; end
        text(ax,q(1),q(2),sprintf(fmt,v),'FontName','Times New Roman','FontSize',fontSize, ...
            'FontWeight','bold','HorizontalAlignment','center','BackgroundColor','w', ...
            'Margin',.6,'Clipping','on');
        placed=[placed; q]; num=num+1; %#ok<AGROW>
        break;
    end
    if num>=maxLabels, break; end
end
end

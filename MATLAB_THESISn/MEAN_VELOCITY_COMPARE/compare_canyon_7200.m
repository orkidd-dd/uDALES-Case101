%% Single-file MATLAB workflow: new isothermal baseline versus digitized WT data
% Adapted from D:/Desktop/engineering/FLUILDCITY/code/test.m.
% All helpers are at the end of this file; no udbase or other .m file required.
% Stored variances are interval statistics, not reconstructed pooled variance.
clear;clc;close all;
%% EDIT SETTINGS HERE
caseDir='D:/udales_data/u-dales-output/101_HW1_isothermal295_7200s';
timeWindow=[1000 7200]; H=16; xc=32; zHorizontal=8;
uStarLES=0.0586; uStarExperiment=0.22; % Experiment CSVs already normalized.
figureSizeCm=[16 13.85]; exportDPI=600; % Same physical size for both thesis figures.
spatialMode='ymean'; % 'ymean' or 'y32' (interpolated y=32 m)
sourceDir='D:/Desktop/engineering/FLUILDCITY/code';
velocityDir='E:/codex/migration-2026-10-28/D-Codex/2026-10-23/upwind_mean_velocity_0_140_240/data';
outDir=fullfile(fileparts(mfilename('fullpath')),'figures');
if ~isfolder(outDir),mkdir(outDir);end
file=fullfile(caseDir,'tdump.101.nc');assert(isfile(file));
%% READ COORDINATES AND SELECT STATISTICS
t=double(ncread(file,'time'));t=t(:);left=[0;t(1:end-1)];
assert(all(diff(t)>0));actualEnd=min(timeWindow(2),t(end));
weights=max(0,min(t,actualEnd)-max(left,timeWindow(1))); assert(sum(weights)>0);
xt=double(ncread(file,'xt'));xm=double(ncread(file,'xm'));
yt=double(ncread(file,'yt'));zt=double(ncread(file,'zt'));zm=double(ncread(file,'zm'));
xt=xt(:);xm=xm(:);yt=yt(:);zt=zt(:);zm=zm(:);
ix=find(xt>xc-H/2 & xt<xc+H/2);
im=find(xm>=xc-H/2 & xm<=xc+H/2);
iz=find(zt<=H+.5);iw=find(zm<=H+.5);
x=xt(ix);xu=xm(im);z=zt(iz);zw=zm(iw);
fprintf('Reading %d interval records; requested %.0f-%.0f s; actual end %.6f s\n',nnz(weights),timeWindow,t(end));
fprintf('Spatial averaging: %s; candidate LES u*=%.6f m/s\n',spatialMode,uStarLES);
U=zeros(numel(im),numel(iz));W=zeros(numel(ix),numel(iw));K=zeros(numel(ix),numel(iz));
for it=find(weights>0).'
 u=readSlab(file,'ut','xm','zt',im,iz,it,yt,spatialMode);
 w=readSlab(file,'wt','xt','zm',ix,iw,it,yt,spatialMode);
 uu=readSlab(file,'upuptc','xt','zt',ix,iz,it,yt,spatialMode);
 ww=readSlab(file,'wpwptc','xt','zt',ix,iz,it,yt,spatialMode);
 assert(all(isfinite(u),'all') && all(isfinite(w),'all') && all(isfinite(uu),'all') && all(isfinite(ww),'all'),'Nonfinite values in canyon slab');
 assert(min(uu,[],'all')>-1e-6 && min(ww,[],'all')>-1e-6,'Negative stored variance');
 q=weights(it)/sum(weights);U=U+q*u;W=W+q*w;K=K+q*.5*(uu+ww);
end
%% NATIVE-GRID INTERPOLATION: U AND W ARE ON DIFFERENT GRIDS
uv=interp1(xu,U,xc,'linear').';
wv=interp1(zw,W.',zHorizontal,'linear').';
kv=interp1(x,K,xc,'linear').';
kh=interp1(z,K.',zHorizontal,'linear').';
% Restrict plots to interior samples; no wall extrapolation or artificial zeros.
zUse=z/H<=1;xUse=abs((x-xc)/H)<.5;
%% EXPERIMENTS: KEEP USER-SUPPLIED DIGITIZED VALUES UNCHANGED
eh=readmatrix(fullfile(sourceDir,'horizontal_TKE.csv'));
ev=readmatrix(fullfile(sourceDir,'Vertical_TKE.csv'));
eHx=eh(:,3);eHk=eh(:,4);eVk=ev(:,3);eVz=ev(:,4);
valid=isfinite(eHx)&isfinite(eHk);eHx=eHx(valid);eHk=eHk(valid);
valid=isfinite(eVz)&isfinite(eVk);eVz=eVz(valid);eVk=eVk(valid);
eU=readmatrix(fullfile(velocityDir,'experiment_vertical_0K.csv'));
eW=readmatrix(fullfile(velocityDir,'experiment_horizontal_0K.csv'));
fprintf('Peak LES k2D: vertical %.6f, horizontal %.6f m2/s2\n',max(kv(zUse)),max(kh(xUse)));
%% TKE: TWO PANELS WITHOUT AMBIGUOUS OVERLAID AXES
f=figure('Visible','off','Color','w','Units','centimeters','Position',[2 2 24 14]);
a=axes(f,'Position',[.10 .24 .35 .62]);b=axes(f,'Position',[.60 .24 .35 .62]);hold(a,'on');hold(b,'on');
plot(a,kv(zUse)/uStarLES^2,z(zUse)/H,'k-','LineWidth',1.6);
plot(a,eVk,eVz,'ko','MarkerFaceColor','k','MarkerSize',3);
plot(b,(x(xUse)-xc)/H,kh(xUse)/uStarLES^2,'r--','LineWidth',1.6);
plot(b,eHx,eHk,'r^','MarkerFaceColor','r','MarkerSize',3);
ylim(a,[0 1]);xlim(b,[-.5 .5]);fmt(a);fmt(b);
xlabel(a,'k2D / (u*)^2','Interpreter','none');ylabel(a,'z/H');
xlabel(b,'(x-xc)/H');ylabel(b,'k2D / (u*)^2','Interpreter','none');
title(a,'(a) Vertical profile at x = 32 m');title(b,'(b) Horizontal profile at z = 8 m');
legend(a,{'uDALES','Digitized experiment'},'Location','best','Box','off');legend(b,{'uDALES','Digitized experiment'},'Location','best','Box','off');
footer(f,timeWindow,actualEnd,uStarLES,spatialMode);emit(f,outDir,'TKE_two_panel');
%% TKE COMPOSITE: MATCH THE THESIS VELOCITY FIGURE SIZE AND TYPOGRAPHY
% No zero values are appended to TKE; all LES values remain extracted data.
f=figure('Visible','off','Color','w','Units','centimeters','Position',[2 2 figureSizeCm]);
a=axes(f,'Position',[.12 .13 .76 .75]);hold(a,'on');
p1=plot(a,(x(xUse)-xc)/H,kh(xUse)/uStarLES^2,'r--','LineWidth',1.35);
p2=plot(a,eHx,eHk,'r^','MarkerFaceColor','r','MarkerSize',3.5);
set(a,'XAxisLocation','top','YAxisLocation','left','XColor','k','YColor','k');
xlim(a,[-.5 .5]);xticks(a,-.5:.1:.5);ylim(a,[0 4.2]);yticks(a,0:.5:4);
xlabel(a,'x/H','FontWeight','normal','FontSize',11,'Color','r');
ylabel(a,'TKE / u*²','Interpreter','none','FontWeight','bold','FontSize',12);grid(a,'on');
b=axes(f,'Position',a.Position,'Color','none','XAxisLocation','bottom','YAxisLocation','right');hold(b,'on');
p3=plot(b,kv(zUse)/uStarLES^2,z(zUse)/H,'k-','LineWidth',1.35);
p4=plot(b,eVk,eVz,'ko','MarkerFaceColor','k','MarkerSize',3);
xlim(b,[0 13]);xticks(b,0:2:12);ylim(b,[0 1.04]);yticks(b,0:.1:1);
xlabel(b,'TKE / u*²','Interpreter','none','FontWeight','bold','FontSize',12);
ylabel(b,'z/H','FontWeight','normal','FontSize',11);grid(b,'off');
set([a b],'FontName','Arial','FontSize',9,'LineWidth',.65,'GridAlpha',.14,'GridColor',[.4 .4 .4],'TickDir','in');
lgH=legend(a,[p1 p2],{'TKE Horiz (Sim)','TKE Horiz (Exp)'},'FontSize',7,'Box','on');
lgV=legend(b,[p3 p4],{'TKE Vert (Sim)','TKE Vert (Exp)'},'FontSize',7,'Box','on');
% Separate legend boxes, positioned in empty space without covering data.
lgH.Units='normalized';lgH.Position=[.42 .69 .27 .065];
lgV.Units='normalized';lgV.Position=[.42 .615 .27 .065];
a.Position=[.12 .13 .76 .75];b.Position=a.Position;
annotation(f,'textbox',[.12 .955 .76 .035],'String','Normalized TKE comparison','FontName','Arial','FontSize',10,'FontWeight','bold','HorizontalAlignment','center','EdgeColor','none');
set(f,'PaperUnits','centimeters','PaperSize',figureSizeCm,'PaperPosition',[0 0 figureSizeCm]);
drawnow;print(f,fullfile(outDir,'TKE_paper_style.png'),'-dpng',sprintf('-r%d',exportDPI));
print(f,fullfile(outDir,'TKE_paper_style.pdf'),'-dpdf','-painters');
savefig(f,fullfile(outDir,'TKE_paper_style.fig'));close(f);
fid=fopen(fullfile(outDir,'TKE_paper_caption.txt'),'w');
fprintf(fid,'Figure X. Normalized horizontal (red) and vertical (black) profiles of two-component resolved turbulent kinetic energy for the isothermal H/W=1 case. Lines represent uDALES and symbols represent digitized wind-tunnel reference curves. Both LES profiles are extracted from tdump.101.nc after averaging in y and over 1000-7106 s. The vertical profile is taken at the canyon centre and the horizontal profile at z/H=0.5. TKE=0.5*(u-prime variance + w-prime variance), normalized by the square of each dataset''s reference friction velocity; the candidate LES value is 0.0586 m/s. No boundary zeros are appended to TKE. The experimental vertical TKE digitization and exact measurement height remain to be confirmed before definitive validation.\n');fclose(fid);

%% VELOCITY COMPARISON USING THE SAME CANDIDATE REFERENCE VELOCITY
f=figure('Visible','off','Color','w','Units','centimeters','Position',[2 2 24 14]);
a=axes(f,'Position',[.10 .24 .35 .62]);b=axes(f,'Position',[.60 .24 .35 .62]);hold(a,'on');hold(b,'on');
plot(a,uv(zUse)/uStarLES,z(zUse)/H,'k-','LineWidth',1.6);plot(a,eU(:,2),eU(:,1),'k:','LineWidth',1.3);
plot(b,(x-xc)/H,wv/uStarLES,'r--','LineWidth',1.6);plot(b,eW(:,1),eW(:,2),'k:','LineWidth',1.3);
ylim(a,[0 1]);xlim(b,[-.5 .5]);fmt(a);fmt(b);
xlabel(a,'Mean U / u*','Interpreter','none');ylabel(a,'z/H');xlabel(b,'(x-xc)/H');ylabel(b,'Mean W / u*','Interpreter','none');
title(a,'(a) Streamwise velocity at x = 32 m');title(b,'(b) Vertical velocity at z = 8 m');
legend(a,{'uDALES','Digitized experiment'},'Location','best','Box','off');legend(b,{'uDALES','Digitized experiment'},'Location','best','Box','off');
footer(f,timeWindow,actualEnd,uStarLES,spatialMode);emit(f,outDir,'velocity_two_panel');
%% COMPOSITE VELOCITY FIGURE: THESIS STYLE (independent axes, fixed image ratio)
% Added endpoint values are prescribed no-slip boundary conditions, not sampled cells.
addNoSlipEndpoints=true; % Set false to show only extracted interior samples.
uvPlot=uv(zUse)/uStarLES;zPlot=z(zUse)/H;
wvPlot=wv/uStarLES;xPlot=(x-xc)/H;
if addNoSlipEndpoints
 uvPlot=[0;uvPlot];zPlot=[0;zPlot];
 wvPlot=[0;wvPlot;0];xPlot=[-.5;xPlot;.5];
end
f=figure('Visible','off','Color','w','Units','centimeters','Position',[2 2 figureSizeCm]);
a=axes(f,'Position',[.12 .13 .76 .75]);hold(a,'on');
p1=plot(a,xPlot,wvPlot,'r--','LineWidth',1.35);
vvalid=all(isfinite(eW),2);vv=eW(vvalid,:);jj=unique(round(linspace(1,size(vv,1),min(65,size(vv,1)))));
p2=plot(a,vv(jj,1),vv(jj,2),'r^','MarkerFaceColor','r','MarkerSize',3.5);
set(a,'XAxisLocation','top','YAxisLocation','left','XColor','r','YColor','r');
xlim(a,[-.5 .5]);xticks(a,-.5:.1:.5);ylim(a,[-11 11]);yticks(a,-10:2:10);
xlabel(a,'x / H','FontWeight','bold','FontSize',12);
ylabel(a,'V / u_*','Interpreter','tex','FontWeight','bold','FontSize',12);grid(a,'on');
b=axes(f,'Position',a.Position,'Color','none','XAxisLocation','bottom','YAxisLocation','right');hold(b,'on');
p3=plot(b,uvPlot,zPlot,'k-','LineWidth',1.35);
vvalid=all(isfinite(eU),2);vv=eU(vvalid,:);jj=unique(round(linspace(1,size(vv,1),min(55,size(vv,1)))));
p4=plot(b,vv(jj,2),vv(jj,1),'ko','MarkerFaceColor','k','MarkerSize',3);
xlim(b,[-12 12]);xticks(b,-10:5:10);ylim(b,[0 1.04]);yticks(b,0:.1:1);
xlabel(b,'U / u_*','Interpreter','tex','FontWeight','bold','FontSize',12);
ylabel(b,'z / H','FontWeight','bold','FontSize',12);grid(b,'on');
set([a b],'FontName','Arial','FontSize',9,'LineWidth',.65,'GridAlpha',.14,'TickDir','in');
% Compact legend in open upper-left area avoids covering current U curves.
legend(b,[p3 p4 p1 p2],{'U (Sim)','U (Exp)','V (Sim)','V (Exp)'},'Location','northwest','FontSize',8,'Box','on');
a.Position=[.12 .13 .76 .75];b.Position=a.Position;
set(f,'PaperUnits','centimeters','PaperSize',figureSizeCm,'PaperPosition',[0 0 figureSizeCm]);
drawnow;print(f,fullfile(outDir,'velocity_paper_style.png'),'-dpng',sprintf('-r%d',exportDPI));
print(f,fullfile(outDir,'velocity_paper_style.pdf'),'-dpdf','-painters');
savefig(f,fullfile(outDir,'velocity_paper_style.fig'));close(f);
fid=fopen(fullfile(outDir,'velocity_paper_caption.txt'),'w');
fprintf(fid,'Figure X. Normalized mean streamwise (U, black) and vertical (V, red) velocity profiles for the isothermal H/W=1 case. Lines represent uDALES and symbols represent digitized wind-tunnel data. Simulation profiles are averaged along y and over 1000-7106 s; U is extracted at the canyon centre and V at z/H=0.5. The candidate simulation friction velocity is 0.0586 m/s. Prescribed zero-velocity values are appended at the ground endpoint of U and the wall endpoints of V; interior simulation samples are unchanged. The top coordinate x/H is centred on the canyon, i.e. x=(x_model-32 m).\n');fclose(fid);

%% SAVE NUMBERS AND SETTINGS
writetable(table(z(zUse)/H,uv(zUse),uv(zUse)/uStarLES,kv(zUse),kv(zUse)/uStarLES^2,'VariableNames',{'zH','U_m_s','U_over_ustar','k2D_m2_s2','k2D_over_ustar2'}),fullfile(outDir,'simulation_vertical.csv'));
writetable(table((x-xc)/H,wv,wv/uStarLES,kh,kh/uStarLES^2,'VariableNames',{'xH','W_m_s','W_over_ustar','k2D_m2_s2','k2D_over_ustar2'}),fullfile(outDir,'simulation_horizontal.csv'));
writetable(table(eHx,eHk,'VariableNames',{'xH','k2D_over_ustar2'}),fullfile(outDir,'experiment_TKE_horizontal.csv'));
writetable(table(eVz,eVk,'VariableNames',{'zH','k2D_over_ustar2'}),fullfile(outDir,'experiment_TKE_vertical.csv'));
save(fullfile(outDir,'comparison_data.mat'),'U','W','K','x','xu','z','zw','uv','wv','kv','kh','weights','t','timeWindow','actualEnd','uStarLES','uStarExperiment','eHx','eHk','eVz','eVk','eU','eW','file','spatialMode');
fid=fopen(fullfile(outDir,'notes.txt'),'w');fprintf(fid,'Candidate u*=%.6f m/s. Requested %.0f-%.0f s, actual end %.6f s.\nY treatment: %s.\nOriginal test.m used nearest y=32; this uses configurable y-mean or y interpolation.\nTKE data: original user horizontal_TKE.csv / Vertical_TKE.csv; unchanged normalized values.\nPrior archive records flag uncertainty in absolute vertical experimental TKE digitization. No offset or clipping applied.\nHorizontal measurement height follows prior project setting z/H=0.5; exact experimental sampling should be confirmed.\nVelocity experiment uses archived digitized 0K curves.\nTKE is weighted mean of stored within-interval variances, not reconstructed long-window variance.\nNo formal validation pass/fail is inferred.\n',uStarLES,timeWindow,actualEnd,spatialMode);fclose(fid);
fprintf('Saved all figures to %s\n',outDir);
%% LOCAL HELPERS: ALL CODE REMAINS IN THIS FILE
function B=readSlab(f,n,xn,zn,ix,iz,it,yt,mode)
 v=ncinfo(f,n);dims={v.Dimensions.Name};start=ones(1,4);count=[v.Dimensions.Length];
 dx=find(strcmp(dims,xn));dy=find(strcmp(dims,'yt'));dz=find(strcmp(dims,zn));dt=find(strcmp(dims,'time'));
 start(dx)=ix(1);count(dx)=numel(ix);start(dz)=iz(1);count(dz)=numel(iz);start(dt)=it;count(dt)=1;
 if strcmp(mode,'y32'),lo=find(yt<=32,1,'last');hi=find(yt>=32,1);start(dy)=lo;count(dy)=hi-lo+1;end
 A=permute(double(ncread(f,n,start,count)),[dx dy dz dt]);
 for k=1:numel(v.Attributes)
  if any(strcmp(v.Attributes(k).Name,{'_FillValue','missing_value'})),A(A==double(v.Attributes(k).Value))=NaN;end
 end
 if strcmp(mode,'ymean'),B=squeeze(mean(A,2,'omitnan'));
 else
  if lo==hi,B=squeeze(A(:,1,:));else,q=(32-yt(lo))/(yt(hi)-yt(lo));B=squeeze((1-q)*A(:,1,:)+q*A(:,end,:));end
 end
end
function fmt(a)
 grid(a,'on');box(a,'on');set(a,'FontName','Arial','FontSize',10,'GridAlpha',.15,'LineWidth',.8);
end
function footer(f,window,te,us,mode)
 annotation(f,'textbox',[.03 .01 .94 .065],'String',sprintf('Requested %.0f-%.0f s; available to %.2f s; %s; candidate LES u* = %.4f m/s.\nExperiment: digitized reference curves; no rescaling or fitted offset.',window,te,mode,us),'EdgeColor','none','HorizontalAlignment','center','FontSize',8);
end
function emit(f,out,name)
 drawnow;exportgraphics(f,fullfile(out,[name '.png']),'Resolution',300);exportgraphics(f,fullfile(out,[name '.pdf']),'ContentType','vector');savefig(f,fullfile(out,[name '.fig']));close(f);
end

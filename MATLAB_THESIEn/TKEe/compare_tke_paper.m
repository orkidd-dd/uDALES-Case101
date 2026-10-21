clear; clc;
% Put this script, plot_tke_publication.m and all experimental CSVs together.
% udbase must already be on the MATLAB path.
here=fileparts(mfilename('fullpath'));
addpath(here);
outdir=fullfile(here,'TKE_paper_figures');
if ~exist(outdir,'dir'), mkdir(outdir); end

%% Settings
files={ ...
    'H:\udales_data\u-dales-output\101_128', ...
    'H:\udales_data\u-dales-output\101_70k', ...
    'H:\udales_data\u-dales-output\101_170k', ...
    'H:\udales_data\u-dales-output\101_240k'};
temps=[0 70 170 240];
H=16; xc=32; yc=32; zTarget=H/2;
timeWindow=[600 3500];
ustar=[0.14 0.14 0.14 0.14]; % simulation friction velocity, m/s
% CSV experiment values are ALREADY normalized: do not divide them again.
% Keep the same physical definition of friction velocity as the experiment.
recomputeSimulation=true; % false: replot saved profiles without loading tdump
cache=fullfile(outdir,'simulation_profiles.mat');

%% Experimental curves: native PDF path nodes, no smoothing
E=struct([]);
for k=1:4
    th=readtable(fullfile(here,sprintf('tke_horizontal_%dK_nodes.csv',temps(k))));
    tv=readtable(fullfile(here,sprintf('tke_vertical_%dK_nodes.csv',temps(k))));
    E(k).xH=th.x_over_H;
    E(k).horizontal=th.TKE_over_ustar2;
    E(k).zH=tv.z_over_H;
    E(k).vertical=tv.TKE_over_ustar2;
end

%% Simulation profiles, case-specific coordinates and times
settings=struct('files',{files},'temps',temps,'H',H,'xc',xc,'yc',yc, ...
    'zTarget',zTarget,'timeWindow',timeWindow,'ustar',ustar);
if recomputeSimulation
    assert(exist('udbase','file')~=0,'Add udbase to your MATLAB path first.');
    assert(all(isfinite(ustar) & ustar>0),'ustar must be positive.');
    S=struct([]);
    for k=1:4
        s=udbase(101,files{k});
        xt=s.xt(:); yt=s.yt(:); zt=s.zt(:);
        assert(all(diff(xt)>0) && all(diff(yt)>0) && all(diff(zt)>0), ...
            'Coordinates must be increasing: %d K.',temps(k));
        assert(xc>=xt(1) && xc<=xt(end) && yc>=yt(1) && yc<=yt(end) ...
            && zTarget>=zt(1) && zTarget<=zt(end),'Profile location is outside grid.');
        time=ncread(fullfile(files{k},'tdump.101.nc'),'time'); time=time(:);
        idx=find(time>=timeWindow(1) & time<=timeWindow(2));
        assert(numel(idx)>=2,'Fewer than two selected times: %d K.',temps(k));
        assert(all(diff(time)>0),'Time must increase monotonically.');
        dt=diff(time(idx));
        if max(abs(dt-median(dt)))>1e-5*max(1,abs(median(dt)))
            warning('Unequal output spacing: %d K. Using equal record weights.',temps(k));
        end
        nx=numel(xt); ny=numel(yt); nz=numel(zt);
        % Retain the supplied definition: 1/2 (<u''u''> + <w''w''>).
        % Average components separately to avoid creating another 4-D array.
        U=s.load_stat_t('upuptc');
        checkDimensions(U,nx,ny,nz,numel(time),temps(k));
        Umean=mean(U(:,:,:,idx),4,'omitnan'); clear U;
        W=s.load_stat_t('wpwptc');
        checkDimensions(W,nx,ny,nz,numel(time),temps(k));
        Wmean=mean(W(:,:,:,idx),4,'omitnan'); clear W;
        F=0.5*(Umean+Wmean); clear Umean Wmean;
        % Interpolate the y-plane rather than selecting the nearest cell.
        Y=reshape(permute(F,[2 1 3]),ny,[]);
        Fxz=reshape(interp1(yt,Y,yc,'linear'),nx,nz);
        clear F Y;
        vv=interp1(xt,Fxz,xc,'linear').';
        hh=interp1(zt,Fxz.',zTarget,'linear').';
        Lx=xt>xc-H/2 & xt<xc+H/2;
        Lz=zt>=0 & zt<=H;
        assert(any(Lx) && any(Lz),'No cells in canyon.');
        S(k).label=sprintf('%d K',temps(k));
        S(k).xH=(xt(Lx)-xc)/H;
        S(k).horizontal=hh(Lx)/ustar(k)^2;
        S(k).zH=zt(Lz)/H;
        S(k).vertical=vv(Lz)/ustar(k)^2;
        S(k).records=numel(idx);
        S(k).selectedTime=[time(idx(1)) time(idx(end))];
        S(k).ustar=ustar(k);
        fprintf('%s: %d records, %.3f to %.3f s\n',S(k).label,numel(idx), ...
            time(idx(1)),time(idx(end)));
        writetable(table(S(k).xH,S(k).horizontal,'VariableNames', ...
            {'x_over_H','TKE_over_ustar2'}), ...
            fullfile(outdir,sprintf('simulation_horizontal_%dK.csv',temps(k))));
        writetable(table(S(k).zH,S(k).vertical,'VariableNames', ...
            {'z_over_H','TKE_over_ustar2'}), ...
            fullfile(outdir,sprintf('simulation_vertical_%dK.csv',temps(k))));
        clear s Fxz;
    end
    save(cache,'S','settings');
else
    C=load(cache,'S','settings');
    assert(isequaln(C.settings,settings), ...
        'Settings differ from cache. Set recomputeSimulation=true.');
    S=C.S;
end

%% Export publication figures and descriptive pointwise differences
figs=plot_tke_publication(S,E,outdir,true);
% Interpolate experiment onto simulation locations; do not extrapolate.
% RMSE/MAE below are unweighted across simulation nodes, not significance tests.
rows=cell(8,5); r=0;
for k=1:4
    for direction=1:2
        if direction==1
            q=S(k).zH; v=S(k).vertical; qe=E(k).zH; ve=E(k).vertical;
            name='Vertical';
        else
            q=S(k).xH; v=S(k).horizontal; qe=E(k).xH; ve=E(k).horizontal;
            name='Horizontal';
        end
        ref=interp1(qe,ve,q,'linear',NaN);
        ok=isfinite(ref) & isfinite(v); err=v(ok)-ref(ok);
        r=r+1;
        if isempty(err), rmse=NaN; mae=NaN;
        else, rmse=sqrt(mean(err.^2)); mae=mean(abs(err)); end
        rows(r,:)={temps(k),name,sum(ok),rmse,mae};
    end
end
metrics=cell2table(rows,'VariableNames',{'DeltaT_K','Direction','N','RMSE','MAE'});
writetable(metrics,fullfile(outdir,'comparison_metrics.csv')); disp(metrics);
fprintf('Saved figures and profiles to:\n%s\n',outdir);

function checkDimensions(A,nx,ny,nz,nt,temp)
assert(size(A,1)==nx && size(A,2)==ny && size(A,3)==nz && size(A,4)==nt, ...
    'Expected [nx ny nz nt] for %d K; got [%s].',temp,num2str(size(A)));
end

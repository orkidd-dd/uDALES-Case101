% Digitized experimental mean velocities, Figure 11(a).
% V in the paper is the VERTICAL velocity (usually w in uDALES).
clear; clc;
here=fileparts(mfilename('fullpath'));
colors=[.1 .1 .1;0 .35 .65;.80 .35 0;.45 .20 .60];
styles={':','-','--','-.'}; temps=[0 70 170 240];
f=figure('Color','w','Units','centimeters','Position',[2 2 18 8.5]);
tl=tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');
for panel=1:2
    a=nexttile(tl); hold(a,'on');
    for k=1:4
        if panel==1
            T=readtable(fullfile(here,sprintf('U_vertical_%dK_nodes.csv',temps(k))));
            valid=T.z_over_H>=0&T.z_over_H<=1;
            plot(a,T.U_over_ustar(valid),T.z_over_H(valid), ...
                'Color',colors(k,:),'LineStyle',styles{k},'LineWidth',1.25, ...
                'DisplayName',sprintf('%d K',temps(k)));
        else
            T=readtable(fullfile(here,sprintf('V_horizontal_%dK_nodes.csv',temps(k))));
            plot(a,T.x_over_H,T.V_over_ustar, ...
                'Color',colors(k,:),'LineStyle',styles{k},'LineWidth',1.25, ...
                'DisplayName',sprintf('%d K',temps(k)));
        end
    end
    if panel==1
        xlabel(a,'$U/u_*$','Interpreter','latex'); ylabel(a,'$z/H$','Interpreter','latex');
        xlim(a,[-5.2 5.2]); ylim(a,[0 1]);
        title(a,'(a) Horizontal velocity at x/H = 0','FontWeight','normal');
    else
        xlabel(a,'$x/H$','Interpreter','latex'); ylabel(a,'$V/u_*$','Interpreter','latex');
        xlim(a,[-.5 .5]); ylim(a,[-5 5]);
        title(a,'(b) Vertical velocity at z/H = 0.5','FontWeight','normal');
    end
    set(a,'FontName','Times New Roman','FontSize',10,'LineWidth',.75, ...
        'TickDir','in','Box','on');
    legend(a,'Location','best','Box','off','FontSize',8);
end
exportgraphics(f,fullfile(here,'mean_velocity_experiment.pdf'),'ContentType','vector');
exportgraphics(f,fullfile(here,'mean_velocity_experiment.png'),'Resolution',600);
savefig(f,fullfile(here,'mean_velocity_experiment.fig'));

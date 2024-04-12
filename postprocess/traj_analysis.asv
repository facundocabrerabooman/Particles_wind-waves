clear;clc;close all

% Set path were functions will be read from
addpath(genpath('/Users/fcb/Documents/GitHub/Particles_wind-waves/'));

fname = 'basecase';

folderin = '/Users/fcb/Library/CloudStorage/GoogleDrive-facundo@pdx.edu/My Drive/Particles_Waves&Wind/PTV/exports/basecase/';
folderout = folderin;
cd(folderin)

Fs=500; % Frame rate

%% Concatenate data

if pi==pi
    trajs_conc = [];
    
    load("trajs_basecase.mat")
    trajs_conc = [trajs_conc tracklong];
    

Ine=find(arrayfun(@(X)(~isempty(X.Vx)),trajs_conc)==1);
trajs_conc = trajs_conc(Ine);

end
clear tracklong traj_ddt
save('traj_conc_tracers_ddt','trajs_conc','Ine','-v7.3')

%% Get Mean Velocity 

for i=1:numel(trajs_conc)
    trajs_conc(i).X = trajs_conc(i).x;
    trajs_conc(i).Y = trajs_conc(i).y;
    trajs_conc(i).Z = trajs_conc(i).z;
end

[U, mBdt, bins] = track2meanDxDt3DProfile(trajs_conc,'Xf',2:2:50,[20 20 20],1,1,'x','cart');
[V, ~, ~] = track2meanDxDt3DProfile(trajs_conc,'Yf',2:2:50,[20 20 20],1,1,'y','cart');
[W, ~, ~] = track2meanDxDt3DProfile(trajs_conc,'Zf',2:2:50,[20 20 20],1,1,'z','cart');

[X,Y,Z]=meshgrid(bins{1},bins{2},bins{3});

U=U.*Fs;
V=V.*Fs;
W=W.*Fs;

%save('output_Vel_meanfields','U','V','W','mBdt','bins','X','Y','Z')


trajs_conc_minus_mean_field = find_closest_bin(trajs_conc, X, Y, Z, U, V, W);

%save('trajs_conc_minus_mean_field','trajs_conc_minus_mean_field','-v7.3')
%% nice colors
mycolormap = mycolor('#063970','#e28743');%('#063970','#eeeee4','#e28743')
color3 = [mycolormap(1,:);mycolormap((size(mycolormap,1)+1)/2,:);mycolormap(end,:)];
color1 = '#476d76';
%% Plot 3D Trajectories

figure(10); clf; hold on; grid on; box on

%streamw_vel_mean = mean(vertcat(trajs_conc_new_axis.Vx));
streamw_vel_mean = mean(vertcat(trajs_conc.Vx));

for i=1:numel(trajs_conc)

    %scatter3(trajs_conc_new_axis(i).Xf,trajs_conc_new_axis(i).Yf,trajs_conc_new_axis(i).Zf, 10, trajs_conc_new_axis(i).Vx - streamw_vel_mean, 'filled');
    scatter3(trajs_conc(i).Xf,trajs_conc(i).Yf,trajs_conc(i).Zf, 10, trajs_conc(i).Vy, 'filled');
    
end
xlabel('X normalized (streamwise)')
ylabel('Y normalized (vertical - antigravity)')
zlabel('Z normalized (spanwise)')
 
colorbar

savefig_FC('traj3D',8,6,'pdf')
savefig_FC('traj3D',8,6,'fig')
%% Plot Position Histograms

figure(2); clf; hold on; grid on; box on

    histogram(vertcat(trajs_conc.Xf),'FaceColor',color3(1,:))
    histogram(vertcat(trajs_conc.Yf),'FaceColor',color3(2,:))
    histogram(vertcat(trajs_conc.Zf),'FaceColor',color3(3,:))

    legend({'X ','Y','Z'})

savefig_FC('histogram_pos',8,6,'pdf')
savefig_FC('histogram_pos',8,6,'fig')
%% Plot Velocity Histograms

figure(3); clf; hold on; grid on; box on

    % histogram(vertcat(trajs_conc.Vx),'FaceColor','r')
    % histogram(vertcat(trajs_conc.Vy),'FaceColor','g')
    % histogram(vertcat(trajs_conc.Vz),'FaceColor','b')

    histogram(vertcat(trajs_conc.Vx),'FaceColor',color3(1,:))
    histogram(vertcat(trajs_conc.Vy),'FaceColor',color3(2,:))
    histogram(vertcat(trajs_conc.Vz),'FaceColor',color3(3,:))

    legend({'Vx','Vy','Vz'})

savefig_FC('histogram_vel',8,6,'pdf')
savefig_FC('histogram_vel',8,6,'fig')


%% 2 times - 1 particle statistics (Lagrangian statistics)
%% Mean Square Separation
MSD(1) = structFunc_struct(trajs_conc,'Xf',2);
MSD(2) = structFunc_struct(trajs_conc,'Yf',2);
MSD(3) = structFunc_struct(trajs_conc,'Zf',2);

save('output_post_processing.mat','MSD','-append')
%%
figure;
loglog(MSD(1).tau/Fs,MSD(1).mean,'d-',MarkerSize=8,Color=color3(1,:),LineWidth=2);hold on
loglog(MSD(2).tau/Fs,MSD(2).mean,'d-',MarkerSize=8,Color=color3(2,:),LineWidth=2);
loglog(MSD(3).tau/Fs,MSD(3).mean,'d-',MarkerSize=8,Color=color3(3,:),LineWidth=2);

xMSD = linspace(1,100,1000)/Fs;
loglog(xMSD,2e5*xMSD.^2,'--',Color=color1,LineWidth=2)

set(gca,FontSize=15)
legend('MSDx','MSDy','MSDz',Location='best',FontSize=12)
title('$MSD$','interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$MSD(m^2)$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$\tau(s)$','interpreter','latex',FontWeight='bold',FontSize=24)
text(2e-3,3,'$\tau^2$','interpreter','latex',FontWeight='bold',FontSize=18)
grid on
axis padded


folderout = 'MSS/';
mkdir(folderout)
savefig_custom([folderout 'MSS'],8,6,'pdf')
savefig_custom([folderout 'MSS'],8,6,'fig')
%% Longitudinal S2

S2L(1)= structFunc_struct(trajs_conc,'Vx',2);
S2L(2)= structFunc_struct(trajs_conc,'Vy',2);
S2L(3)= structFunc_struct(trajs_conc,'Vz',2);

save('output_post_processing.mat','S2L','-append')
%%
% figure;loglog(S2Lx.tau,S2Lx.mean./S2Lx.tau/Fs/2)
figure;
loglog(S2L(1).tau/Fs,S2L(1).mean,'d-',MarkerSize=8,Color=color3(1,:),LineWidth=2);hold on
loglog(S2L(2).tau/Fs,S2L(2).mean,'d-',MarkerSize=8,Color=color3(2,:),LineWidth=2);
loglog(S2L(3).tau/Fs,S2L(3).mean,'d-',MarkerSize=8,Color=color3(3,:),LineWidth=2);

xS2L = linspace(1,9,100)/Fs;
loglog(xS2L,6e8*xS2L.^2,'--',Color=color1,LineWidth=2)
xS2L = linspace(10,200,100)/Fs;
loglog(xS2L,2e6*xS2L.^1,'--',Color=color1,LineWidth=2)
% xS2L = linspace(100,300,100);
% loglog(xS2L,8e4*xS2L.^0,'--',Color=color1,LineWidth=2)

set(gca,FontSize=15)
legend('$S_2^L(x)$','$S_2^L(y)$','$S_2^L(z)$','interpreter','latex',Location='best',FontSize=12)
title('$S_2^L$','interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$S_2^L$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$\tau/s$','interpreter','latex',FontWeight='bold',FontSize=24)
text(8e-4,1e3,'$\tau^2$','interpreter','latex',FontWeight='bold',FontSize=18)
text(1e-2,4e4,'$\tau$','interpreter','latex',FontWeight='bold',FontSize=18)
grid on
axis padded

folderout = 'S2L/';
mkdir(folderout)
savefig_custom([folderout 'S2L'],8,6,'pdf')
savefig_custom([folderout 'S2L'],8,6,'fig')

%% Velocity and Acceleration Correlations
if 1==pi
disp('corr')

Ine=find(arrayfun(@(X)(~isempty(X.Vx)),trajs_conc_minus_mean_field)==1);


Ruu(1) = xcorr_struct(trajs_conc_minus_mean_field(Ine),'Vx',1);
Ruu(2) = xcorr_struct(trajs_conc_minus_mean_field(Ine),'Vy',1);
Ruu(3) = xcorr_struct(trajs_conc_minus_mean_field(Ine),'Vz',1);

Raa(1) = xcorr_struct(trajs_conc_minus_mean_field(Ine),'Ax',1);
Raa(2) = xcorr_struct(trajs_conc_minus_mean_field(Ine),'Ay',1);
Raa(3) = xcorr_struct(trajs_conc_minus_mean_field(Ine),'Az',1);

%%% another option: 
% n=1;
% [Rvx,D2x,Nptsvx,Ntrackvx]=lagstats_tracks(trajs_conc,'Vx',n,'tabsframes');
% [Rvy,D2y,Nptsvy,Ntrackvy]=lagstats_tracks(trajs_conc,'Vy',n,'tabsframes');
% [Rvz,D2z,Nptsvz,Ntrackvz]=lagstats_tracks(trajs_conc,'Vz',n,'tabsframes');
% 
% [Rax,~,Nptsax,Ntrackax]=lagstats_tracks(traj_conc_0_f,'Ax',n,'tabsframes');
% [Ray,~,Nptsay,Ntrackay]=lagstats_tracks(traj_conc_0_f,'Ay',n,'tabsframes');
% [Raz,~,Nptsaz,Ntrackaz]=lagstats_tracks(traj_conc_0_f,'Az',n,'tabsframes');

%% Correlation fit
Fs=2990;
Ruufit(1) = correlationFit(Ruu(1),Fs,1,100,'V');
Ruufit(2) = correlationFit(Ruu(2),Fs,1,100,'V');
Ruufit(3) = correlationFit(Ruu(3),Fs,1,100,'V');

Raafit(1) = correlationFit(Raa(1),Fs,1,100,'A');
Raafit(2) = correlationFit(Raa(2),Fs,1,100,'A');
Raafit(3) = correlationFit(Raa(3),Fs,1,100,'A');
%%
figure;
% main plot: zoom in

plot(Ruu(1).tau/Fs,Ruu(1).mean/Ruu(1).mean(1),'d',MarkerSize=8,Color=color3(1,:),LineWidth=2);hold on
plot(Ruu(2).tau/Fs,Ruu(2).mean/Ruu(2).mean(1),'d',MarkerSize=8,Color=color3(2,:),LineWidth=2);
plot(Ruu(3).tau/Fs,Ruu(3).mean/Ruu(3).mean(1),'d',MarkerSize=8,Color=color3(3,:),LineWidth=2);

plot(Raa(1).tau/Fs,Raa(1).mean/Raa(1).mean(1),'^',MarkerSize=8,Color=color3(1,:),LineWidth=2);hold on
plot(Raa(2).tau/Fs,Raa(2).mean/Raa(2).mean(1),'^',MarkerSize=8,Color=color3(2,:),LineWidth=2);
plot(Raa(3).tau/Fs,Raa(3).mean/Raa(3).mean(1),'^',MarkerSize=8,Color=color3(3,:),LineWidth=2);

plot(Ruufit(1).x,Ruufit(1).yfit,'-',Color=color3(1,:),LineWidth=2);hold on
plot(Ruufit(2).x,Ruufit(2).yfit,'-',Color=color3(2,:),LineWidth=2)
plot(Ruufit(3).x,Ruufit(3).yfit,'-',Color=color3(3,:),LineWidth=2)

plot(Raafit(1).x,Raafit(1).yfit,'-',Color=color3(1,:),LineWidth=2)
plot(Raafit(2).x,Raafit(2).yfit,'-',Color=color3(2,:),LineWidth=2)
plot(Raafit(3).x,Raafit(3).yfit,'-',Color=color3(3,:),LineWidth=2)


set(gca,FontSize=15)
legend('$R_{uu}(x)1$','$R_{uu}(y)1$','$R_{uu}(z)1$','$R_{aa}(x)1$','$R_{aa}(y)1$','$R_{aa}(z)1$','interpreter','latex',Location='best',FontSize=12)
title('$R_{uu}, R_{aa}$','interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$R_{uu}, R_{aa}$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$\tau$/s','interpreter','latex',FontWeight='bold',FontSize=24)
grid on
axis tight
xlim([0 5e-2]) 
ylim([-0.1 1.1])


% add inset: zoom out 
axes('Position',[0.4 0.5 0.3 0.2]);
plot(Ruu(1).tau/Fs,Ruu(1).mean/Ruu(1).mean(1),'d',MarkerSize=3,Color=color3(1,:),LineWidth=1);hold on
plot(Ruu(2).tau/Fs,Ruu(2).mean/Ruu(2).mean(1),'d',MarkerSize=3,Color=color3(2,:),LineWidth=1);
plot(Ruu(3).tau/Fs,Ruu(3).mean/Ruu(3).mean(1),'d',MarkerSize=3,Color=color3(3,:),LineWidth=1);

plot(Raa(1).tau/Fs,Raa(1).mean/Raa(1).mean(1),'^',MarkerSize=3,Color=color3(1,:),LineWidth=1);
plot(Raa(2).tau/Fs,Raa(2).mean/Raa(2).mean(1),'^',MarkerSize=3,Color=color3(2,:),LineWidth=1);
plot(Raa(3).tau/Fs,Raa(3).mean/Raa(3).mean(1),'^',MarkerSize=3,Color=color3(3,:),LineWidth=1);

plot(Ruufit(1).x,Ruufit(1).yfit,'-',Color=color3(1,:),LineWidth=1);hold on
plot(Ruufit(2).x,Ruufit(2).yfit,'-',Color=color3(2,:),LineWidth=1)
plot(Ruufit(3).x,Ruufit(3).yfit,'-',Color=color3(3,:),LineWidth=1)

plot(Raafit(1).x,Raafit(1).yfit,'-',Color=color3(1,:),LineWidth=1)
plot(Raafit(2).x,Raafit(2).yfit,'-',Color=color3(2,:),LineWidth=1)
plot(Raafit(3).x,Raafit(3).yfit,'-',Color=color3(3,:),LineWidth=1)
set(gca,FontSize=12)
grid on
axis tight

ylim([-0.1 1.1])
xlim([0 0.2])


folderout = 'corr/';
mkdir(folderout)
savefig_custom([folderout 'corr'],8,6,'pdf')
savefig_custom([folderout 'corr'],8,6,'fig')

save('output_post_processing.mat','Ruu','Raa','Ruufit','Raafit','-append')
end
%% Eulerian 2-point statistics
%clearvars -except trajs_conc Ine Ine Fs color1 color3 folderin folderout mycolormap

%for j=1:numel(trajs_conc); trajs_conc(j).Tf = trajs_conc(j).t_sec_abs; end % rename Tf field

tic  
[eulerStats, pair] = twoPointsEulerianStats_Mica_Speedup(trajs_conc,[0.5 40],40,'off');
toc
save('output_post_processing.mat','eulerStats','pair','-append')

%clearvars -except eulerStats pair trajs_conc Ine Fs folderin folderout color3 color1
%% Plot
figure;
loglog(eulerStats.r,eulerStats.S2x,'d-',MarkerSize=8,Color=color3(1,:),LineWidth=2);hold on
loglog(eulerStats.r,eulerStats.S2y,'d-',MarkerSize=8,Color=color3(2,:),LineWidth=2);
loglog(eulerStats.r,eulerStats.S2z,'d-',MarkerSize=8,Color=color3(3,:),LineWidth=2);

rS2E = linspace(0.5,40,100);
loglog(rS2E,6e3*rS2E.^1,'--',Color=color1,LineWidth=2)

set(gca,FontSize=15)
legend('$S_2^E(x)$','$S_2^E(y)$','$S_2^E(z)$','interpreter','latex',Location='best',FontSize=12)
title('$S_2^E$','interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$S_2^E$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$r/mm$','interpreter','latex',FontWeight='bold',FontSize=24)
text(4,4e4,'$r^1$','interpreter','latex',FontWeight='bold',FontSize=18)
grid on
axis padded

folderout = 'S2euler/';
mkdir(folderout)
savefig_custom([folderout 'S2e'],8,6,'pdf')
savefig_custom([folderout 'S2e'],8,6,'fig')

figure;hold on;
semilogy(eulerStats.r,(eulerStats.Splong{2}./2.1).^(3/2)./eulerStats.r'./1e6,'o-')
%semilogy(eulerStats.r,eulerStats.Splong{2},'o-')
grid;
xlabel('r (mm)','Interpreter','latex');
ylabel('$(S_2^E^\parallel / C_2 r^{2/3})^{3/2}$','Interpreter','latex');
set(gca,'FontSize',24);
set(gca,'Xscale','log','Yscale','log');
title('Compensated Eulerian S2ps');
fname = 'compensated_S2_epsilon';
 

folderout = 'S2euler/';
mkdir(folderout)
savefig_custom([folderout 'S2el_compensated_epsilon'],8,6,'pdf')
savefig_custom([folderout 'S2el_compensated_epsilon'],8,6,'fig')


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Added on 10/04/2023
%% don't forget to change the figure saving parts ..
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot VAt (to check stationary)

figure
t=tiledlayout(4,1,'TileSpacing','tight');
nexttile;
plot(eulerStats.Vmoy,'d-',MarkerSize=3,Color=color1(1,:),LineWidth=1);hold on
plot(eulerStats.VmoyX,'-',MarkerSize=3,Color=color3(1,:),LineWidth=1);
plot(eulerStats.VmoyY,'-',MarkerSize=3,Color=color3(2,:),LineWidth=1);
plot(eulerStats.VmoyZ,'-',MarkerSize=3,Color=color3(3,:),LineWidth=1);
set(gca,FontSize=15)
legend('$\langle \sqrt{x^2+y^2+z^2} \rangle$','$\langle x \rangle$','$\langle y \rangle$','$\langle z \rangle$','interpreter','latex',Location='best',FontSize=10);
title('$|V|, \sigma_V, |A|,\sigma_A$','interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$|V|$','interpreter','latex',FontWeight='bold',FontSize=18)
axis tight
xticklabels([])
grid on

nexttile;
plot(eulerStats.Vstd,'d-',MarkerSize=3,Color=color1(1,:),LineWidth=1);hold on
plot(eulerStats.VstdX,'-',MarkerSize=3,Color=color3(1,:),LineWidth=1);
plot(eulerStats.VstdY,'-',MarkerSize=3,Color=color3(2,:),LineWidth=1);
plot(eulerStats.VstdZ,'-',MarkerSize=3,Color=color3(3,:),LineWidth=1);
set(gca,FontSize=15)
ylabel('$\sigma_V$','interpreter','latex',FontWeight='bold',FontSize=18)
axis tight
xticklabels([])
grid on

nexttile;
plot(eulerStats.Amoy,'d-',MarkerSize=3,Color=color1(1,:),LineWidth=1);hold on
plot(eulerStats.AmoyX,'-',MarkerSize=3,Color=color3(1,:),LineWidth=1);
plot(eulerStats.AmoyY,'-',MarkerSize=3,Color=color3(2,:),LineWidth=1);
plot(eulerStats.AmoyZ,'-',MarkerSize=3,Color=color3(3,:),LineWidth=1);
set(gca,FontSize=15)
ylabel('$|A|$','interpreter','latex',FontWeight='bold',FontSize=18)
axis tight
xticklabels([])
grid on

nexttile;
plot(eulerStats.Astd,'d-',MarkerSize=3,Color=color1(1,:),LineWidth=1);hold on
plot(eulerStats.AstdX,'-',MarkerSize=3,Color=color3(1,:),LineWidth=1);
plot(eulerStats.AstdY,'-',MarkerSize=3,Color=color3(2,:),LineWidth=1);
plot(eulerStats.AstdZ,'-',MarkerSize=3,Color=color3(3,:),LineWidth=1);
set(gca,FontSize=15)
ylabel('$\sigma_A$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$t/s$','interpreter','latex',FontWeight='bold',FontSize=18)
axis tight
xticks(0:Fs/10:size(eulerStats.Astd,2))
xticklabels(num2cell([0:Fs/10:size(eulerStats.Astd,2)]/Fs))
grid on


linkaxes(t.Children,'x')

% figname = ['./' fout '/VAt'];
% savefig(figname)
% saveas(gcf,figname,'png')
% saveas(gcf,figname,'pdf')

%% plot Spn_abs
color5 = [mycolormap(1,:);mycolormap(round((size(mycolormap,1)+1)/4),:);mycolormap(round(2*(size(mycolormap,1)+1)/4),:);mycolormap(round(3*(size(mycolormap,1)+1)/4),:);mycolormap(end,:)];

figure;
loglog(eulerStats.r,eulerStats.SplongAbs{1,1},'d-',MarkerSize=8,Color=color5(1,:),LineWidth=2);hold on
loglog(eulerStats.r,eulerStats.SplongAbs{1,2},'d-',MarkerSize=8,Color=color5(2,:),LineWidth=2);
loglog(eulerStats.r,eulerStats.SplongAbs{1,3},'d-',MarkerSize=8,Color=color5(3,:),LineWidth=2);
loglog(eulerStats.r,eulerStats.SplongAbs{1,4},'d-',MarkerSize=8,Color=color5(4,:),LineWidth=2);
loglog(eulerStats.r,eulerStats.SplongAbs{1,5},'d-',MarkerSize=8,Color=color5(5,:),LineWidth=2);

rSplong = linspace(0.4,100,100);
loglog(rSplong,6e1*rSplong.^(1/3),'--',Color=color1,LineWidth=2)
loglog(rSplong,5e3*rSplong.^(2/3),'--',Color=color1,LineWidth=2)
loglog(rSplong,9e5*rSplong.^(3/3),'--',Color=color1,LineWidth=2)
loglog(rSplong,1e8*rSplong.^(4/3),'--',Color=color1,LineWidth=2)
loglog(rSplong,3e10*rSplong.^(5/3),'--',Color=color1,LineWidth=2)

set(gca,FontSize=15)
legend('$S_1^{\parallel}$','$S_2^{\parallel}$','$S_3^{\parallel}$','$S_4^{\parallel}$','$S_5^{\parallel}$','interpreter','latex',Location='best',FontSize=12)
title('$S_n^{\parallel}$','interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$S_n^{\parallel}$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$r/mm$','interpreter','latex',FontWeight='bold',FontSize=24)
text(8,5e12,'$r^{5/3}$','interpreter','latex',FontWeight='bold',FontSize=18)
text(8,2e10,'$r^{4/3}$','interpreter','latex',FontWeight='bold',FontSize=18)
text(8,5e7,'$r^{3/3}$','interpreter','latex',FontWeight='bold',FontSize=18)
text(8,3e5,'$r^{2/3}$','interpreter','latex',FontWeight='bold',FontSize=18)
text(8,5e2,'$r^{1/3}$','interpreter','latex',FontWeight='bold',FontSize=18)
grid on
axis padded

% figname = ['./' fout '/SplongAbs'];
% savefig(figname)
% saveas(gcf,figname,'png')
% saveas(gcf,figname,'pdf')

%% plot Spn
figure;
loglog(eulerStats.r,eulerStats.Splong{1,1},'d-',MarkerSize=8,Color=color5(1,:),LineWidth=1);hold on
loglog(eulerStats.r,eulerStats.Splong{1,2},'d-',MarkerSize=8,Color=color5(2,:),LineWidth=1);
loglog(eulerStats.r,eulerStats.Splong{1,3},'d-',MarkerSize=8,Color=color5(3,:),LineWidth=1);
loglog(eulerStats.r,eulerStats.Splong{1,4},'d-',MarkerSize=8,Color=color5(4,:),LineWidth=1);
loglog(eulerStats.r,eulerStats.Splong{1,5},'d-',MarkerSize=8,Color=color5(5,:),LineWidth=1);

rSplong = linspace(0.4,100,100);
loglog(rSplong,6e1*rSplong.^(1/3),'--',Color=color1,LineWidth=1)
loglog(rSplong,5e3*rSplong.^(2/3),'--',Color=color1,LineWidth=1)
loglog(rSplong,9e5*rSplong.^(3/3),'--',Color=color1,LineWidth=1)
loglog(rSplong,1e8*rSplong.^(4/3),'--',Color=color1,LineWidth=1)
loglog(rSplong,3e10*rSplong.^(5/3),'--',Color=color1,LineWidth=1)

set(gca,FontSize=15)
legend('$S_1^{\parallel}$','$S_2^{\parallel}$','$S_3^{\parallel}$','$S_4^{\parallel}$','$S_5^{\parallel}$','interpreter','latex',Location='best',FontSize=12)
title('$S_n^{\parallel}$','interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$S_n^{\parallel}$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$r/mm$','interpreter','latex',FontWeight='bold',FontSize=24)
text(8,5e12,'$r^{5/3}$','interpreter','latex',FontWeight='bold',FontSize=18)
text(8,2e10,'$r^{4/3}$','interpreter','latex',FontWeight='bold',FontSize=18)
text(8,5e7,'$r^{3/3}$','interpreter','latex',FontWeight='bold',FontSize=18)
text(8,3e5,'$r^{2/3}$','interpreter','latex',FontWeight='bold',FontSize=18)
text(8,5e2,'$r^{1/3}$','interpreter','latex',FontWeight='bold',FontSize=18)
grid on
axis padded

% figname = ['./' fout '/Splong'];
% savefig(figname)
% saveas(gcf,figname,'png')
% saveas(gcf,figname,'pdf')

%% plot Sau
figure
semilogx(eulerStats.r,eulerStats.Sau,'d-',MarkerSize=8,Color=color3(1,:),LineWidth=1);hold on
semilogx(eulerStats.r,eulerStats.Saulong,'d-',MarkerSize=8,Color=color3(3,:),LineWidth=1);hold on

set(gca,FontSize=15)
legend('$S_{au}$','$S_{au}^{\parallel}$','interpreter','latex',Location='best',FontSize=12)
title('$S_{au}, S_{au}^{\parallel}$','interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$S_{au}, S_{au}^{\parallel}$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$r/mm$','interpreter','latex',FontWeight='bold',FontSize=24)
grid on
axis padded

% figname = ['./' fout '/Sau'];
% savefig(figname)
% saveas(gcf,figname,'png')
% saveas(gcf,figname,'pdf')

%% plot S2E
figure;
loglog(eulerStats.r,eulerStats.S2x,'d-',MarkerSize=8,Color=color3(1,:),LineWidth=1);hold on
loglog(eulerStats.r,eulerStats.S2y,'d-',MarkerSize=8,Color=color3(2,:),LineWidth=1);
loglog(eulerStats.r,eulerStats.S2z,'d-',MarkerSize=8,Color=color3(3,:),LineWidth=1);

loglog(eulerStats.r,7e3*eulerStats.r.^(2/3),'--',Color=color1,LineWidth=1)

set(gca,FontSize=15)
legend('$S_2^E(x)$','$S_2^E(y)$','$S_2^E(z)$','interpreter','latex',Location='best',FontSize=12)
title('$S_2^E$','interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$S_2^E$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$r/mm$','interpreter','latex',FontWeight='bold',FontSize=24)
text(3,3e4,'$r^{2/3}$','interpreter','latex',FontWeight='bold',FontSize=18)
grid on
axis padded

% figname = ['./' fout '/S2E'];
% savefig(figname)
% saveas(gcf,figname,'png')
% saveas(gcf,figname,'pdf')
%% plot epsilon
Ckolomogrov = 2.1;
figure;
loglog(eulerStats.r,(eulerStats.Splong{1,2}./Ckolomogrov).^(1.5)./eulerStats.r','d-',MarkerSize=8,Color=color3(1,:),LineWidth=2);

hold on
%figure
loglog(eulerStats.r,abs(eulerStats.Splong{1,3})./(4/5*eulerStats.r)','d-',MarkerSize=8,Color=color3(2,:),LineWidth=2);
loglog(eulerStats.r,abs(eulerStats.Sau)./2,'d-',MarkerSize=8,Color=color3(3,:),LineWidth=2);

set(gca,FontSize=15)
legend('$(S_2^{\parallel}/C_k)^{3/2}\cdot r^{-1}$','$|S_3^{\parallel}|\cdot (4/5r)^{-1}$','$|S_{au}|/2$','interpreter','latex',Location='best',FontSize=12)
title('$\epsilon$','interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$\epsilon$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$r/mm$','interpreter','latex',FontWeight='bold',FontSize=24)
grid on;axis padded


% figname = ['./' fout '/Epsilon'];
% savefig(figname)
% saveas(gcf,figname,'png')
% saveas(gcf,figname,'pdf')



%% plot Ruu(r) --- the data could be wrong, have to check the code 'twoPointEulerianStats_Mica'
figure;
plot(eulerStats.Ruur,eulerStats.Ruu,'-',MarkerSize=8,Color=color3(1,:),LineWidth=1);hold on
set(gca,FontSize=15)
% legend('$Ruu(r)$','interpreter','latex',Location='best',FontSize=12)
title('$Ruu(r)$','interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$Ruu(r)$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$Ruu_r$','interpreter','latex',FontWeight='bold',FontSize=24)
grid on
axis padded


figname = ['./' fout '/Ruu'];
savefig(figname)
saveas(gcf,figname,'png')
saveas(gcf,figname,'pdf')
%% plot PSD(k) --- the data could be wrong, have to check the code 'twoPointEulerianStats_Mica'
figure
loglog(eulerStats.PSDk,eulerStats.PSD,'-',MarkerSize=8,Color=color3(1,:),LineWidth=1);hold on

set(gca,FontSize=15)
% legend('$S_2^E(x)$','$S_2^E(y)$','$S_2^E(z)$','interpreter','latex',Location='best',FontSize=12)
title('$PSD$','interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$PSD$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$k$','interpreter','latex',FontWeight='bold',FontSize=24)
grid on
axis padded

figname = ['./' fout '/PSD'];
savefig(figname)
saveas(gcf,figname,'png')
saveas(gcf,figname,'pdf')
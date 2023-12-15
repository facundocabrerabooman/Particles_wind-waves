clear;clc;close all

% Set path were functions will be read from
addpath(genpath('/Users/FC/Documents/GitHub/Particle-laden-turbulence'));

% Set as current directory the folder with the data
folderin = '/Users/FC/Library/CloudStorage/GoogleDrive-facundo@pdx.edu/My Drive/Drop Tower Multiphase Flow Project/Data to Shake-the-Box/PostProcessing/all_concatenated/';
folderout = '/Users/FC/Library/CloudStorage/GoogleDrive-facundo@pdx.edu/My Drive/Drop Tower Multiphase Flow Project/Data to Shake-the-Box/PostProcessing/ddt_all_concatenated/';
cd(folderout)

Fs=2990; % Frame rate

%% Set cool colors for plots
mycolormap = mycolor('#063970','#e28743');%('#063970','#eeeee4','#e28743')
color3 = [mycolormap(1,:);mycolormap((size(mycolormap,1)+1)/2,:);mycolormap(end,:)];
color1 = '#476d76';
%% Load data

%%%%%%%%% If want to concatenate data use this
if 1==pi
dconc = [];
load('/Users/FC/Library/CloudStorage/GoogleDrive-facundo@pdx.edu/My Drive/Drop Tower Multiphase Flow Project/Data to Shake-the-Box/PostProcessing/july7a/july7a_tracers.mat')
dconc = vertcat(dconc,d); 
load('/Users/FC/Library/CloudStorage/GoogleDrive-facundo@pdx.edu/My Drive/Drop Tower Multiphase Flow Project/Data to Shake-the-Box/PostProcessing/july7b/july7b_tracers.mat')
dconc = vertcat(dconc,d); 
load('/Users/FC/Library/CloudStorage/GoogleDrive-facundo@pdx.edu/My Drive/Drop Tower Multiphase Flow Project/Data to Shake-the-Box/PostProcessing/july7c/july7c_tracers.mat')
dconc = vertcat(dconc,d); 

d=dconc; 
save('july7_conc','d')
clear d
end
%%%%%%%%%

fname = 'july7_conc';
%load(fname);

%% Track Particles (i.e. go from particle positions to trajectories)

maxdist = 1;  
lmin=10;
flag_pred=0;
npriormax=4;
porder=3;
flag_conf=1;
numFrames = 5e6;

[traj,tracks]=track3d_fc_stb(folderout,fname,maxdist,lmin,flag_pred,npriormax,porder,flag_conf, numFrames, Fs);


save('output_post_processing.mat','traj','tracks')
clearvars -except traj Fs folderin folderout color3 color1
%% Only keep long tracks -- redundant if using track3d_fc_stb.m
L = arrayfun(@(X)(numel(X.x)),traj);
Ilong = find(L>=10);
%% Find proper filter width
[s(1), m(1), w]=findFilterWidth_PTV(traj(Ilong),'x');
[s(2), m(2), w]=findFilterWidth_PTV(traj(Ilong),'y');
[s(3), m(3), w]=findFilterWidth_PTV(traj(Ilong),'z');
%%
figure;
yyaxis left
loglog(w,s(1).vx,'d-',MarkerSize=8,Color=color3(1,:),LineWidth=2);hold on;
loglog(w,s(2).vx,'d-',MarkerSize=8,Color=color3(2,:),LineWidth=2);
loglog(w,s(3).vx,'d-',MarkerSize=8,Color=color3(3,:),LineWidth=2);
hold off

yyaxis right
loglog(w,s(1).ax,'^-',MarkerSize=8,Color=color3(1,:),LineWidth=2);hold on;
loglog(w,s(2).ax,'^-',MarkerSize=8,Color=color3(2,:),LineWidth=2);
loglog(w,s(3).ax,'^-',MarkerSize=8,Color=color3(3,:),LineWidth=2);

plot([4 4],ylim,'--',Color=color1,LineWidth=2)

set(gca,FontSize=15)
yyaxis left
legend('$V_x$','$V_y$','$V_z$','$A_x$','$A_y$','$A_z$','interpreter','latex',Location='southwest',FontSize=12);
title('$std.(w)$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$fliter\ width\ w$','interpreter','latex',FontWeight='bold',FontSize=18)

yyaxis left
ylabel('$\sigma_{v}$','interpreter','latex',FontWeight='bold',FontSize=24)
yyaxis right
ylabel('$\sigma_{a}$','interpreter','latex',FontWeight='bold',FontSize=24)

grid on
axis padded

folderout = 'filter/';
mkdir(folderout)
savefig_custom([folderout 'filter_check'],8,6,'pdf')
savefig_custom([folderout 'filter_check'],8,6,'fig')
save('output_post_processing.mat','s','m','w','-append')

clearvars -except traj Fs folderin folderout Ilong color3 color1
%%  Estimate filtered tracks, velocities and accelerations with optimal filter width
wopt = 4;
lopt = 12;

w_acc = 10;
l_acc = 30;

%%%
%tracklong=calcVelLEM(traj(Ilong),wopt,lopt,Fs); % does not give you time

%[~, tracklong]=compute_vel_acc_traj(traj(Ilong),Fs,wopt,lopt);

tracklong=calcVelLEM(traj,wopt,lopt,Fs, w_acc, l_acc);

Ine=find(arrayfun(@(X)(~isempty(X.Vx)),tracklong)==1);
Ine_acc=find(arrayfun(@(X)(~isempty(X.Ax)),tracklong)==1);


save([folderout filesep 'output_post_processing.mat'],'Ine','Ine_acc','tracklong','-append')
clearvars -except tracklong Ine Fs folderin folderout  color3 color1

%% Split Drop Tower data 

[traj_dec, traj_ddt, traj_fullg] = split_ddt_data(folderin, folderout);

save([folderout filesep 'output_post_processing.mat'],'traj_dec','traj_ddt','traj_fullg','-append')



% filtered data goes up to 4.05s versus 7s before...



%load([folderin filesep 'output_post_processing'],'traj_dec','traj_ddt','traj_fullg')


%% If data is already split
load('output_post_processing_split.mat')

tracklong = traj_ddt;
Ine=find(arrayfun(@(X)(~isempty(X.Vx)),tracklong)==1);
Ine_acc=find(arrayfun(@(X)(~isempty(X.Ax)),tracklong)==1);


%% 1 time - 1 particle statistics
%% Calculate & plot velocity and acceleration pdfs
pdfV(1) = mkpdf5(tracklong(Ine),'Vx',256,10);
pdfV(2) = mkpdf5(tracklong(Ine),'Vy',256,10);
pdfV(3) = mkpdf5(tracklong(Ine),'Vz',256,10);

pdfA(1) = mkpdf5(tracklong(Ine_acc),'Ax',256,20);
pdfA(2) = mkpdf5(tracklong(Ine_acc),'Ay',256,20);
pdfA(3) = mkpdf5(tracklong(Ine_acc),'Az',256,20);

save('output_post_processing_ddt.mat','pdfV','pdfA','tracklong','Ine','Ine_acc')
%% Plot Normalized PDFs
figure;
semilogy(pdfV(1).xpdfn,pdfV(1).pdfn,'d-',MarkerSize=5,Color=color3(1,:),LineWidth=2);hold on;
semilogy(pdfV(2).xpdfn,pdfV(2).pdfn,'d-',MarkerSize=5,Color=color3(2,:),LineWidth=2);
semilogy(pdfV(3).xpdfn,pdfV(3).pdfn,'d-',MarkerSize=5,Color=color3(3,:),LineWidth=2);

semilogy(pdfA(1).xpdfn,pdfA(1).pdfn,'^-',MarkerSize=5,Color=color3(1,:),LineWidth=2);
semilogy(pdfA(2).xpdfn,pdfA(2).pdfn,'^-',MarkerSize=5,Color=color3(2,:),LineWidth=2);
semilogy(pdfA(3).xpdfn,pdfA(3).pdfn,'^-',MarkerSize=5,Color=color3(3,:),LineWidth=2);

xpdf=linspace(-5,5,1024);
plot(xpdf,normpdf(xpdf,0,1),Color=color1,LineWidth=2);

set(gca,FontSize=15)
legend('$V_x$','$V_y$','$V_z$','$A_x$','$A_y$','$A_z$','interpreter','latex',Location='best',FontSize=12);
title('$PDF$','interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$PDF(V,A)$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$V, A$','interpreter','latex',FontWeight='bold',FontSize=18)
grid on
axis padded

% text(5,1,['MeanAX = ' num2str(pdfA(1).mean)])
% text(5,0.6,['MeanAY = ' num2str(pdfA(2).mean)])
% text(5,0.3,['MeanAZ = ' num2str(pdfA(3).mean)])
% 
% text(5,0.1,['MeanVX = ' num2str(pdfV(1).mean)])
% text(5,0.05,['MeanVY = ' num2str(pdfV(2).mean)])
% text(5,0.03,['MeanVZ = ' num2str(pdfV(3).mean)])

% add subfigure
axes('Position',[0.22 0.62 0.22 0.22]);
semilogy(pdfV(1).xpdfn,pdfV(1).pdfn,'d-',MarkerSize=2,Color=color3(1,:),LineWidth=2);hold on;
semilogy(pdfV(2).xpdfn,pdfV(2).pdfn,'d-',MarkerSize=2,Color=color3(2,:),LineWidth=2);
semilogy(pdfV(3).xpdfn,pdfV(3).pdfn,'d-',MarkerSize=2,Color=color3(3,:),LineWidth=2);

semilogy(pdfA(1).xpdfn,pdfA(1).pdfn,'^-',MarkerSize=2,Color=color3(1,:),LineWidth=2);
semilogy(pdfA(2).xpdfn,pdfA(2).pdfn,'^-',MarkerSize=2,Color=color3(2,:),LineWidth=2);
semilogy(pdfA(3).xpdfn,pdfA(3).pdfn,'^-',MarkerSize=2,Color=color3(3,:),LineWidth=2);

xpdf=linspace(-5,5,1024);
plot(xpdf,normpdf(xpdf,0,1),'k',LineWidth=2);
grid on
set(gca,FontSize=12)
xlim([-5 5])

folderout = 'pdfs/';
mkdir(folderout)
savefig_custom([folderout 'PDFs'],8,6,'pdf')
savefig_custom([folderout 'PDFs'],8,6,'fig')

%% Table with moments of distribution
maketable(pdfA,pdfV,folderout)
%% 
% figure;
% 
% semilogy(pdfV(1).xpdf,pdfV(1).pdf,'d-',MarkerSize=5,Color=color3(1,:),LineWidth=2);hold on;
% semilogy(pdfV(2).xpdf,pdfV(2).pdf,'d-',MarkerSize=5,Color=color3(2,:),LineWidth=2);
% semilogy(pdfV(3).xpdf,pdfV(3).pdf,'d-',MarkerSize=5,Color=color3(3,:),LineWidth=2);
% 
% xpdfn.V1 = linspace(-1.2,1.0,1024)*1e3;
% xpdfn.V2 = linspace(-1.0,1.0,1024)*1e3;
% xpdfn.V3 = linspace(-1.1,0.8,1024)*1e3;
% semilogy(xpdfn.V1,normpdf(xpdfn.V1,pdfV(1).mean,pdfV(1).std),'--',Color=color3(1,:),LineWidth=2);
% semilogy(xpdfn.V2,normpdf(xpdfn.V2,pdfV(2).mean,pdfV(2).std),'--',Color=color3(2,:),LineWidth=2);
% semilogy(xpdfn.V3,normpdf(xpdfn.V3,pdfV(3).mean,pdfV(3).std),'--',Color=color3(3,:),LineWidth=2);
% 
% % lavision output
% % m/s to mm/s
% % load('LavisionOutput\tracers.mat')
% [pdfVL.x,xpdfVL.x] = hist(d(:,5),256); 
% [pdfVL.y,xpdfVL.y] = hist(d(:,6),256); 
% [pdfVL.z,xpdfVL.z] = hist(d(:,7),256); 
% semilogy(xpdfVL.x.*1e3,pdfVL.x/sum(pdfVL.x),'-',Color=color3(1,:),LineWidth=2); hold on
% semilogy(xpdfVL.y.*1e3,pdfVL.y/sum(pdfVL.y),'-',Color=color3(2,:),LineWidth=2); 
% semilogy(xpdfVL.z.*1e3,pdfVL.z/sum(pdfVL.z),'-',Color=color3(3,:),LineWidth=2); 
% 
% set(gca,FontSize=15)
% legend('$V_x$','$V_y$','$V_z$','$Gfitx$','$Gfity$','$Gfitz$','$V_xL$','$V_yL$','$V_zL$','interpreter','latex',Location='best',FontSize=12);
% title('$PDFn$','interpreter','latex',FontWeight='bold',FontSize=18)
% ylabel('$PDF(V)$','interpreter','latex',FontWeight='bold',FontSize=18)
% xlabel('$V(mm/s)$','interpreter','latex',FontWeight='bold',FontSize=18)
% grid on
% axis padded

%%
% figure
% semilogy(pdfA(1).xpdf,pdfA(1).pdf,'^-',MarkerSize=5,Color=color3(1,:),LineWidth=2);hold on
% semilogy(pdfA(2).xpdf,pdfA(2).pdf,'^-',MarkerSize=5,Color=color3(2,:),LineWidth=2);
% semilogy(pdfA(3).xpdf,pdfA(3).pdf,'^-',MarkerSize=5,Color=color3(3,:),LineWidth=2);
% 
% xpdfn.A1 = linspace(-0.8e5,0.8e5,1024);
% xpdfn.A2 = linspace(-0.8e5,0.8e5,1024);
% xpdfn.A3 = linspace(-0.8e5,0.8e5,1024);
% semilogy(xpdfn.A1,normpdf(xpdfn.A1,pdfA(1).mean,pdfA(1).std),'--',MarkerSize=5,Color=color3(1,:),LineWidth=2);
% semilogy(xpdfn.A2,normpdf(xpdfn.A2,pdfA(2).mean,pdfA(2).std),'--',MarkerSize=5,Color=color3(2,:),LineWidth=2);
% semilogy(xpdfn.A3,normpdf(xpdfn.A3,pdfA(3).mean,pdfA(3).std),'--',MarkerSize=5,Color=color3(3,:),LineWidth=2);
% 
% 
% % lavision output
% % m/s to mm/s
% [pdfAL.x,xpdfAL.x] = hist(d(:,8),256); 
% [pdfAL.y,xpdfAL.y] = hist(d(:,9),256); 
% [pdfAL.z,xpdfAL.z] = hist(d(:,10),256); 
% semilogy(xpdfAL.x*1e3,pdfAL.x/sum(pdfAL.x),'-',Color=color3(1,:),LineWidth=2); hold on
% semilogy(xpdfAL.y*1e3,pdfAL.y/sum(pdfAL.y),'-',Color=color3(2,:),LineWidth=2); 
% semilogy(xpdfAL.z*1e3,pdfAL.z/sum(pdfAL.z),'-',Color=color3(3,:),LineWidth=2); 
% 
% set(gca,FontSize=15)
% legend('$A_x$','$A_y$','$A_z$','$Gfitx$','$Gfity$','$Gfitz$','$A_xL$','$A_yL$','$A_zL$','interpreter','latex',Location='best',FontSize=12);
% title('$PDFn$','interpreter','latex',FontWeight='bold',FontSize=18)
% ylabel('$PDF(A)n$','interpreter','latex',FontWeight='bold',FontSize=18)
% xlabel('$A(mm/s^2)$','interpreter','latex',FontWeight='bold',FontSize=18)
% grid on
% axis padded
% % xlim([-5 5])
% 
% % savefig('PDFs')
% % saveas(gcf,'PDFs','png')

%% 2 times - 1 particle statistics (Lagrangian statistics)

%% Mean Square Separation
MSD(1) = structFunc_struct(tracklong(Ine),'Xf',2);
MSD(2) = structFunc_struct(tracklong(Ine),'Yf',2);
MSD(3) = structFunc_struct(tracklong(Ine),'Zf',2);

save('output_post_processing_ddt.mat','MSD','-append')
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

S2L(1)= structFunc_struct(tracklong(Ine),'Vx',2);
S2L(2)= structFunc_struct(tracklong(Ine),'Vy',2);
S2L(3)= structFunc_struct(tracklong(Ine),'Vz',2);

save('output_post_processing_ddt.mat','S2L','-append')
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
Ruu(1) = xcorr_struct(tracklong(Ine),'Vx',1);
Ruu(2) = xcorr_struct(tracklong(Ine),'Vy',1);
Ruu(3) = xcorr_struct(tracklong(Ine),'Vz',1);

Raa(1) = xcorr_struct(tracklong(Ine_acc),'Ax',1);
Raa(2) = xcorr_struct(tracklong(Ine_acc),'Ay',1);
Raa(3) = xcorr_struct(tracklong(Ine_acc),'Az',1);

%%% another option: 
% n=1;
% [Rvx,D2x,Nptsvx,Ntrackvx]=lagstats_tracks(tracklong(Ine),'Vx',n,'tabsframes');
% [Rvy,D2y,Nptsvy,Ntrackvy]=lagstats_tracks(tracklong(Ine),'Vy',n,'tabsframes');
% [Rvz,D2z,Nptsvz,Ntrackvz]=lagstats_tracks(tracklong(Ine),'Vz',n,'tabsframes');
% 
% [Rax,~,Nptsax,Ntrackax]=lagstats_tracks(traj_conc_0_f,'Ax',n,'tabsframes');
% [Ray,~,Nptsay,Ntrackay]=lagstats_tracks(traj_conc_0_f,'Ay',n,'tabsframes');
% [Raz,~,Nptsaz,Ntrackaz]=lagstats_tracks(traj_conc_0_f,'Az',n,'tabsframes');

%% Correlation fit
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
xlim([0 12e-3])
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
xlim([0 1e-1])


folderout = 'corr/';
mkdir(folderout)
savefig_custom([folderout 'corr'],8,6,'pdf')
savefig_custom([folderout 'corr'],8,6,'fig')

save('output_post_processing_ddt.mat','Ruu','Raa','Ruufit','Raafit','-append')
%% Eulerian 2-point statistics
clearvars -except tracklong Ine Fs Ine_acc

% if using data that has acceleration: (remember acc. is filtered more than vel)
for j=1:numel(tracklong); tracklong(j).Tf = tracklong(j).Tf_acc; end % rename Tf field

tic  
[eulerStats, pair] = twoPointsEulerianStats_Mica_Speedup(tracklong(Ine_acc),[0.5 40],40,'off');
toc
save('output_post_processing_ddt.mat','eulerStats','pair','-append')

clearvars -except eulerStats pair tracklong Ine Ine_acc Fs folderin folderout color3 color1
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
%%
%% Added on 10/04/2023
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

folderout = 'Vat/';
mkdir(folderout)
savefig_custom([folderout 'Vat'],8,6,'pdf')
savefig_custom([folderout 'Vat'],8,6,'fig')

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


folderout = 'euler/';
mkdir(folderout)
savefig_custom([folderout 'SplongAbs'],8,6,'pdf')
savefig_custom([folderout 'SplongAbs'],8,6,'fig')

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

folderout = 'euler/';
mkdir(folderout)
savefig_custom([folderout 'Splong'],8,6,'pdf')
savefig_custom([folderout 'Splong'],8,6,'fig')

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

folderout = 'euler/';
mkdir(folderout)
savefig_custom([folderout 'Sau'],8,6,'pdf')
savefig_custom([folderout 'Sau'],8,6,'fig')

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

folderout = 'euler/';
mkdir(folderout)
savefig_custom([folderout 'S2E'],8,6,'pdf')
savefig_custom([folderout 'S2E'],8,6,'fig')
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


folderout = 'euler/';
mkdir(folderout)
savefig_custom([folderout 'epsilon'],8,6,'pdf')
savefig_custom([folderout 'epsilon'],8,6,'fig')



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


folderout = 'euler/';
mkdir(folderout)
savefig_custom([folderout 'Ruu'],8,6,'pdf')
savefig_custom([folderout 'Ruu'],8,6,'fig')
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

folderout = 'euler/';
mkdir(folderout)
savefig_custom([folderout 'PSD'],8,6,'pdf')
savefig_custom([folderout 'PSD'],8,6,'fig')




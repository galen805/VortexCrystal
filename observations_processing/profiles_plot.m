% profiles and parameters for N and S pole
clear;clc;close all

% set up the figures

% main figure with vorticity and velocity
f1 = figure('pos',[200,200,1200,1400]);
t1=tiledlayout(f1,4,2,'TileSpacing','compact','Padding','compact');

% supplement figure for shape param
fb = figure('pos',[200,200,1200,600]);
tb=tiledlayout(fb,1,2,'TileSpacing','compact','Padding','compact');


%% south pole
load filtered_data_south.mat
xmin = min(Xg,[],'all');
xmax = max(Xg,[],'all');
ymin = min(Yg,[],'all');
ymax = max(Yg,[],'all');
[Ny, Nx] = size(Xg);
dx = Xg(1,2) - Xg(1,1);   % assume uniform grid
dy = Yg(2,1) - Yg(1,1);

speed = sqrt(u_specfilt.^2 + v_specfilt.^2);
f  = -3.5e-4;

%% find vortex centers
% center is defined by grad(psi)=0 which is exactly min(speed)

center_guess=[-600,1550;...
    -7400,2600;...
    -5300,-4700;...
    4000,-5000;...
    7000,2400]*1e3;
Rmax = 500e3;   % search radius
numK=size(center_guess,1);

for k=1:numK
    x0 = center_guess(k,1);
    y0 = center_guess(k,2);

    % distance from guess point to every grid point
    dist2 = (Xg - x0).^2 + (Yg - y0).^2;
    mask  = dist2 <= Rmax^2;   % within 100 km

    speed_cand = speed(mask);
    [~, ii]  = min(speed_cand);          % minimum speed within 100 km

    ind_all  = find(mask);
    ind_min  = ind_all(ii);

    [i_cen, j_cen] = ind2sub(size(speed), ind_min);
    
    centers(k,:) = [i_cen, j_cen];

    i0 = centers(k,1);
    j0 = centers(k,2);

    centers_x0(k) = Xg(i0,j0);
    centers_y0(k) = Yg(i0,j0);

    % polar position
    centers_r0(k)     = hypot(centers_x0(k), centers_y0(k));
    centers_theta0(k) = atan2d(centers_y0(k), centers_x0(k));

end

fprintf('\nGeometry parameters:\n')
fprintf('r0_pc = %.2f\n',centers_r0(1))
fprintf('r0_cpc = %.2f\n',mean(centers_r0(2:end)))


figure(f1)
nexttile(t1,1)
hold on
levels=[-1:0.1:1];
contourf(Xg, Yg,vort_specfilt_calc/f,levels,'edgecolor','none');
contour(Xg,Yg,psi_specfilt_calc,'linewidth',1,'linecolor','k')
plot(centers_x0,centers_y0,'o','linewidth',2,'markersize',5,'color','m','MarkerFaceColor','m')
text(centers_x0(1)+300e3,centers_y0(1)+100e3,'SPC','fontsize',12,'color','m')
text(centers_x0(2)+300e3,centers_y0(2)+100e3,'CPC1','fontsize',12,'color','m')
text(centers_x0(3)+300e3,centers_y0(3)+100e3,'CPC2','fontsize',12,'color','m')
text(centers_x0(4)+300e3,centers_y0(4)+100e3,'CPC3','fontsize',12,'color','m')
text(centers_x0(5)+300e3,centers_y0(5)+100e3,'CPC4','fontsize',12,'color','m')

discreteColorbar(gca,levels,cmocean('balance'),'eastoutside','\zeta/f',4);
axis equal
% xlim([xmin,xmax])
% ylim([ymin,ymax])
xlabel('x (km)');
ylabel('y (km)');
title('(a) SP Rossby number');
set(gca,'fontsize',14);

nexttile(t1,3)
hold on
levels = [-100:10:100];
contourf(Xg, Yg, speed,levels,'edgecolor','none');
contour(Xg,Yg,psi_specfilt_calc,'linewidth',1,'linecolor','k')
plot(centers_x0,centers_y0,'o','linewidth',2,'markersize',5,'color','m','MarkerFaceColor','m')
text(centers_x0(1)+300e3,centers_y0(1)+100e3,'SPC','fontsize',12,'color','m')
text(centers_x0(2)+300e3,centers_y0(2)+100e3,'CPC1','fontsize',12,'color','m')
text(centers_x0(3)+300e3,centers_y0(3)+100e3,'CPC2','fontsize',12,'color','m')
text(centers_x0(4)+300e3,centers_y0(4)+100e3,'CPC3','fontsize',12,'color','m')
text(centers_x0(5)+300e3,centers_y0(5)+100e3,'CPC4','fontsize',12,'color','m')

discreteColorbar(gca,levels,cmocean('balance'),'eastoutside','m/s',2);
axis equal
% xlim([xmin,xmax])
% ylim([ymin,ymax])
xlabel('x (km)');
ylabel('y (km)');
title('(c) SP speed');
set(gca,'fontsize',14);

%% derive profiles about each center
% parameters
dr   = 100e3;        % ring width [m]
Rmax = 5000e3;      % max radius [m]
rb = 0:dr:Rmax;
rc = rb(1:end-1) + dr/2;

circulation_invert

%% Derive observed Vm and Rm for each cyclone

fitpars = nan(numK,2);   % [Vm Rm]

for k = 1:numK

    rdat = rc_keep{k};
    vdat = ut_inv{k};

    valid_prof = isfinite(rdat) & isfinite(vdat) & rdat >= 0;
    rdat = rdat(valid_prof);
    vdat = vdat(valid_prof);

    % Observed maximum tangential velocity and its radius
    [Vm_obs,ind] = max(abs(vdat));
    Rm_obs = rdat(ind);

    fitpars(k,:) = [Vm_obs,Rm_obs];
end

% Shared observed properties for CPCs; independent properties for PC
Vm_cpc = mean(fitpars(2:numK,1),'omitnan');
Rm_cpc = mean(fitpars(2:numK,2),'omitnan');

Vm_pc = fitpars(1,1);
Rm_pc = fitpars(1,2);

fprintf('\nObserved velocity parameters:\n')
fprintf('Vm_cpc = %.4f\n',Vm_cpc)
fprintf('Rm_cpc = %.2f\n',Rm_cpc)
fprintf('Vm_pc  = %.4f\n',Vm_pc)
fprintf('Rm_pc  = %.2f\n',Rm_pc)


%% Find best-fit steepness parameters using signed velocity RMSE

valid = isfinite(Xg) & isfinite(Yg) & ...
        isfinite(u_specfilt) & isfinite(v_specfilt);

xdat = Xg(valid);
ydat = Yg(valid);

u_obs = -1*u_specfilt(valid);
v_obs = -1*v_specfilt(valid);

% Define grid search
b_pc_vec  = 0.25:0.05:3;
b_cpc_vec = 0.25:0.05:3;

bfit

fprintf('\nSigned-velocity grid-search best fit:\n')
fprintf('b_pc  = %.2f\n',b_pc_best)
fprintf('b = %.2f\n',b_cpc_best)
fprintf('Minimum signed velocity RMSE = %.5f m/s\n',min_rmse)


%% Plot RMSE solution space
figure(fb)
nexttile(tb,1)
hold on
levels=[10:1:20];
cmap=cmocean('thermal',length(levels));
cmap=cmap(1:(length(levels)-1),:);
contourf(b_pc_vec,b_cpc_vec,rmse_map,levels,'linecolor','none')
plot(b_pc_best,b_cpc_best,'kp','markersize',20,'markerfacecolor','w','linewidth',2)
xlabel('b_{pc}')
ylabel('b_{cpc}')
axis equal
title('(a) SP solution space for b')
grid on
box on
set(gca,'fontsize',14)

cb1 = discreteColorbar(gca, levels, cmap, 'eastoutside', 'RMSE (m/s)', 4);


%% plot profiles and best fit

Vcyclone = @(Vm,Rm,r,b) Vm./Rm .* r .* ...
    exp((1./b) .* (1 - (r./Rm).^b));

Zcyclone = @(Vm,Rm,r,b) 2*Vm./Rm .* ...
    (1 - 0.5*(r./Rm).^b) .* ...
    exp((1./b) .* (1 - (r./Rm).^b));

% Best-fit analytical profiles
rfit = linspace(0,Rmax,600);
ut_pc_fit   = Vcyclone(Vm_pc, Rm_pc, rfit, b_pc_best);
ut_cpc_fit  = Vcyclone(Vm_cpc,Rm_cpc,rfit, b_cpc_best);
vort_pc_fit  = Zcyclone(Vm_pc, Rm_pc, rfit, b_pc_best);
vort_cpc_fit = Zcyclone(Vm_cpc,Rm_cpc,rfit, b_cpc_best);

figure(f1)
nexttile(t1,2)
hold on
for k=1:numK
    plot(rc_keep{k},vort_mean{k},'linewidth',2)
end
plot(rfit,vort_pc_fit,'k-','linewidth',2)
plot(rfit,vort_cpc_fit,'k--','linewidth',2)
legend('SPC','CPC1','CPC2','CPC3','CPC4','SPC fit','CPC fit','fontsize',12)
xlabel('Radius (m)')
ylabel('1/s')
title('(b) SP cyclonic vorticity');
grid on
box on
set(gca,'fontsize',14);

nexttile(t1,4)
hold on
for k=1:numK
    plot(rc_keep{k},ut_inv{k},'linewidth',2)
end
plot(rfit,ut_pc_fit,'k-','linewidth',2)
plot(rfit,ut_cpc_fit,'k--','linewidth',2)
legend('SPC','CPC1','CPC2','CPC3','CPC4','SPC fit','CPC fit','fontsize',12)
xlabel('Radius (m)')
ylabel('m/s')
title('(d) SP cyclonic velocity');
grid on
box on
set(gca,'fontsize',14);

%% north pole
clearvars -except f1 fb t1 tb
% load the data
load filtered_data_north.mat
xmin = min(Xg,[],'all');
xmax = max(Xg,[],'all');
ymin = min(Yg,[],'all');
ymax = max(Yg,[],'all');
[Ny, Nx] = size(Xg);
dx = Xg(1,2) - Xg(1,1);   % assume uniform grid
dy = Yg(2,1) - Yg(1,1);

speed = sqrt(u_specfilt.^2 + v_specfilt.^2);
f  = 3.5e-4;

%% find vortex centers
% center is defined by grad(psi)=0 which is exactly min(speed)

center_guess=[240,-24;...
    -7072,-1576;...
    7872,1624]*1e3;
Rmax = 500e3;   % search radius
numK=size(center_guess,1);

for k=1:numK
    x0 = center_guess(k,1);
    y0 = center_guess(k,2);

    % distance from guess point to every grid point
    dist2 = (Xg - x0).^2 + (Yg - y0).^2;
    mask  = dist2 <= Rmax^2;   % within 100 km

    speed_cand = speed(mask);
    [~, ii]  = min(speed_cand);          % minimum speed within 100 km

    ind_all  = find(mask);
    ind_min  = ind_all(ii);

    [i_cen, j_cen] = ind2sub(size(speed), ind_min);
    
    centers(k,:) = [i_cen, j_cen];

    i0 = centers(k,1);
    j0 = centers(k,2);

    centers_x0(k) = Xg(i0,j0);
    centers_y0(k) = Yg(i0,j0);

    % polar position
    centers_r0(k)     = hypot(centers_x0(k), centers_y0(k));
    centers_theta0(k) = atan2(centers_y0(k), centers_x0(k));

end

fprintf('\nGeometry parameters:\n')
fprintf('r0_pc = %.2f km\n',centers_r0(1)/1e3)
fprintf('r0_cpc = %.2f km\n',mean(centers_r0(2:end))/1e3)

figure(f1)
nexttile(t1,5,[1,1])
hold on
levels=[-1:0.1:1];
contourf(Xg, Yg,vort_specfilt_calc/f,levels,'edgecolor','none');
contour(Xg,Yg,psi_specfilt_calc,'linewidth',1,'linecolor','k')
plot(centers_x0,centers_y0,'o','linewidth',2,'markersize',5,'color','m','MarkerFaceColor','m')
text(centers_x0(1)+300e3,centers_y0(1)+100e3,'NPC','fontsize',12,'color','m')
text(centers_x0(2)+300e3,centers_y0(2)+100e3,'CPC1','fontsize',12,'color','m')
text(centers_x0(3)+300e3,centers_y0(3)+100e3,'CPC2','fontsize',12,'color','m')

discreteColorbar(gca,levels,cmocean('balance'),'eastoutside','\zeta/f',4);
axis equal
% xlim([xmin,xmax])
% ylim([ymin,ymax])
xlabel('x (km)');
ylabel('y (km)');
title('(e) NP Rossby number');
set(gca,'fontsize',14);


nexttile(t1,7,[1,1])
hold on
levels = [-100:10:100];
contourf(Xg, Yg, speed,levels,'edgecolor','none');
contour(Xg,Yg,psi_specfilt_calc,'linewidth',1,'linecolor','k')
plot(centers_x0,centers_y0,'o','linewidth',2,'markersize',5,'color','m','MarkerFaceColor','m')
text(centers_x0(1)+300e3,centers_y0(1)+100e3,'NPC','fontsize',12,'color','m')
text(centers_x0(2)+300e3,centers_y0(2)+100e3,'CPC1','fontsize',12,'color','m')
text(centers_x0(3)+300e3,centers_y0(3)+100e3,'CPC2','fontsize',12,'color','m')

discreteColorbar(gca,levels,cmocean('balance'),'eastoutside','m/s',4);
axis equal
% xlim([xmin,xmax])
% ylim([ymin,ymax])
xlabel('x (km)');
ylabel('y (km)');
title('(g) NP speed');
set(gca,'fontsize',14);


%% invert circulation to derive profiles about each center
% parameters
dr   = 100e3;        % ring width [m]
Rmax = 5000e3;      % max radius [m]
rb = 0:dr:Rmax;
rc = rb(1:end-1) + dr/2;

circulation_invert


%% Derive observed Vm and Rm for each cyclone

fitpars = nan(numK,2);   % [Vm Rm]

for k = 1:numK

    rdat = rc_keep{k};
    vdat = ut_inv{k};

    valid_prof = isfinite(rdat) & isfinite(vdat) & rdat >= 0;
    rdat = rdat(valid_prof);
    vdat = vdat(valid_prof);

    % Observed maximum tangential velocity and its radius
    [Vm_obs,ind] = max(abs(vdat));
    Rm_obs = rdat(ind);

    fitpars(k,:) = [Vm_obs,Rm_obs];
end

% Shared observed properties for CPCs; independent properties for PC
Vm_cpc = mean(fitpars(2:numK,1),'omitnan');
Rm_cpc = mean(fitpars(2:numK,2),'omitnan');

Vm_pc = fitpars(1,1);
Rm_pc = fitpars(1,2);

fprintf('\nObserved velocity parameters:\n')
fprintf('Vm_cpc = %.4f m/s\n',Vm_cpc)
fprintf('Rm_cpc = %.2f km\n',Rm_cpc/1e3)
fprintf('Vm_pc  = %.4f m/s\n',Vm_pc)
fprintf('Rm_pc  = %.2f km\n',Rm_pc/1e3)


%% Find best-fit steepness parameters using signed velocity RMSE

valid = isfinite(Xg) & isfinite(Yg) & ...
        isfinite(u_specfilt) & isfinite(v_specfilt);

xdat = Xg(valid);
ydat = Yg(valid);

u_obs = u_specfilt(valid);
v_obs = v_specfilt(valid);

% Define grid search
b_pc_vec  = 0.25:0.05:3;
b_cpc_vec = 0.25:0.05:3;

bfit

fprintf('\nSigned-velocity grid-search best fit:\n')
fprintf('b_pc  = %.2f\n',b_pc_best)
fprintf('b = %.2f\n',b_cpc_best)
fprintf('Minimum signed velocity RMSE = %.5f m/s\n',min_rmse)


%% Plot RMSE solution space

figure(fb)
nexttile(tb,2)
hold on
levels=[10:1:20];
cmap=cmocean('thermal',length(levels));
cmap=cmap(1:(length(levels)-1),:);
contourf(b_pc_vec,b_cpc_vec,rmse_map,levels,'linecolor','none')
plot(b_pc_best,b_cpc_best,'kp','markersize',20,'markerfacecolor','w','linewidth',2)
xlabel('b_{pc}')
ylabel('b_{cpc}')
axis equal
title('(b) NP solution space for b')
grid on
box on
set(gca,'fontsize',14)

cb = discreteColorbar(gca, levels, cmap, 'eastoutside', 'RMSE (m/s)', 2);



%% plot profiles and best fit

Vcyclone = @(Vm,Rm,r,b) Vm./Rm .* r .* ...
    exp((1./b) .* (1 - (r./Rm).^b));

Zcyclone = @(Vm,Rm,r,b) 2*Vm./Rm .* ...
    (1 - 0.5*(r./Rm).^b) .* ...
    exp((1./b) .* (1 - (r./Rm).^b));

% Best-fit analytical profiles
rfit = linspace(0,Rmax,600);
ut_pc_fit   = Vcyclone(Vm_pc, Rm_pc, rfit, b_pc_best);
ut_cpc_fit  = Vcyclone(Vm_cpc,Rm_cpc,rfit, b_cpc_best);
vort_pc_fit  = Zcyclone(Vm_pc, Rm_pc, rfit, b_pc_best);
vort_cpc_fit = Zcyclone(Vm_cpc,Rm_cpc,rfit, b_cpc_best);

figure(f1)
nexttile(t1,6)
hold on
for k=1:numK
    plot(rc_keep{k},vort_mean{k},'linewidth',2)
end
plot(rfit,vort_pc_fit,'k-','linewidth',2)
plot(rfit,vort_cpc_fit,'k--','linewidth',2)
legend('NPC','CPC1','CPC2','NPC fit','CPC fit','fontsize',12)
xlabel('Radius (m)')
ylabel('1/s')
title('(f) NP cyclonic vorticity');
grid on
box on
set(gca,'fontsize',14);

nexttile(t1,8)
hold on
for k=1:numK
    plot(rc_keep{k},ut_inv{k},'linewidth',2)
end
plot(rfit,ut_pc_fit,'k-','linewidth',2)
plot(rfit,ut_cpc_fit,'k--','linewidth',2)
legend('NPC','CPC1','CPC2','NPC fit','CPC fit','fontsize',12)
xlabel('Radius (m)')
ylabel('m/s')
title('(h) NP cyclonic velocity');
grid on
box on
set(gca,'fontsize',14);



%% size and export
% set(f1,'pos',[0,0,900,1000])
set(f1,'pos',[0,0,1200,1350])
set(f1,'visible','off')
exportgraphics(f1,'plots/alldata.png','resolution',500)

set(fb,'pos',[0,0,1200,600])
exportgraphics(fb,'plots/bfit.png','resolution',500)


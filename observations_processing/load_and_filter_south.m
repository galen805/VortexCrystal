% load and filter the data
clear;clc;%close all

% SOUTH POLE
f  = -3.5e-4;

load PJ4_south_r4_1961_1961_lon_lat_vx_vy_vor_v0.mat

u = double(u);
v = double(v);
vorticity = double(vorticity); % we will overwrite this
lat = double(lat);
lon0 = double(lon);

%% plot raw data
figure;
set(gcf,'pos',[200,200,1200,600])

tiledlayout(1,2,'TileSpacing','compact','Padding','compact')

nexttile
hold on
dat = vorticity/f;
m_proj('stereographic','lat',-90,'long',0,'radius',12);
m_pcolor(lon,lat,dat,'linestyle','none');
set(gca,'fontsize',14)
cmocean('balance')
h=colorbar('southoutside');
clim([-1,1]);
h.Label.String ='Rossby number \zeta/f';
m_grid('FontSize',14,'FontWeight','bold','tickdir','out',...
    'xtick',12,'ytick',6,...
    'XaxisLocation','top',...
    'YaxisLocation','right',...
    'linest','-','box','on');
set(gca,'fontsize',14)
drawnow;

nexttile
hold on
dat = sqrt(u.^2+v.^2);
m_proj('stereographic','lat',-90,'long',0,'radius',12);
m_pcolor(lon,lat,dat,'linestyle','none');
cmocean('balance')
h=colorbar('southoutside');
clim([-100,100]);
h.Label.String ='Speed m/s';
m_grid('FontSize',14,'FontWeight','bold','tickdir','out',...
    'xtick',12,'ytick',6,...
    'XaxisLocation','top',...
    'YaxisLocation','right',...
    'linest','-','box','on');
set(gca,'fontsize',14)
drawnow;

exportgraphics(gcf,'plots/rawdata_S.png','resolution',500)

%% grid data
% Rotate the longitude coordinate
alpha = 48;
lon = lon0 - alpha;

R = 7.1492e7;

phi = deg2rad(lat);
lam = deg2rad(lon);

% South-pole distance from the pole: positive outward
theta = phi + pi/2;          % angular distance from -90 degrees
rho   = R * theta;            % physical distance [m]

% Local Cartesian coordinates:
% +x = east, +y = north
x_grid = rho .* sin(lam);
y_grid = rho .* cos(lam);

% Convert source velocity components using the same rotation convention
u_east  =  -1*(u .* cosd(alpha) + v .* sind(alpha));
v_north = -u .* sind(alpha) + v .* cosd(alpha);

%% trim away invalid data
% Require all interpolated quantities to be finite
valid = isfinite(x_grid) & isfinite(y_grid) & ...
        isfinite(u_east) & isfinite(v_north) & ...
        isfinite(vorticity);

% % Determine the South Pole domain from valid samples
% xmin = min(x_grid(valid),[],'omitnan');
% xmax = max(x_grid(valid),[],'omitnan');
% ymin = min(y_grid(valid),[],'omitnan');
% ymax = max(y_grid(valid),[],'omitnan');

% hard coded adjustments to the above
% ymin = -4000e3; % conservative
ymin = -5300e3;
ymax = 4500e3;
xmin = -10000e3;
xmax = 11000e3;

% Apply the rectangular trimming bounds
gridmask = valid & ...
           x_grid >= xmin & x_grid <= xmax & ...
           y_grid >= ymin & y_grid <= ymax;

x_grid(~gridmask) = NaN;
y_grid(~gridmask) = NaN;
u_east(~gridmask) = NaN;
v_north(~gridmask) = NaN;
vorticity(~gridmask) = NaN;

%% interpolate onto grid
% create regular grid
L_cut = 1600e3; % filter scale
dx = 50e3;
dy = dx;

% Create regular grid
pad = L_cut;  % pad by 1 * L_cut (you can increase)
xgridfilt = (xmin-pad) : dx : (xmax+pad);
ygridfilt = (ymin-pad) : dy : (ymax+pad);
[Xg_pad, Yg_pad] = meshgrid(xgridfilt, ygridfilt);
[Ny, Nx] = size(Xg_pad);

% quantities to interpolate in vector form
xv = x_grid(:);
yv = y_grid(:);
uv = u_east(:);
vv = v_north(:);
vortv = vorticity(:);
valid = isfinite(uv) & isfinite(vv) & isfinite(vortv);

xv(~valid)=[];
yv(~valid)=[];
uv(~valid)=[];
vv(~valid)=[];
vortv(~valid)=[];

% Construct the mask for interpolation
Fdomain = scatteredInterpolant( ...
    xv, yv, ones(size(xv)), ...
    'natural', 'none');
insideMask_pad = isfinite(Fdomain(Xg_pad,Yg_pad));
extrapMask_pad = ~insideMask_pad;

% Interpolate data
Fu = scatteredInterpolant(xv, yv, uv,'natural','nearest');
Fv = scatteredInterpolant(xv, yv, vv, 'natural','nearest');
Fvort = scatteredInterpolant(xv, yv, vortv, 'natural','nearest');
Ugrid = Fu(Xg_pad, Yg_pad);
Vgrid = Fv(Xg_pad, Yg_pad);
Vortgrid = Fvort(Xg_pad, Yg_pad);

%% plot interpolated data
figure
set(gcf,'Position',[200,200,1200,1000])
tiledlayout(2,1,'TileSpacing','compact','Padding','compact')

nexttile
dat=Vortgrid/f;
dat(extrapMask_pad)=NaN;
% levels=[-1:0.1:1];
% contourf(Xg_pad, Yg_pad, dat,levels,'edgecolor','none');
% discreteColorbar(gca,levels,cmocean('balance'),'eastoutside','\zeta/f',4);
pcolor(Xg_pad, Yg_pad, dat,'linestyle','none');
cmocean('balance')
h=colorbar('eastoutside');
clim([-1,1]);
h.Label.String ='\zeta/f';
axis equal
xlim([xmin,xmax])
ylim([ymin,ymax])
xlabel('x (km)');
ylabel('y (km)');
title('(a) SP Rossby number');
set(gca,'fontsize',14);

nexttile
dat= sqrt(Ugrid.^2 + Vgrid.^2);
dat(extrapMask_pad)=NaN;
% levels = [-100:10:100];
% contourf(Xg_pad, Yg_pad, dat,levels,'edgecolor','none');
% discreteColorbar(gca,levels,cmocean('balance'),'eastoutside','m/s',2);
pcolor(Xg_pad, Yg_pad, dat,'linestyle','none');
cmocean('balance')
h=colorbar('eastoutside');
clim([-100,100]);
h.Label.String ='m/s';
axis equal
xlim([xmin,xmax])
ylim([ymin,ymax])
xlabel('x (km)');
ylabel('y (km)');
title('(b) SP Speed');
set(gca,'fontsize',14);

exportgraphics(gcf,'plots/interpolated_raw_S.png','resolution',500)

%% spectral filter on regular grid
% Hann taper
wx = hann(Nx);
wy = hann(Ny);
W = wy * wx.';          % Ny x Nx
taper_factor = 0.2;
W = (1 - taper_factor) + taper_factor * W;

Utaper = Ugrid .* W;
Vtaper = Vgrid .* W;

% Wavenumber grid
kx = (2*pi/(Nx*dx)) * ifftshift(-floor(Nx/2):ceil(Nx/2)-1);
ky = (2*pi/(Ny*dy)) * ifftshift(-floor(Ny/2):ceil(Ny/2)-1);
[KX, KY] = meshgrid(kx, ky);
K2 = KX.^2 + KY.^2;

k_cut = 2*pi / L_cut;
n = 4;   % Butterworth order

H = 1 ./ (1 + (sqrt(K2) ./ k_cut).^(2*n));

% FFT, filter velocity
Uhat = fft2(Utaper);
Vhat = fft2(Vtaper);

Uhatf = Uhat .* H;
Vhatf = Vhat .* H;

% --- vorticity and streamfunction in spectral space ---
% zeta_hat = i kx Vhatf - i ky Uhatf
vort_hat = 1i * (KX .* Vhatf - KY .* Uhatf);

% psi_hat = - zeta_hat / K^2  (handle k=0 separately)
K2(1,1) = Inf;                 % avoid division by zero
psi_hat = -vort_hat ./ K2;
psi_hat(1,1) = 0;              % set mean psi = 0 (arbitrary constant)


% get the nondivergent part of the velocity from psi
Upsi_hat = -1i * KY .* psi_hat;   % U = -i ky psi_hat
Vpsi_hat =  1i * KX .* psi_hat;   % V =  i kx psi_hat

% back to physical space (still tapered)
Uf_taper_obs     = real(ifft2(Uhatf));
Vf_taper_obs      = real(ifft2(Vhatf));
Uf_taper      = real(ifft2(Upsi_hat));   % nondivergent velocity
Vf_taper      = real(ifft2(Vpsi_hat));   % nondivergent velocity
vort_taper    = real(ifft2(vort_hat));
psi_taper     = real(ifft2(psi_hat));

% Remove taper
epsW = 1e-10;
u_specfilt_pad         = Uf_taper      ./ (W + epsW);
v_specfilt_pad         = Vf_taper      ./ (W + epsW);
u_specfilt_pad_obs         = Uf_taper_obs      ./ (W + epsW);
v_specfilt_pad_obs         = Vf_taper_obs      ./ (W + epsW);
vort_specfilt_calc_pad = vort_taper    ./ (W + epsW);
psi_specfilt_calc_pad  = psi_taper     ./ (W + epsW);

% Remove values located outside the interpolation domain
u_specfilt_pad(extrapMask_pad)         = NaN;
v_specfilt_pad(extrapMask_pad)         = NaN;
u_specfilt_pad_obs(extrapMask_pad)     = NaN;
v_specfilt_pad_obs(extrapMask_pad)     = NaN;
vort_specfilt_calc_pad(extrapMask_pad) = NaN;
psi_specfilt_calc_pad(extrapMask_pad)  = NaN;

% remove padding
ix = (xgridfilt >= xmin) & (xgridfilt <= xmax);
iy = (ygridfilt >= ymin) & (ygridfilt <= ymax);

Xg = Xg_pad(iy, ix);
Yg = Yg_pad(iy, ix);
u_specfilt    = u_specfilt_pad(iy, ix);
v_specfilt    = v_specfilt_pad(iy, ix);
u_specfilt_obs    = u_specfilt_pad_obs(iy, ix);
v_specfilt_obs    = v_specfilt_pad_obs(iy, ix);
vort_specfilt_calc = vort_specfilt_calc_pad(iy, ix);
psi_specfilt_calc  = psi_specfilt_calc_pad(iy, ix);

save('filtered_data_south.mat', ...
    'Xg', 'Yg', ...        % trimmed grid
    'u_specfilt', 'v_specfilt', ...          % filtered nondivergent velocities
    'u_specfilt_obs', 'v_specfilt_obs',...  % filtered observed velocities
    'vort_specfilt_calc', 'psi_specfilt_calc', ...     % filtered vorticity & streamfunction
    '-v7.3');


%% plot filtered data
% calculate the speed
speed = sqrt(u_specfilt.^2 + v_specfilt.^2);

figure
set(gcf,'Position',[200,200,1000,800])
tiledlayout(2,1,'TileSpacing','compact','Padding','compact')

nexttile
hold on
levels=[-1:0.1:1];
contourf(Xg, Yg,vort_specfilt_calc/f,levels,'edgecolor','none');
contour(Xg,Yg,psi_specfilt_calc,'linewidth',1,'linecolor','k')
discreteColorbar(gca,levels,cmocean('balance'),'eastoutside','\zeta/f',4);
axis equal
% xlim([xmin,xmax])
% ylim([ymin,ymax])
xlabel('x (km)');
ylabel('y (km)');
title('Rossby number');
set(gca,'fontsize',14);

nexttile
hold on
levels = [-100:10:100];
contourf(Xg, Yg, speed,levels,'edgecolor','none');
contour(Xg,Yg,psi_specfilt_calc,'linewidth',1,'linecolor','k')
discreteColorbar(gca,levels,cmocean('balance'),'eastoutside','m/s',2);
axis equal
% xlim([xmin,xmax])
% ylim([ymin,ymax])
xlabel('x (km)');
ylabel('y (km)');
title('Speed');
set(gca,'fontsize',14);

exportgraphics(gcf,'plots/filtered_S.png','resolution',500)


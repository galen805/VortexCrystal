% contour dynamics main script
% simulates vortex patch interactions on an f-plane
close all; clear; clc;

%% setup
% runname
runName = 'testStable';
figName = strcat(runName,'.png');
videoName = strcat(runName,'.mp4');

% switches / flags
writeVideoOn   = 1;  % <--- set true to show/write MP4
saveStatesOn   = 1;   % <--- save (x,y,t) snapshots every plotfreq
makeStaticPlot = 1;   % <--- after run, make static alpha-time plot

% Runtime parameters
Tf = 10;               % dimensional time in days
plotfreq = 5;        % save state every plotfreq steps
stretchFactorMax = 5;  % stop if any segment exceeds this * its initial length

%% initialize vortex patches
% evenly distributed ring vortices

% STABLE EXAMPLE
n_vort = 5;      % number on the ring
N = 100;         % nodes per vortex patch
R = 1;           % initial displacement
Gam = 0.5;       % L/Ld
r0 = 0.3;       % vortex radius
Qvort = 1;       % vortex strength

% % UNSTABLE EXAMPLE
% n_vort = 5;      % number on the ring
% N = 100;         % nodes per vortex patch
% R = 1;           % initial displacement
% Gam = 1;         % L/Ld
% r0 = 0.5;        % vortex radius
% Qvort = 1;       % vortex strength

%
theta = linspace(0, 2*pi, N+1); theta(end) = [];
numPatch = n_vort;
x = zeros(numPatch, N);
y = zeros(numPatch, N);
Q = Qvort * ones(numPatch, 1);

for p = 1:n_vort
    phi0 = (p-1) * 2*pi / n_vort;
    xc = R * cos(phi0);
    yc = R * sin(phi0);

    x(p,:) = xc + r0*cos(theta);
    y(p,:) = yc + r0*sin(theta);
end


%% bookkeeping
% record initial max segment length (per patch)
maxSeg0 = zeros(numPatch,1);
minSeg0 = zeros(numPatch,1);
for j = 1:numPatch
    xj = x(j,:); yj = y(j,:);
    xjp = [xj, xj(1)];
    yjp = [yj, yj(1)];
    segLen = hypot(diff(xjp), diff(yjp));   % N segments
    maxSeg0(j) = max(segLen);
    minSeg0(j) = min(segLen);
end

% pass parameters for velocity field calculation
params.numPatch = numPatch;
params.N        = N;
params.Gam      = Gam;
params.Q        = Q;

%% set up plotting (ONLY if we will show preview or write video)
% performance preference: no live figure updates unless video is being written
if writeVideoOn
    showPreviewOn  = true;
else
    showPreviewOn = false;
end

if showPreviewOn || writeVideoOn
    figure(1);
    set(gcf,'pos',[200,200,600,600])
    axis equal
    xlim([-R, R]*2)
    ylim([-R, R]*2)
    grid on
    box on
    hold on
    set(gca,'fontsize',16)
end

%% video writer (optional)
if writeVideoOn
    vobj = VideoWriter(videoName, 'MPEG-4');
    vobj.FrameRate = 24;
    open(vobj);
end

%% main loop

% estimate timestep
[uTmp, vTmp] = velocity_field(x, y, params);
umax = max(hypot(uTmp(:), vTmp(:)));
dsmin = min(minSeg0);
dt = 0.2 * dsmin / umax; % CFL

t  = 0:dt:Tf;
numSteps = length(t);

% STATE STORAGE (every plotfreq)
if saveStatesOn
    snapIdx = 1:plotfreq:numSteps;
    nsnap   = numel(snapIdx);
    x_snap  = NaN(numPatch, N, nsnap);
    y_snap  = NaN(numPatch, N, nsnap);
    t_snap  = NaN(1, nsnap);        % nondimensional time
    tdim_snap = NaN(1, nsnap);      % days
    sCount  = 0;
end

fprintf('starting run: dt=%.2f, Tf=%d, numSteps=%d\n\n\n',dt,Tf,numSteps)
for k = 1:numSteps

    % diagnostics / saving / optional preview & video
    if rem(k-1, plotfreq) == 0

        % save state (before any break)
        if saveStatesOn
            sCount = sCount + 1;
            x_snap(:,:,sCount) = x;
            y_snap(:,:,sCount) = y;
            t_snap(sCount)     = t(k);
        end

        % track area
        Atot = 0;
        for i = 1:numPatch
            Atot = Atot + polyarea([x(i,:), x(i,1)], [y(i,:), y(i,1)]);
        end
        fprintf('step %d/%d, t = %.2f, Atot = %.6g \n',k, numSteps, t(k), Atot);

        % stretch safety check
        maxSegNow = zeros(numPatch,1);
        for j = 1:numPatch
            xj = x(j,:); yj = y(j,:);
            xjp = [xj, xj(1)];
            yjp = [yj, yj(1)];
            segLen = hypot(diff(xjp), diff(yjp));
            maxSegNow(j) = max(segLen);
        end

        stretchRatio = max(maxSegNow ./ maxSeg0);
        if stretchRatio > stretchFactorMax
            warning('Stopping: max segment stretched by %.3g (> %.3g).', ...
                stretchRatio, stretchFactorMax);
            break
        end

        % optional live preview (disabled when writeVideoOn == false)
        if showPreviewOn
            cla
            for i = 1:numPatch
                plot([x(i,:), x(i,1)], [y(i,:), y(i,1)], '-', 'LineWidth', 2);
            end
            title(sprintf('t = %.2f', t(k)));
            axis equal
            xlim([-R, R]*2)
            ylim([-R, R]*2)
            drawnow
        end

        % optional video frame (requires a figure)
        if writeVideoOn
            if ~showPreviewOn
                % draw without showing frequent UI updates: still need a figure to capture frames
                clf(1);
                figure(1);
                hold on; box on; grid on; axis equal
                xlim([-R, R]*2)
                ylim([-R, R]*2)
                for i = 1:numPatch
                    plot([x(i,:), x(i,1)], [y(i,:), y(i,1)], '-', 'LineWidth', 2);
                end
                title(sprintf('t = %.2f', t(k)));
                drawnow
            end
            writeVideo(vobj, getframe(gcf));
        end
    end

    % RK2 (midpoint) step: X_{n+1} = X_n + dt * F(X_n + 0.5 dt F(X_n))
    [u1, v1] = velocity_field(x, y, params);            % stage 1
    xmid = x + 0.5*dt*u1;
    ymid = y + 0.5*dt*v1;

    [u2, v2] = velocity_field(xmid, ymid, params);      % stage 2
    x = x + dt*u2;
    y = y + dt*v2;

end

% trim snapshots if we broke early
if saveStatesOn
    x_snap = x_snap(:,:,1:sCount);
    y_snap = y_snap(:,:,1:sCount);
    t_snap = t_snap(1:sCount);
end

% close video if used
if writeVideoOn
    close(vobj);
    fprintf('Video written to %s\n', videoName);
end

% close preview figure if it exists
if showPreviewOn || writeVideoOn
    close(1);
end

%% STATIC SUMMARY PLOT:
% centroids + large direction arrow + initial marker
if saveStatesOn && makeStaticPlot

    nsnap = size(x_snap,3);
    if nsnap < 2
        warning('Not enough snapshots to make static plot.');
    else
        cmap = lines(numPatch);
        patchColors = cmap;


        % compute centroids
        xc = NaN(numPatch, nsnap);
        yc = NaN(numPatch, nsnap);
        for s = 1:nsnap
            for i = 1:numPatch
                xc(i,s) = mean(x_snap(i,:,s), 'omitnan');
                yc(i,s) = mean(y_snap(i,:,s), 'omitnan');
            end
        end

        % keyframes for contours
        nKey = min(5, nsnap);
        keyIdx = unique(round(linspace(1, nsnap, nKey)));

        figure;
        set(gcf,'pos',[200,200,600,600])
        hold on; box on; grid on; axis equal
        xlim([-R, R]*2)
        ylim([-R, R]*2)
        title('Vortex patch evolution');
        set(gca,'fontsize',16);

        % keyframe contours (thin, semi-transparent)
        for ii = 1:numel(keyIdx)
            s = keyIdx(ii);
            a = 0.1 + 0.9*(ii-1)/max(1,(numel(keyIdx)-1));
            for i = 1:numPatch
                xx = [x_snap(i,:,s), x_snap(i,1,s)];
                yy = [y_snap(i,:,s), y_snap(i,1,s)];
                h = plot(xx, yy, '-', 'LineWidth', 1);
                h.Color = [patchColors(i,:), a];

                if ii==numel(keyIdx)
                    h.LineWidth = 2;
                end
            end
        end

        % centroid trajectories
        for i = 1:numPatch
            plot(xc(i,:), yc(i,:), '-', ...
                'LineWidth', 4, 'Color', patchColors(i,:));
        end

        % final state
        for i = 1:numPatch
            plot(xc(i,end), yc(i,end), '.', ...
                'MarkerSize', 20, ...
                'LineWidth', 4, ...
                'Color', patchColors(i,:), ...
                'MarkerFaceColor', 'none');
        end

        % annotation
        text(0.02,0.98, ...
            sprintf('t = %.2f → %.2f', ...
            t_snap(1), t_snap(end)), ...
            'Units','normalized', ...
            'HorizontalAlignment','left', ...
            'VerticalAlignment','top', ...
            'FontSize',14,'BackgroundColor','none');

        set(gcf,'renderer','painters')
        exportgraphics(gcf,figName,'resolution',300)
    end
end

%% save snapshots to .mat for later analysis
if saveStatesOn
    save('contour_snapshots.mat','x_snap','y_snap','t_snap','numPatch','N','plotfreq');
end


%% velocity field function
function [u, v] = velocity_field(x, y, params)
% Vectorized velocity evaluation for contour dynamics
%
% x, y : (numPatch x N)
% u, v : (numPatch x N)

numPatch = params.numPatch;
N        = params.N;
Gam      = params.Gam;
Q        = params.Q(:);              % force column

% Build source segment midpoints and segment vectors
xm = zeros(numPatch, N);
ym = zeros(numPatch, N);
dX = zeros(numPatch, N);
dY = zeros(numPatch, N);

for j = 1:numPatch
    xj = x(j,:);
    yj = y(j,:);
    xjp = [xj, xj(1)];
    yjp = [yj, yj(1)];

    xA = xjp(1:N);
    yA = yjp(1:N);
    xB = xjp(2:N+1);
    yB = yjp(2:N+1);

    xm(j,:) = 0.5*(xA + xB);
    ym(j,:) = 0.5*(yA + yB);
    dX(j,:) = (xB - xA);
    dY(j,:) = (yB - yA);
end

% Flatten sources (S = numPatch*N)
Xs  = reshape(xm, 1, []);
Ys  = reshape(ym, 1, []);
dXs = reshape(dX, [], 1);   % Sx1
dYs = reshape(dY, [], 1);   % Sx1

% Repeat patch strengths per segment (Sx1)
Qseg = repelem(Q, N);

% Segment weights (Sx1)
wX = (Qseg/(2*pi)) .* dXs;
wY = (Qseg/(2*pi)) .* dYs;

% Flatten targets (M = numPatch*N)
Xt = reshape(x, [], 1);     % Mx1
Yt = reshape(y, [], 1);     % Mx1

% Kernel evaluation: distances (MxS)
DX = Xt - Xs;               % MxS (implicit expansion)
DY = Yt - Ys;               % MxS
Rf = sqrt(DX.^2 + DY.^2);

K0 = besselk(0, Gam * Rf);   % MxS

% Induced velocities (Mx1)
Ut = K0 * wX;               % (MxS)*(Sx1)
Vt = K0 * wY;

% Reshape back to (numPatch x N)
u = reshape(Ut, numPatch, N);
v = reshape(Vt, numPatch, N);

end

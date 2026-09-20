function cb = discreteColorbar(ax, levels, cmap, location, label, labelFreq)
%DISCRETECOLORBAR  Discrete colorbar tied to specific axes
%
%   cb = discreteColorbar(ax, levels, cmap, location, labelFreq)
%
%   ax        : target axes handle
%   levels    : monotonic vector of bin edges
%   cmap      : Mx3 colormap array (do not need to specify M)
%   location  : colorbar location
%   label     : colorbar label string
%   labelFreq : label every Nth edge (default 1)

% you can then use the handle to set further properties


% usage example
% levels=-10:1:10;
% contourf(peaks,levels)
% cb=discreteColorbar(gca,levels, cmocean('thermal'),'eastoutside','units',2);
% set(gca,'fontsize',20)


n = numel(levels) - 1;

% Resample to one color per bin
cmap = interp1(linspace(0,1,size(cmap,1)), cmap, linspace(0,1,n));
colormap(ax, cmap)

% Apply to specific axes
colormap(ax, cmap)

% axis limits
clim(ax, [levels(1) levels(end)])

% Create colorbar
cb = colorbar(ax, 'Location', location);
cb.Label.String = label;

% add ticks
cb.Ticks = levels;

% Label control
labels = repmat({''}, size(levels));
idx = 1:labelFreq:numel(levels);
labels(idx) = arrayfun(@(x) sprintf('%g',x), levels(idx), ...
                       'UniformOutput', false);

cb.TickLabels = labels;
cb.Ruler.TickLabelRotation=0;

end



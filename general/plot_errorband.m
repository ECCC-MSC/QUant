function plot_errorband(xx, yy, st_err, plot_args, use_transparency, alpha)
%PLOT_ERRORBAND PLOT a line with an error polygon made with FILL
%
%     plot_errorband(xx, yy, st_err, plot_args, use_transparency, alpha)
%
% This function is useful for plots where the x-axis is densely sampled and many
% error-bars would look ugly. Instead a blob of error goop is drawn behind the
% line. This is colored fainter than the line and is optionally transparent.
%
% Inputs:
%                xx Nx1 or 1xN   x-data
%                yy Nx1 or 1xN   y-data
%            st_err Nx1 or 1xN   half-width of error band, or cell array of two
%                                such vectors for upper band and lower band
%         plot_args cell/string  Optional arg(s) passed to PLOT. Often a color.
%  use_transparency    1x1       Default false, band color is mixed with
%                                white, but is opaque
%             alpha    1x1       Default 0.3. How dark/opaque the band is.

% Why the ugly cell array thing in st_err instead of Nx2 arrays? Because what if
% I input a 2x2 array? I don't want to really force the user to insert vectors
% the right way around, but I don't like interfaces that are ambiguous.

% This function doesn't return a handle. More thought and work would be required
% if one wanted to interact with the resulting line+band object. I'm not sure
% what the best way of doing this would be. If one were going all out, then the
% transparency and alpha options should probably be specified in the same way as
% MATLAB's plot arguments.

% Iain Murray, June 2009

% Option processing, sanitization and setting of defaults
if ~exist('plot_args', 'var')
    plot_args = {};
else
    if ~iscell(plot_args)
        plot_args = {plot_args};
    end
end
xx = xx(:)';
yy = yy(:)';
if iscell(st_err)
    up = st_err{1}(:)';
    down = st_err{2}(:)';
else
    up = st_err(:)';
    down = st_err(:)';
end
if ~exist('alpha', 'var')
    alpha = 0.3;
end

% Plot the line
hh = plot(xx, yy, plot_args{:});

% Grab color from plot (allows user not to specify color)
% and sort out the color for the patch
plot_color = get(hh, 'Color');
if exist('use_transparency', 'var') && use_transparency
    opts = {'EdgeAlpha', alpha, 'FaceAlpha', alpha};
    band_color = plot_color;
else
    white = ones(size(plot_color));
    band_color = alpha*plot_color + (1-alpha)*white;
    opts = {};
end

% Here is the contour to fill
xbits = [xx, flipdim(xx, 2)];
ybits = [yy + up, flipdim(yy - down, 2)];

% When actually filling do a dance that will maintain the user's HOLD state
% Note remembering with ISHOLD isn't good enough as it doesn't distinguish
% between "hold on" and "hold all"
old_NextPlot = get(gcf, 'NextPlot');
set(gcf, 'NextPlot', 'add');
fill(xbits , ybits, band_color, 'EdgeColor', band_color, opts{:});
set(gcf, 'NextPlot', old_NextPlot);

% Put the line back above the fill. We drew it first to grab the color, which is
% simpler to code than working out what the color would be without drawing.
uistack(hh);
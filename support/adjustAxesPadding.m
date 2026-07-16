function adjustAxesPadding(ax, padLeft, padRight, padBottom, padTop)
%adjustAxesPadding Expand axes limits with independent padding for each side.
%
%   adjustAxesPadding(ax, padLeft, padRight, padBottom, padTop)
%
%   Each padding value is a *fraction* of the total data range.
%   Examples:
%       adjustAxesPadding(gca, 0.05, 0.05, 0.05, 0.05);   % symmetric 5%
%       adjustAxesPadding(gca, 0.10, 0.02, 0.08, 0.15);   % per-side padding

    if nargin < 5
        error('Need 5 arguments: ax, padLeft, padRight, padBottom, padTop.');
    end

    ch = ax.Children;

    % Initialize bounds
    xmin = inf; xmax = -inf;
    ymin = inf; ymax = -inf;

    for k = 1:numel(ch)
        obj = ch(k);

        % LINE, SCATTER, etc (XData / YData)
        if isprop(obj, 'XData') && isprop(obj, 'YData')
            x = obj.XData(:);
            y = obj.YData(:);
            xmin = min(xmin, min(x));
            xmax = max(xmax, max(x));
            ymin = min(ymin, min(y));
            ymax = max(ymax, max(y));
        end

        % TEXT objects (use extent box)
        if isa(obj, 'matlab.graphics.primitive.Text')
            ext = obj.Extent; % [x, y, width, height]
            xmin = min(xmin, ext(1));
            xmax = max(xmax, ext(1) + ext(3));
            ymin = min(ymin, ext(2));
            ymax = max(ymax, ext(2) + ext(4));
        end

        % PATCH (Vertices)
        if isprop(obj, 'Vertices')
            verts = obj.Vertices;
            if ~isempty(verts)
                xmin = min(xmin, min(verts(:,1)));
                xmax = max(xmax, max(verts(:,1)));
                ymin = min(ymin, min(verts(:,2)));
                ymax = max(ymax, max(verts(:,2)));
            end
        end
    end

    % Compute ranges
    dx = xmax - xmin;
    dy = ymax - ymin;

    % Apply per-side padding
    newXmin = xmin - padLeft;
    newXmax = xmax + padRight;
    newYmin = ymin - padBottom;
    newYmax = ymax + padTop;
    % REVIEW PROPOSAL MAIN-09 (ADD): uncomment to make padding values match
    % the documented fractional interpretation.
    % newXmin=xmin-padLeft*dx;   newXmax=xmax+padRight*dx;
    % newYmin=ymin-padBottom*dy; newYmax=ymax+padTop*dy;

    % Apply limits
    ax.XLim = [newXmin newXmax];
    ax.YLim = [newYmin newYmax];
end

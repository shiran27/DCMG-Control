function saveFigureHighRes(baseName, varargin)
%SAVEFIGUREHIGHRES Save a figure as high-res PNG and FIG, with padding.
%
%   saveFigureHighRes('myPlot')
%   saveFigureHighRes('myPlot','Width',6,'Height',4,'FontSize',14,'DPI',300)
%   saveFigureHighRes('myPlot','Padding',20)   % padding in points

    % ---- Defaults ----
    defaults.Figure   = gcf;
    defaults.Width    = 6;           % in Units
    defaults.Height   = 4;           % in Units
    defaults.Units    = 'inches';    % 'inches','centimeters',...
    defaults.DPI      = 300;
    defaults.FontSize = 14;
    defaults.Padding  = 20;          % points around content

    % ---- Parse inputs ----
    p = inputParser;
    addRequired(p,'baseName',@(x)ischar(x) || isstring(x));
    addParameter(p,'Figure',   defaults.Figure);
    addParameter(p,'Width',    defaults.Width,  @isnumeric);
    addParameter(p,'Height',   defaults.Height, @isnumeric);
    addParameter(p,'Units',    defaults.Units,  @ischar);
    addParameter(p,'DPI',      defaults.DPI,    @isnumeric);
    addParameter(p,'FontSize', defaults.FontSize,@isnumeric);
    addParameter(p,'Padding',  defaults.Padding,@isnumeric);

    parse(p,baseName,varargin{:});
    R = p.Results;

    fig = R.Figure;

    % ---- Set figure size on screen ----
    set(fig,'Units',R.Units);
    pos = get(fig,'Position');
    pos(3) = R.Width;
    pos(4) = R.Height;
    set(fig,'Position',pos);

    % % ---- Font size for all text/axes ----
    % if ~isempty(R.FontSize)
    %     objs = findall(fig,'-property','FontSize');
    %     set(objs,'FontSize',R.FontSize);
    % end

    % ---- File names ----
    baseName = char(baseName);
    % REVIEW PROPOSAL MAIN-08 (ADD): allow a clean checkout to save figures.
    % if ~isfolder('Results'), mkdir('Results'); end
    pngName  = ['Results/', baseName '.png'];
    figName  = ['Results/', baseName '.fig'];

    % ---- Save high-res PNG (with padding so labels aren’t cut) ----
    % Works in R2020a+; for older versions we can fall back to print.
    try
        exportgraphics(fig, pngName, ...
            'Resolution',     R.DPI, ...
            'ContentType',    'image', ...
            'BackgroundColor','white', ...
            'Padding',        R.Padding);
    catch
        % Fallback for older MATLAB
        set(fig,'PaperUnits',R.Units, ...
                'PaperPosition',[0 0 R.Width R.Height], ...
                'PaperPositionMode','manual');
        print(fig, pngName, '-dpng', sprintf('-r%d',R.DPI));
    end

    % ---- Save FIG ----
    savefig(fig, figName);

    fprintf('Saved figure as: %s and %s\n', pngName, figName);
end

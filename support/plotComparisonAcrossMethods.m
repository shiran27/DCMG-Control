function plotComparisonAcrossMethods(Results,axisLimits)
%plotComparisonAcrossMethods  Compare mean per-unit trajectories across methods.
%   Results is a struct array with fields:
%     Results(k).method
%     Results(k).t
%     Results(k).avg.V_pu, Results(k).avg.It_pu, Results(k).avg.I_pu, Results(k).avg.Vt_pu
%
%   (Optionally) Results(k).std.* can exist, but not used here.

    if nargin < 1 || isempty(Results)
        warning('plotComparisonAcrossMethods:Empty', 'Results is empty. Nothing to plot.');
        return;
    end

    % --- Sort by method name (optional but makes legend stable) ---
    % [~,ord] = sort(string({Results.method}));
    % Results = Results(ord);

    % --- Event lines (same logic you used: fractions of t(end)) ---
    t_f = Results(1).t(end);
    % REVIEW PROPOSAL MAIN-10 (ADD): comparison requires compatible horizons.
    % assert(all(arrayfun(@(r) abs(r.t(end)-t_f)<=1e-12,Results)), ...
    %     'main:IncompatibleResults','Compared methods use different horizons.');
    eventTimes = [0.05, 0.2, 0.4, 0.6, 0.65, 0.8] * t_f;

    % -------------------- Figure 1: States --------------------
    figure('Name','Comparison: Avg Per-Unit States','Color','w');

    % (1) Avg V
    ax1 = subplot(2,1,1); hold(ax1,'on'); grid(ax1,'on');
    for k = 1:numel(Results)
        plot(ax1, Results(k).t, Results(k).avg.V_pu, 'LineWidth', 1.4, ...
            'DisplayName', Results(k).method);
    end
    yline(ax1, 1, 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
    for te = eventTimes
        xline(ax1, te, 'k--', 'LineWidth', 0.75, 'HandleVisibility', 'off');
    end
    ylabel(ax1, '$\overline{V}$ [pu]', 'Interpreter','latex');
    legend(ax1, 'show', 'Location','south');  % change to 'best' if you prefer
    axis(axisLimits(1,:));
    % REVIEW PROPOSAL MAIN-11 (ADD): explicitly target the intended axes.
    % axis(ax1,axisLimits(1,:));

    % (2) Avg It
    ax2 = subplot(2,1,2); hold(ax2,'on'); grid(ax2,'on');
    for k = 1:numel(Results)
        plot(ax2, Results(k).t, Results(k).avg.It_pu, 'LineWidth', 1.4, ...
            'DisplayName', Results(k).method);
    end
    yline(ax2, 1, 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
    for te = eventTimes
        xline(ax2, te, 'k--', 'LineWidth', 0.75, 'HandleVisibility', 'off');
    end
    ylabel(ax2, '$\overline{I_t}$ [pu]', 'Interpreter','latex');
    xlabel(ax2, 'Time [s]');
    legend(ax2, 'show', 'Location','south');
    axis(axisLimits(2,:));
    % axis(ax2,axisLimits(2,:)); % REVIEW PROPOSAL MAIN-11

    % -------------------- Figure 2: Inputs --------------------
    figure('Name','Comparison: Avg Per-Unit Inputs','Color','w');

    % (1) Avg I
    ax3 = subplot(2,1,1); hold(ax3,'on'); grid(ax3,'on');
    for k = 1:numel(Results)
        plot(ax3, Results(k).t, Results(k).avg.I_pu, 'LineWidth', 1.4, ...
            'DisplayName', Results(k).method);
    end
    yline(ax3, 1, 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
    for te = eventTimes
        xline(ax3, te, 'k--', 'LineWidth', 0.75, 'HandleVisibility', 'off');
    end
    ylabel(ax3, '$\overline{I}$ [pu]', 'Interpreter','latex');
    legend(ax3, 'show', 'Location','south');
    axis(axisLimits(3,:));
    % axis(ax3,axisLimits(3,:)); % REVIEW PROPOSAL MAIN-11

    % (2) Avg Vt
    ax4 = subplot(2,1,2); hold(ax4,'on'); grid(ax4,'on');
    for k = 1:numel(Results)
        plot(ax4, Results(k).t, Results(k).avg.Vt_pu, 'LineWidth', 1.4, ...
            'DisplayName', Results(k).method);
    end
    yline(ax4, 1, 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
    for te = eventTimes
        xline(ax4, te, 'k--', 'LineWidth', 0.75, 'HandleVisibility', 'off');
    end
    ylabel(ax4, '$\overline{V_t}$ [pu]', 'Interpreter','latex');
    xlabel(ax4, 'Time [s]');
    legend(ax4, 'show', 'Location','south');
    axis(axisLimits(4,:));
    % axis(ax4,axisLimits(4,:)); % REVIEW PROPOSAL MAIN-11


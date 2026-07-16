function trajDataPU = drawStateTrajectories(t,X,U,Xs,Us,N,saveImages,figTitle)
%DRAWSTATETRAJECTORIES Plot and return per-unit DCMG states and inputs.
%   trajDataPU = drawStateTrajectories(...) separates the stacked state and
%   input trajectories, normalizes them by their operating points, and
%   optionally saves state and input figures under Results/.

    
    t_f = t(end);
    eventTimes = [];
    for eventTime = [0.05, 0.2, 0.4, 0.6, 0.65, 0.8]
        eventTimes = [eventTimes, eventTime*t_f];
    end
    
    It = zeros(numel(t), N);
    V = zeros(numel(t), N);
    I = zeros(numel(t), N);
    Vt = zeros(numel(t), N);

    Its = zeros(numel(t), N);
    Vs = zeros(numel(t), N);
    Is = zeros(numel(t), N);
    Vts = zeros(numel(t), N);
    
    for k = 1:N
        V(:,k) = X(:, 2*k-1);
        It(:,k) = X(:, 2*k);
        I(:,k) = U(:, 2*k-1);
        Vt(:,k) = U(:, 2*k);

        Vs(:,k) = Xs(:, 2*k-1);
        Its(:,k) = Xs(:, 2*k);
        Is(:,k) = Us(:, 2*k-1);
        Vts(:,k) = Us(:, 2*k);
    end

    % Per-Unit Values
    V_pu = V ./ Vs;
    It_pu = It ./ Its;
    I_pu = I ./ Is;
    % REVIEW PROPOSAL MAIN-07 (ADD): line-current equilibrium values may be
    % zero or very small. This guarded override prevents Inf/NaN plots. A DG
    % rated-current base is preferable if absolute per-unit meaning is needed.
    % Ibase = sign(Is).*max(abs(Is),1e-6); Ibase(Ibase==0)=1e-6;
    % I_pu = I./Ibase;
    Vt_pu = Vt ./ Vts;

    % Data storage
    trajDataPU.V_pu = V_pu;
    trajDataPU.It_pu = It_pu;
    trajDataPU.I_pu = I_pu;
    trajDataPU.Vt_pu = Vt_pu;

    % Plotting state and input trajectories
    figure('Name','Per-Unit States over time');
    subplot(2,1,1);
    plot(t, V_pu, 'LineWidth', 1.25);
    yline(1,'k--','LineWidth',0.8); grid on;
    ylabel('$V_i\ $ [pu]',Interpreter='latex'); 
    hold on
    for te = eventTimes
        xline(te, 'k--', 'LineWidth', 0.75);
    end

    subplot(2,1,2);
    plot(t, It_pu, 'LineWidth', 1.25);
    yline(1,'k--','LineWidth',0.8); grid on;
    ylabel('$I_{ti}\ $ [pu]',Interpreter='latex'); xlabel('Time [s]'); 
    % legend(arrayfun(@(k) sprintf('DG %d',k), 1:N, 'UniformOutput',false), Location="southoutside", Orientation="horizontal");
    hold on
    for te = eventTimes
        xline(te, 'k--', 'LineWidth', 0.75);
    end
    

    if saveImages
        saveFigureHighRes([figTitle,'-State'], ...
            'Width', 6, 'Height', 6, 'Units', 'inches', ...
            'DPI', 600, 'FontSize', 8);
    end
    
    figure('Name','Per-Unit Controls over time');
    subplot(2,1,1);
    plot(t, I_pu, 'LineWidth', 1.25);
    yline(1,'k--','LineWidth',0.8); grid on;
    ylabel('$I_i\ $ [pu]',Interpreter='latex'); 
    hold on
    for te = eventTimes
        xline(te, 'k--', 'LineWidth', 0.75);
    end
    axis([-inf,inf,-1,3])
    

    subplot(2,1,2);
    plot(t, Vt_pu, 'LineWidth', 1.25);
    yline(1,'k--','LineWidth',0.8); grid on;
    ylabel('$V_{ti}\ $ [pu]',Interpreter='latex'); xlabel('Time [s]'); 
    % legend(arrayfun(@(k) sprintf('DG %d',k), 1:N, 'UniformOutput',false), Location="southoutside", Orientation="horizontal");
    hold on
    for te = eventTimes
        xline(te, 'k--', 'LineWidth', 0.75);
    end
    

    if saveImages
        saveFigureHighRes([figTitle,'-Input'], ...
            'Width', 6, 'Height', 6, 'Units', 'inches', ...
            'DPI', 600, 'FontSize', 8);
    end
end

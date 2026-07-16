function Results = appendResult(Results, methodName, t, trajDataPU)
%APPENDRESULT Add or replace one named trajectory result.
    R = struct();
    R.method = char(methodName);
    R.t = t(:);

    % Raw per-unit (Nt x N)
    R.raw = trajDataPU;

    % Mean/std across DGs
    R.avg.V_pu  = mean(trajDataPU.V_pu,  2);
    R.avg.It_pu = mean(trajDataPU.It_pu, 2);
    R.avg.I_pu  = mean(trajDataPU.I_pu,  2);
    R.avg.Vt_pu = mean(trajDataPU.Vt_pu, 2);

    R.std.V_pu  = std(trajDataPU.V_pu,  0, 2);
    R.std.It_pu = std(trajDataPU.It_pu, 0, 2);
    R.std.I_pu  = std(trajDataPU.I_pu,  0, 2);
    R.std.Vt_pu = std(trajDataPU.Vt_pu, 0, 2);

    % Overwrite if same method already exists (useful when re-running one section)
    if ~isempty(Results)
        idx = find(strcmp({Results.method}, R.method), 1);
        if ~isempty(idx)
            Results(idx) = R;
            return;
        end
    end

    Results = [Results; R];
end

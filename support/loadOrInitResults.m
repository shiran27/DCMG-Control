function Results = loadOrInitResults(resultsFile, resetResults)
%LOADORINITRESULTS Load an existing result array or return an empty one.
    if nargin < 2, resetResults = 0; end

    if resetResults
        Results = struct([]);
        return;
    end

    if exist(resultsFile, "file")
        S = load(resultsFile, "Results");
        if isfield(S, "Results")
            Results = S.Results;
        else
            Results = struct([]);
        end
    else
        Results = struct([]);
    end
end

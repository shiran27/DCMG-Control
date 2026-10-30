function ops = dcmgSolverOptions()
%DCMGSOLVEROPTIONS Use an available, tested SDP solver without global path edits.
persistent solver
if isempty(solver)
    root = fileparts(fileparts(mfilename('fullpath')));
    if isfolder(fullfile(root,'vendor','sedumi'))
        addpath(fullfile(root,'vendor','sedumi'));
    end
    assert(exist('sdpvar','file')==2,'DCMG:MissingYALMIP', ...
        'Add YALMIP to the MATLAB path before running this experiment.');
    candidates = {'mosek','sedumi','sdpt3'};
    for k=1:numel(candidates)
        v=sdpvar(2,2,'symmetric');
        trial=optimize(v>=eye(2),trace(v),sdpsettings( ...
            'solver',candidates{k},'verbose',0,'warning',0));
        if trial.problem==0
            solver=candidates{k};
            break
        end
    end
    assert(~isempty(solver),'DCMG:NoSDPSolver', ...
        'No working SDP solver. See vendor/README.md for SeDuMi installation.');
    fprintf('DCMG SDP solver: %s\n',solver);
end
ops=sdpsettings('solver',solver,'verbose',0,'warning',0);
if strcmp(solver,'sedumi')
    ops.sedumi.eps=1e-9;
    ops.sedumi.maxiter=200;
end
end

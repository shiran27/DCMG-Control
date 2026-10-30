function report = run_all_paper_results(options)
%RUN_ALL_PAPER_RESULTS Rebuild proposed methods, both baselines, and checks.
%   report = run_all_paper_results;
% main.mlx calls this complete reproduction entry point.
arguments
    options.MakePlots (1,1) logical = true
end
root=fileparts(mfilename('fullpath'));
addpath(root,fullfile(root,'support'),fullfile(root,'tests'));
core=run_paper_results(MakePlots=options.MakePlots);
literature=run_literature_pi_comparison(MakePlots=options.MakePlots);
checks=run_all_regressions();
assert(checks.Passed,'DCMG:FullComparisonFailed','Complete comparison regression failed.');
report=struct('Core',core,'LiteraturePI',literature,'Checks',checks, ...
    'Summary',literature.Summary);
save(fullfile(root,'Results','LiteraturePI','suite_checks.mat'),'checks');
fprintf('Complete five-controller comparison and all regression checks passed.\n');
end

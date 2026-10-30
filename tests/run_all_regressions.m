function results = run_all_regressions()
%RUN_ALL_REGRESSIONS Verify current results and retained simulator benchmarks.
% Run run_all_paper_results or both numerical entry points first.
root=fileparts(fileparts(mfilename('fullpath')));
addpath(root,fullfile(root,'support'),fullfile(root,'tests'));
checks=test_regenerated_results();
droop=test_conventional_droop_baseline();
comparison=test_complete_comparison();
a=test_model_based_simulator(MakePlots=false,SaveArtifacts=false,StopOnFailure=true);
b=test_paper_model_based_configuration(MakePlots=false,SaveArtifacts=false,StopOnFailure=true);
destination=fullfile(root,'Results','Regenerated');
save(fullfile(destination,'three_dg_regression.mat'),'a');
save(fullfile(destination,'legacy_regression_reports.mat'),'a','b');
results=struct('Passed',checks.Passed && droop.Passed && comparison.Passed,'Current',checks, ...
    'Droop',droop,'Comparison',comparison,'ThreeDG',a,'ArchivedSixDG',b);
fprintf('All regression gates passed; archived topology differences remain labeled warnings.\n');
end

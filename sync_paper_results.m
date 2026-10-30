function sync_paper_results(paperDirectory)
%SYNC_PAPER_RESULTS Copy generated publication assets to the manuscript repo.
%   sync_paper_results('/path/to/Data_Driven_Design_ACC2025')
% Run run_all_paper_results first. Rebuild exports from the verified records
% before copying, including after a numerical run with MakePlots=false.
arguments
    paperDirectory (1,1) string
end
root=fileparts(mfilename('fullpath'));
addpath(fullfile(root,'support'));
source=fullfile(root,'Results','Regenerated');
assert(isfile(fullfile(paperDirectory,'main.tex')),'DCMG:PaperDirectory', ...
    'The selected manuscript directory does not contain main.tex.');
s=load(fullfile(source,'paper_results.mat'),'report');
assert(s.report.Accepted && s.report.ImplementationChecks.Passed, ...
    'DCMG:UnverifiedResults','The regenerated results did not pass verification.');
literatureRecord=fullfile(root,'Results','LiteraturePI','literature_pi_results.mat');
assert(isfile(literatureRecord),'DCMG:MissingPIComparison', ...
    'Run run_all_paper_results before synchronizing the complete paper.');
p=load(literatureRecord,'report');
verifyCompleteComparison(s.report,p.report);
% Never publish stale plots or manually altered CSVs beside fresh MAT data.
writePaperResultsText(s.report,source);
plotPaperResults(s.report,source);
writetable(s.report.Summary,fullfile(source,'metrics.csv'));
writetable(s.report.DroopSensitivity.Table,fullfile(source,'droop_sensitivity.csv'));
destination=fullfile(paperDirectory,'Results','Regenerated');
if ~isfolder(destination), mkdir(destination); end
files=[dir(fullfile(source,'paper_*.pdf'));dir(fullfile(source,'paper_*.png'))];
names=[{files.name},{'numerical_results.tex','metrics.csv', ...
    'droop_sensitivity.csv','RESULTS_SUMMARY.md','artifact_sources.sha256'}];
for k=1:numel(names), copyfile(fullfile(source,names{k}),destination); end
fprintf('Copied verified numerical-results assets to %s\n',destination);

syncLiteraturePIComparison(paperDirectory);
end

function checks = test_complete_comparison()
%TEST_COMPLETE_COMPARISON Check shared provenance and exported numeric tables.
root=fileparts(fileparts(mfilename('fullpath')));
addpath(root,fullfile(root,'support'));
a=load(fullfile(root,'Results','Regenerated','paper_results.mat'),'report');
b=load(fullfile(root,'Results','LiteraturePI','literature_pi_results.mat'),'report');
checks=verifyCompleteComparison(a.report,b.report);
csv=readtable(fullfile(root,'Results','LiteraturePI','metrics.csv'),'TextType','string');
checks.CSVMatches=isequal(csv.Method,b.report.Summary.Method) && ...
    max(abs(csv{:,2:end}-b.report.Summary{:,2:end}),[],'all')<1e-10;
assert(checks.CSVMatches,'DCMG:ComparisonCSV','Exported metrics differ from the MAT record.');
% A newer core experiment must never be silently paired with an older PI run.
stale=b.report; stale.ReferenceGeneratedAt=stale.ReferenceGeneratedAt-seconds(1);
rejected=false;
try, verifyCompleteComparison(a.report,stale);
catch e, rejected=strcmp(e.identifier,'DCMG:StaleComparison'); end
assert(rejected,'DCMG:StaleComparisonGuard','A stale comparison was not rejected.');
checks.StaleRecordRejected=rejected;
fprintf('All five methods, source generations, trajectories, and CSV values are consistent.\n');
end

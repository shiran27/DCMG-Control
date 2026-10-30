function checks = verifyCompleteComparison(core,literature)
%VERIFYCOMPLETECOMPARISON Reject stale or mismatched comparison records.
assert(core.Accepted && core.ImplementationChecks.Passed && ...
    literature.Accepted && literature.Checks.Passed,'DCMG:UnverifiedComparison', ...
    'Both result records must pass their independent checks.');
checks.SameGeneration=isequal(core.GeneratedAt,literature.ReferenceGeneratedAt);
checks.SameSources=isequal(core.SourceHashes,literature.ReferenceSourceHashes);
checks.SamePlant=isequaln(core.PhysicalParameters,literature.PhysicalParameters);
checks.SameMethods=isequal(core.Methods,literature.Methods(1:numel(core.Methods)));
checks.MaxStateDifference=0; checks.MaxCommandDifference=0;
for field={'Nominal','Disturbance','BoundedNoise','ArchivedInitialState'}
    a=core.(field{1}); b=literature.(field{1});
    for k=1:numel(a)
        assert(isequal(a{k}.Time,b{k}.Time),'DCMG:ComparisonTimeGrid', ...
            'Comparison records use different time grids.');
        checks.MaxStateDifference=max(checks.MaxStateDifference, ...
            max(abs(a{k}.State-b{k}.State),[],'all'));
        checks.MaxCommandDifference=max(checks.MaxCommandDifference, ...
            max(abs(a{k}.Command-b{k}.Command),[],'all'));
    end
end
checks.Passed=checks.SameGeneration && checks.SameSources && checks.SamePlant ...
    && checks.SameMethods && checks.MaxStateDifference<1e-10 ...
    && checks.MaxCommandDifference<1e-10;
assert(checks.Passed,'DCMG:StaleComparison', ...
    'Result records are stale or inconsistent. Run run_all_paper_results.');
end

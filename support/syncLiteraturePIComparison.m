function syncLiteraturePIComparison(paperDirectory)
%SYNCLITERATUREPICOMPARISON Install the verified PI comparison into the paper.
% Adds one generated input after the existing numerical section, preserving
% other methods and paper edits. Repeated calls do not duplicate the input.
arguments
    paperDirectory (1,1) string
end
root=fileparts(fileparts(mfilename('fullpath')));
source=fullfile(root,'Results','LiteraturePI');
s=load(fullfile(source,'literature_pi_results.mat'),'report');
assert(s.report.Accepted && s.report.Checks.Passed,'DCMG:UnverifiedPI');
core=load(fullfile(root,'Results','Regenerated','paper_results.mat'),'report');
verifyCompleteComparison(core.report,s.report);
mainFile=fullfile(paperDirectory,'main.tex');
bibliography=fileread(fullfile(paperDirectory,'References.bib'));
assert(contains(bibliography,'{Nahata2020,'),'DCMG:MissingLiteratureCitation');
main=fileread(mainFile); original=main;
anchor='\input{Results/Regenerated/numerical_results}';
addition='\input{Results/LiteraturePI/comparison}';
assert(contains(main,anchor),'DCMG:PaperResultsAnchor');
% Recreate every export from the same verified record, even after no-plot runs.
writeLiteraturePIComparison(s.report,source);
plotLiteraturePIComparison(s.report,source);
writetable(s.report.Summary,fullfile(source,'metrics.csv'));
writetable(s.report.LiteratureTuning.Candidates,fullfile(source,'tuning.csv'));
destination=fullfile(paperDirectory,'Results','LiteraturePI');
if ~isfolder(destination), mkdir(destination); end
names={'literature_pi_comparison.pdf','literature_pi_comparison.png','comparison.tex', ...
    'metrics.csv','tuning.csv','RESULTS_SUMMARY.md'};
for k=1:numel(names), copyfile(fullfile(source,names{k}),destination); end
if ~contains(main,addition)
    main=strrep(main,anchor,[anchor,newline,addition]);
end
% Reconcile duplicate integration by another task without changing its prose.
while contains(main,[addition,newline,addition])
    main=strrep(main,[addition,newline,addition],addition);
end
if ~strcmp(main,original)
    fid=fopen(mainFile,'w'); assert(fid>=0); cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
    fprintf(fid,'%s',main);
end
fprintf('Installed the verified PI comparison in %s\n',paperDirectory);
end

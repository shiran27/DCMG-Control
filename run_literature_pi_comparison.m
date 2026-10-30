function report = run_literature_pi_comparison(options)
%RUN_LITERATURE_PI_COMPARISON Reproduce a published passivity-based PI baseline.
%   report=run_literature_pi_comparison;
% Uses verified proposed gains, then recomputes every comparison trajectory.
%   report=run_literature_pi_comparison(RebuildProposed=true);
% Re-synthesizes the proposed controllers first via run_paper_results.
arguments
    options.RebuildProposed (1,1) logical = false
    options.MakePlots (1,1) logical = true
    options.OutputDirectory (1,1) string = ""
end
root=fileparts(mfilename('fullpath')); addpath(root,fullfile(root,'support'));
if options.OutputDirectory=="", options.OutputDirectory=fullfile(root,'Results','LiteraturePI'); end
if ~isfolder(options.OutputDirectory), mkdir(options.OutputDirectory); end
if options.RebuildProposed
    base=run_paper_results(MakePlots=options.MakePlots);
else
    data=load(fullfile(root,'Results','Regenerated','paper_results.mat'),'report');
    base=data.report;
end
assert(base.Accepted && base.ImplementationChecks.Passed,'DCMG:UnverifiedReference');
[net,benchmark]=createPaperDCMicrogrid(LegacyEquilibrium=false,FinalTime=0.12);
assert(norm(net.x_s-base.PhysicalParameters.SteadyState.x_ss,inf)<1e-8, ...
    'DCMG:BaselinePlantMismatch');
for field={'L','C','Rf','RL','Ibar','Vrated','Irated'}
    actual=arrayfun(@(dg)dg.(field{1}),net.DGs);
    expected=base.PhysicalParameters.(field{1});
    assert(isequal(actual(:),expected(:)),'DCMG:BaselinePlantMismatch', ...
        'The reference report and current plant have different %s values.',field{1});
end
assert(norm(net.YBar-base.PhysicalParameters.Laplacian,inf)<1e-12, ...
    'DCMG:BaselinePlantMismatch','The physical line network has changed.');
initialState=base.Nominal{1}.State(1,:)';
peakBudget=max(base.Summary.PeakCommandDeviation(2:3));
[piController,tuning]=tuneLiteraturePI(net,initialState,peakBudget);
report=struct('GeneratedAt',datetime('now','TimeZone','UTC'), ...
    'MatlabVersion',version,'Options',options,'PhysicalParameters',base.PhysicalParameters, ...
    'ReferenceGeneratedAt',base.GeneratedAt,'ReferenceSourceHashes',base.SourceHashes, ...
    'Methods',base.Methods,'Controllers',{base.Controllers},'LiteratureTuning',tuning);
report.LiteratureIndex=numel(report.Controllers)+1; p=report.LiteratureIndex;
report.Controllers{p}=piController; report.Methods(p)="Passivity PI (Nahata et al.)";
metrics=zeros(p,11);
for k=1:p
    controller=report.Controllers{k};
    if k==p, gain=controller; else, gain=controller.CommandGain; end
    report.Nominal{k}=simulatePaperCase(net,gain,InitialState=initialState);
    report.Disturbance{k}=simulatePaperCase(net,gain,Scenario="disturbance");
    report.BoundedNoise{k}=simulatePaperCase(net,gain,Scenario="bounded-noise");
    report.LoadStep{k}=simulatePaperCase(net,gain,Scenario="load-step");
    report.ArchivedInitialState{k}=simulatePaperCase(net,gain, ...
        InitialState=benchmark.InitialState,FinalTime=0.05);
    a=report.Nominal{k}.Metrics; b=report.Disturbance{k}.Metrics;
    d=report.LoadStep{k}.Metrics;
    assert(report.Nominal{k}.Accepted && a.FinalRatio<1e-5,'DCMG:BaselineTrajectory');
    metrics(k,:)=[a.SettlingTime2Percent*1e3,b.VoltageRMS,b.SharingRMS, ...
        b.CommandRMSDeviation,a.PeakCommandDeviation,d.TailVoltageRMS,d.TailSharingRMS, ...
        a.MaxRealPole,controller.DirectedEdges,report.BoundedNoise{k}.Metrics.EmpiricalGain, ...
        a.FinalRatio];
end
report.Summary=[table(report.Methods','VariableNames',{'Method'}), ...
    array2table(metrics,'VariableNames',{'SettlingMs','DisturbanceVoltageRMS', ...
    'DisturbanceSharingRMS','DisturbanceCommandRMS','NominalPeakCommand', ...
    'LoadStepTailVoltageRMS','LoadStepTailSharingRMS','MaxRealPole','DirectedLinks', ...
    'EmpiricalNoiseGain','NominalFinalRatio'})];
report.Checks=verifyLiteraturePI(report);
assert(report.Checks.Passed,'DCMG:LiteratureChecks');
report.Accepted=true;
report.SourceHashes=hashSources(root);
save(fullfile(options.OutputDirectory,'literature_pi_results.mat'),'report','-v7.3');
writetable(report.Summary,fullfile(options.OutputDirectory,'metrics.csv'));
writetable(tuning.Candidates,fullfile(options.OutputDirectory,'tuning.csv'));
disp(report.Summary);
if options.MakePlots
    plotLiteraturePIComparison(report,options.OutputDirectory);
end
writeLiteraturePIComparison(report,options.OutputDirectory);
fprintf('Verified literature PI comparison saved in %s\n',options.OutputDirectory);
end

function hashes=hashSources(root)
names={'run_literature_pi_comparison.m','support/designLiteraturePI.m', ...
    'support/tuneLiteraturePI.m','support/simulatePaperCase.m', ...
    'support/verifyLiteraturePI.m','support/plotLiteraturePIComparison.m', ...
    'support/writeLiteraturePIComparison.m','DG.m','DCMicrogrid.m','TransmissionLine.m', ...
    'support/createPaperDCMicrogrid.m','support/designSteadyState.m'};
hashes=struct();
for k=1:numel(names)
    fid=fopen(fullfile(root,names{k}),'rb'); bytes=fread(fid,Inf,'*uint8'); fclose(fid);
    md=java.security.MessageDigest.getInstance('SHA-256'); md.update(bytes);
    hash=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
    hashes.(matlab.lang.makeValidName(names{k}))=hash;
end
end

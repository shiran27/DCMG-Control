function report = run_paper_results(options)
%RUN_PAPER_RESULTS Rebuild the paper benchmark, synthesize, verify, and plot.
%   report=run_paper_results;
%   report=run_paper_results(MakePlots=false);
% All artifacts are derived from this run, without loading archived gains.
% DD scope: exact CT integral-data local synthesis and known-B global design.
arguments
    options.MakePlots (1,1) logical = true
    options.OutputDirectory (1,1) string = ""
    options.RhoCandidates (1,:) double = [0.15,0.2,0.25,0.28]
    options.DroopFraction (1,1) double = 0.05
    options.DroopSensitivity (1,:) double = [0.01,0.03,0.05,0.10]
    options.RunVariation (1,1) logical = true
end
root=fileparts(mfilename('fullpath')); addpath(root,fullfile(root,'support'));
if options.OutputDirectory=="", options.OutputDirectory=fullfile(root,'Results','Regenerated'); end
if ~isfolder(options.OutputDirectory), mkdir(options.OutputDirectory); end
diaryFile=fullfile(options.OutputDirectory,'run.log');
if isfile(diaryFile), delete(diaryFile); end
diary(diaryFile); diaryCleanup=onCleanup(@()diary('off')); %#ok<NASGU>
rng(20260910,'twister'); ops=dcmgSolverOptions();
[net,benchmark]=createPaperDCMicrogrid(LegacyEquilibrium=false,FinalTime=0.12);
report=struct('GeneratedAt',datetime('now','TimeZone','UTC'),'MatlabVersion',version, ...
    'YalmipVersion',yalmip('version'),'Solver',ops.solver,'Options',options, ...
    'Benchmark',benchmark,'PhysicalParameters',physicalParameters(net), ...
    'Scope','Robust CT integral-data local design; known-B global composition', ...
    'Seed',20260910,'DesignOptions',paperDesignOptions());
report.SourceHashes=sourceHashes(root);
report.Data=collectIntegralData(net);
fprintf('Integral-data checks: %d DGs, %d transitions each.\n',net.N, ...
    report.Data.Options.Samples);
% Tune against certificates only. The same rho, budget, gain limits, graph
% threshold, and performance objective are used for both synthesis methods.
attempts=cell(numel(options.RhoCandidates),2);
scores=inf(numel(options.RhoCandidates),1);
for r=1:numel(options.RhoCandidates)
    net.DesignOptions.Rho=options.RhoCandidates(r);
    for d=1:2
        if d==1, [~,~,c]=net.codesign_MB_DRC(1e-4);
        else, [~,~,c]=net.codesign_DD_DRC(1e-4); end
        if c.problem==0
            c.AggregateVerification=verifyGlobalStorage(net,c);
            if ~c.AggregateVerification.Accepted
                c.problem=4; c.info='Original aggregate storage inequality failed';
            end
        end
        attempts{r,d}=c;
        fprintf('rho=%.3g, DD=%d, problem=%d',options.RhoCandidates(r),d-1,c.problem);
        if c.problem==0
            fprintf(', gamma=%.6g, directed links=%d\n',c.Gamma,c.DirectedEdges);
        else, fprintf(', %s\n',c.info); end
        save(fullfile(options.OutputDirectory,'design_attempts.mat'),'attempts');
    end
    if all(cellfun(@(c)c.problem==0,attempts(r,:)))
        scores(r)=max(attempts{r,1}.Gamma,attempts{r,2}.Gamma);
    end
end
[score,best]=min(scores);
assert(isfinite(score),'DCMG:NoAcceptedPair','No common MB/DD tuning passed all checks.');
report.Tuning=attempts; report.SelectedRho=options.RhoCandidates(best);
MB=attempts{best,1}; DD=attempts{best,2};
% Fresh stabilizing baseline from the existing global Lyapunov design and
% a solver-independent conventional V-I droop baseline from the literature.
[~,~,stabilizing]=net.design_MB_GSC(0,false);
assert(stabilizing.problem==0,'DCMG:BaselineFailure','Stabilizing baseline solve failed.');
KGSC=net.K;
stabilizing.CommandGain=KGSC;
stabilizing.DirectedEdges=countEdges(KGSC,net.N);
stabilizing.Gamma=NaN; stabilizing.Local={};
droop=designConventionalDroop(net,VoltageDroopFraction=options.DroopFraction);
report.Controllers={stabilizing,MB,DD,droop};
report.Methods=["MB stabilizing","MB co-design", ...
    "DD local + known-B global",string(droop.Name)];
gains={KGSC,MB.CommandGain,DD.CommandGain,droop.CommandGain};
numberOfMethods=numel(gains);
report.Nominal=cell(1,numberOfMethods);
report.Disturbance=cell(1,numberOfMethods);
report.BoundedNoise=cell(1,numberOfMethods);
report.ArchivedInitialState=cell(1,numberOfMethods);
perturb=repmat([1.5;0.35],net.N,1).*[1;-1;-1;1;1;1;-1;-1;1;-1;-1;1];
for k=1:numberOfMethods
    report.Nominal{k}=simulatePaperCase(net,gains{k}, ...
        InitialState=net.x_s+perturb,Certificate=report.Controllers{k});
    report.Disturbance{k}=simulatePaperCase(net,gains{k}, ...
        InitialState=net.x_s,Scenario="disturbance",Certificate=report.Controllers{k});
    report.BoundedNoise{k}=simulatePaperCase(net,gains{k}, ...
        InitialState=net.x_s,Scenario="bounded-noise",Certificate=report.Controllers{k});
    report.ArchivedInitialState{k}=simulatePaperCase(net,gains{k}, ...
        InitialState=benchmark.InitialState,FinalTime=0.05,Certificate=report.Controllers{k});
    assert(report.Nominal{k}.Accepted && report.Nominal{k}.Metrics.FinalRatio<1e-5, ...
        'DCMG:TrajectoryFailure','Nominal convergence check failed.');
end
report.DroopSensitivity=droopSensitivityStudy( ...
    net,options.DroopSensitivity,perturb);
report.Summary=summaryTable(report);
disp(report.Summary);
save(fullfile(options.OutputDirectory,'nominal_checkpoint.mat'),'report','-v7.3');
if options.RunVariation
    report.Variation=variationStudy(net,report,gains);
end
report.ImplementationChecks=verifyPaperImplementation(net,report);
assert(report.ImplementationChecks.Passed,'DCMG:ImplementationChecks', ...
    'One or more independent implementation checks failed.');
report.Accepted=true;
save(fullfile(options.OutputDirectory,'paper_results.mat'),'report','-v7.3');
writetable(report.Summary,fullfile(options.OutputDirectory,'metrics.csv'));
writetable(report.DroopSensitivity.Table, ...
    fullfile(options.OutputDirectory,'droop_sensitivity.csv'));
writePaperResultsText(report,options.OutputDirectory);
if options.MakePlots, plotPaperResults(report,options.OutputDirectory); end
fprintf('Accepted results saved in %s\n',options.OutputDirectory);
end

function s=summaryTable(r)
n=numel(r.Methods); links=zeros(n,1); gamma=nan(n,1); pole=zeros(n,1);
settle=zeros(n,1); final=zeros(n,1); voltage=zeros(n,1); sharing=zeros(n,1);
effort=zeros(n,1); noiseGain=zeros(n,1); localMargin=nan(n,1); globalMargin=nan(n,1);
for k=1:n
    c=r.Controllers{k}; m=r.Nominal{k}.Metrics;
    links(k)=c.DirectedEdges; gamma(k)=c.Gamma; pole(k)=m.MaxRealPole;
    settle(k)=1000*m.SettlingTime2Percent; final(k)=m.FinalRatio;
    voltage(k)=r.Disturbance{k}.Metrics.VoltageRMS;
    sharing(k)=r.Disturbance{k}.Metrics.SharingRMS;
    effort(k)=m.PeakCommandDeviation;
    noiseGain(k)=r.BoundedNoise{k}.Metrics.EmpiricalGain;
    if isfield(c,'Local') && ~isempty(c.Local) && ...
            isfield(c.Local{1},'TruthVerification')
        localMargin(k)=min(cellfun(@(l)l.TruthVerification.MinEigenvalue,c.Local));
    end
    if isfield(c,'AggregateVerification')
        globalMargin(k)=c.AggregateVerification.MinEigenvalue;
    end
end
s=table(r.Methods',links,gamma,pole,settle,final,voltage,sharing,effort,noiseGain, ...
    localMargin,globalMargin,'VariableNames',{'Method','DirectedLinks','Gamma', ...
    'MaxRealPole','SettlingMs','FinalErrorRatio','DisturbanceVoltageRMS', ...
    'DisturbanceSharingRMS','PeakCommandDeviation','EmpiricalNoiseGain', ...
    'OriginalLocalMinEig','OriginalGlobalMinEig'});
end

function sensitivity=droopSensitivityStudy(net,fractions,perturbation)
% The grid is declared in the run options and is not optimized on outcomes.
n=numel(fractions);
pole=zeros(n,1); settling=zeros(n,1); voltage=zeros(n,1);
sharing=zeros(n,1); effort=zeros(n,1); accepted=false(n,1);
controllers=cell(n,1);
for k=1:n
    c=designConventionalDroop(net,VoltageDroopFraction=fractions(k));
    nominal=simulatePaperCase(net,c.CommandGain, ...
        InitialState=net.x_s+perturbation,Certificate=c);
    disturbed=simulatePaperCase(net,c.CommandGain, ...
        InitialState=net.x_s,Scenario="disturbance",Certificate=c);
    controllers{k}=c;
    pole(k)=nominal.Metrics.MaxRealPole;
    settling(k)=1000*nominal.Metrics.SettlingTime2Percent;
    voltage(k)=disturbed.Metrics.VoltageRMS;
    sharing(k)=disturbed.Metrics.SharingRMS;
    effort(k)=nominal.Metrics.PeakCommandDeviation;
    accepted(k)=nominal.Accepted && disturbed.Accepted;
end
sensitivity=struct('DeclaredFractions',fractions,'Controllers',{controllers}, ...
    'Table',table(100*fractions(:),pole,settling,voltage,sharing,effort,accepted, ...
    'VariableNames',{'DroopPercent','MaxRealPole','SettlingMs', ...
    'DisturbanceVoltageRMS','DisturbanceSharingRMS', ...
    'PeakCommandDeviation','Accepted'}));
assert(all(accepted),'DCMG:DroopSensitivity', ...
    'At least one declared droop sensitivity case was unstable or nonfinite.');
end

function count=countEdges(K,N)
count=0;
for i=1:N
    for j=1:N
        count=count+(i~=j && norm(K(i,2*j-1:2*j))>1e-9);
    end
end
end

function p=physicalParameters(net)
p=struct('L',[net.DGs.L],'C',[net.DGs.C],'Rf',[net.DGs.Rf], ...
    'RL',[net.DGs.RL],'Ibar',[net.DGs.Ibar],'Vrated',[net.DGs.Vrated], ...
    'Irated',[net.DGs.Irated],'Positions',vertcat(net.DGs.pos), ...
    'LineI',[net.Lines.i],'LineJ',[net.Lines.j],'LineR',[net.Lines.R], ...
    'Laplacian',net.YBar,'SteadyState',net.ss);
end

function h=sourceHashes(root)
files=[dir(fullfile(root,'*.m'));dir(fullfile(root,'support','*.m'))];
h=struct();
for k=1:numel(files)
    fid=fopen(fullfile(files(k).folder,files(k).name),'rb');
    bytes=fread(fid,Inf,'*uint8'); fclose(fid);
    md=java.security.MessageDigest.getInstance('SHA-256'); md.update(bytes);
    hash=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
    h.(matlab.lang.makeValidName(files(k).name))=hash;
end
end

function v=variationStudy(net,report,gains)
% Declared physical variation: fixed gains are an empirical stress test.
% Redesign data are collected in a separate offline commissioning experiment.
[changed,~]=createPaperDCMicrogrid(LegacyEquilibrium=false,FinalTime=0.12);
rng(1618,'twister'); factors=0.9+0.2*rand(changed.N,5);
for i=1:changed.N
    dg=changed.DGs(i); fields={'L','C','Rf','RL','Ibar'};
    for j=1:numel(fields), dg.(fields{j})=dg.(fields{j})*factors(i,j); end
    dg.updateModel();
end
changed.buildSystemMatrices(); changed.solveSteadyState();
v=struct('Factors',factors,'Parameters',physicalParameters(changed), ...
    'Scope','Empirical physical-parameter stress test; no uncertainty certificate', ...
    'Frozen',{{}},'Redesigned',{{}});
x0=changed.x_s+(report.Nominal{1}.State(1,:)'-net.x_s);
for k=1:numel(gains)
    v.Frozen{k}=simulatePaperCase(changed,gains{k},InitialState=x0,Scenario="disturbance");
end
v.Data=collectIntegralData(changed,Seed=2719,InputMatrixNormBound=1.4);
v.RedesignAttempts={}; v.SelectedRho=NaN;
for rho=unique([report.SelectedRho,0.15,0.1])
    % Both methods use the same grid on the changed plant.
    changed.DesignOptions.Rho=rho; pair=cell(1,2);
    for d=1:2
        if d==1, [~,~,c]=changed.codesign_MB_DRC(1e-4);
        else, [~,~,c]=changed.codesign_DD_DRC(1e-4); end
        if c.problem==0
            c.AggregateVerification=verifyGlobalStorage(changed,c);
            if ~c.AggregateVerification.Accepted, c.problem=4; end
        end
        pair{d}=c;
    end
    v.RedesignAttempts{end+1}=pair;
    if all(cellfun(@(c)c.problem==0,pair))
        v.SelectedRho=rho;
        for d=1:2
            v.Redesigned{d}=simulatePaperCase(changed,pair{d}.CommandGain,InitialState=x0, ...
                Scenario="disturbance",Certificate=pair{d});
        end
    end
end
assert(isfinite(v.SelectedRho),'DCMG:VariationRedesign','No common variation tuning passed.');
end

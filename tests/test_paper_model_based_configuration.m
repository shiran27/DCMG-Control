function report = test_paper_model_based_configuration(options)
%TEST_PAPER_MODEL_BASED_CONFIGURATION Exercise the recovered paper benchmark.
%
%   report = test_paper_model_based_configuration()
%   report = test_paper_model_based_configuration(MakePlots=false)
%
% The test reconstructs the exact six-DG, seven-line physical configuration
% used by the current paper figures. Three controller cases are compared:
%   1. Current model-based global stabilizing design.
%   2. Archived model-based dissipative gains used for the paper panel.
%   3. Current model-based dissipative local/global co-design.
%
% Simulation uses the unbounded continuous-time linear error dynamics used
% by the original paper script. The staged data-driven workflow and active
% nonlinear clipping in DG.dynamics are intentionally excluded.

arguments
    options.MakePlots (1,1) logical = true
    options.SaveArtifacts (1,1) logical = true
    options.StopOnFailure (1,1) logical = false
    options.FinalTime (1,1) double {mustBePositive,mustBeFinite} = 0.05
    options.OutputStep (1,1) double {mustBePositive,mustBeFinite} = 1e-4
    options.LinkThreshold (1,1) double {mustBeNonnegative,mustBeFinite} = 1e-4
end

close all;
clc;

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(repoRoot,fullfile(repoRoot,'support'));

methodKeys = ["current-stabilizing","archived-dissipative", ...
    "current-dissipative"];
methodLabels = ["Current MB stabilizing","Archived paper MB dissipative", ...
    "Current MB dissipative"];
caseResults = repmat(emptyCaseResult(),numel(methodKeys),1);
caseNetworks = cell(numel(methodKeys),1);

fprintf('\nP6 MODEL-BASED PAPER-CONFIGURATION TEST\n');
fprintf('Physical benchmark: 6 DGs, 7 lines, final time %.4g s\n', ...
    options.FinalTime);
fprintf('Archived source: Git revision 6a0091d and archive/data/matlab.mat\n\n');

for methodIndex = 1:numel(methodKeys)
    methodKey = methodKeys(methodIndex);
    methodLabel = methodLabels(methodIndex);
    result = emptyCaseResult();
    result.Method = methodLabel;

    fprintf('--- %s ---\n',methodLabel);
    try
        [net,benchmark] = createPaperDCMicrogrid( ...
            FinalTime=options.FinalTime,NoiseStep=options.OutputStep);
        validateBenchmark(net,benchmark);

        [design,Kcommand] = configureController(net,benchmark,methodKey, ...
            options.LinkThreshold);
        result.DesignProblem = design.Problem;
        result.Message = design.Message;

        [t,X,metrics] = simulateLinearPaperCase(net,benchmark,Kcommand, ...
            options.FinalTime,options.OutputStep);

        result.ReferenceTopologyMatch = referenceTopologyMatch( ...
            net,benchmark,methodKey);
        result.Status = classifyResult( ...
            design.Problem,metrics,result.ReferenceTopologyMatch);
        result.MaxRealPole = metrics.MaxRealPole;
        result.InitialError = metrics.ErrorNorm(1);
        result.FinalError = metrics.ErrorNorm(end);
        result.FinalRatio = metrics.FinalRatio;
        result.EquilibriumResidualInf = metrics.EquilibriumResidualInf;
        result.OffDiagonalLinks = nnz(net.commAdj-diag(diag(net.commAdj)));
        result.ReferenceGainDistance = referenceGainDistance( ...
            net,benchmark,methodKey);
        result.Time = t;
        result.State = X;
        result.VoltagePU = metrics.VoltagePU;
        result.CurrentPU = metrics.CurrentPU;
        result.ErrorNorm = metrics.ErrorNorm;
        result.Poles = metrics.Poles;
        result.GlobalGain = net.K;
        result.LocalGain = collectLocalGain(net);
        result.CommunicationAdjacency = net.commAdj;

        net.setStateVector(X(end,:)');
        caseNetworks{methodIndex} = net;

        fprintf('status=%s, max Re(pole)=%.3e, final/error ratio=%.3e\n', ...
            result.Status,result.MaxRealPole,result.FinalRatio);
        fprintf('off-diagonal directed links=%d, equilibrium residual=%.3e\n\n', ...
            result.OffDiagonalLinks,result.EquilibriumResidualInf);
    catch exception
        result.Status = "ERROR";
        result.Message = string(getReport(exception,'extended','hyperlinks','off'));
        fprintf(2,'ERROR: %s\n\n',exception.message);
    end

    caseResults(methodIndex) = result;
end

summary = buildSummary(caseResults);
disp(summary);

report = struct();
report.GeneratedAt = datetime('now','TimeZone','local');
report.MatlabVersion = string(version);
report.Options = options;
report.BenchmarkSourceRevision = "6a0091d";
report.BenchmarkWorkspace = "archive/data/matlab.mat";
report.Summary = summary;
report.Cases = caseResults;
report.Artifacts = struct('MatFile',"",'ComparisonPngFile',"", ...
    'ComparisonFigFile',"",'NetworkPngFile',"",'NetworkFigFile',"");

if options.MakePlots
    comparisonFigure = plotPaperTrajectories(caseResults);
    networkFigure = plotPaperNetworkStates(caseResults,caseNetworks, ...
        options.FinalTime,options.OutputStep);
else
    comparisonFigure = gobjects(0);
    networkFigure = gobjects(0);
end

if options.SaveArtifacts
    resultsDirectory = fullfile(repoRoot,'Results');
    if ~isfolder(resultsDirectory)
        mkdir(resultsDirectory);
    end

    matFile = fullfile(resultsDirectory, ...
        'PaperConfigurationModelBasedComparison.mat');
    report.Artifacts.MatFile = string(matFile);

    if options.MakePlots
        comparisonPng = fullfile(resultsDirectory, ...
            'PaperConfigurationModelBasedComparison.png');
        comparisonFig = fullfile(resultsDirectory, ...
            'PaperConfigurationModelBasedComparison.fig');
        exportgraphics(comparisonFigure,comparisonPng,'Resolution',200, ...
            'BackgroundColor','white');
        savefig(comparisonFigure,comparisonFig);
        report.Artifacts.ComparisonPngFile = string(comparisonPng);
        report.Artifacts.ComparisonFigFile = string(comparisonFig);

        networkPng = fullfile(resultsDirectory, ...
            'PaperConfigurationNetworkStates.png');
        networkFig = fullfile(resultsDirectory, ...
            'PaperConfigurationNetworkStates.fig');
        exportgraphics(networkFigure,networkPng,'Resolution',200, ...
            'BackgroundColor','white');
        savefig(networkFigure,networkFig);
        report.Artifacts.NetworkPngFile = string(networkPng);
        report.Artifacts.NetworkFigFile = string(networkFig);
    end

    save(matFile,'report');
    fprintf('Saved paper-configuration report: %s\n',matFile);
end

failed = summary.Status == "FAIL" | summary.Status == "ERROR";
if any(failed)
    failedMethods = strjoin(summary.Method(failed),', ');
    warning('DCMG:PaperConfigurationIncomplete', ...
        'One or more paper-configuration cases failed: %s',failedMethods);
    if options.StopOnFailure
        error('DCMG:PaperConfigurationFailed', ...
            'The paper-configuration model-based test failed.');
    end
else
    fprintf('All paper-configuration cases completed without failure.\n');
end
end


function validateBenchmark(net,benchmark)
assert(net.N==6 && net.M==7,'DCMG:PaperBenchmarkSize', ...
    'The recovered paper benchmark must contain six DGs and seven lines.');
assert(isequal(vertcat(net.DGs.pos),benchmark.Positions), ...
    'DCMG:PaperBenchmarkPositions','Paper DG positions changed.');
actualEndpoints = [[net.Lines.i]' [net.Lines.j]'];
actualResistance = [net.Lines.R];
assert(isequal(actualEndpoints,benchmark.LineEndpoints), ...
    'DCMG:PaperBenchmarkTopology','Paper physical topology changed.');
assert(max(abs(actualResistance-benchmark.LineResistance))<1e-14, ...
    'DCMG:PaperBenchmarkResistance','Paper line resistances changed.');
assert(norm(net.getStateVector()-benchmark.InitialState,inf)<1e-12, ...
    'DCMG:PaperBenchmarkInitialState','Paper initial state changed.');
end


function [design,Kcommand] = configureController(net,benchmark,methodKey,threshold)
design = struct('Problem',0,'Message',"");
net.K = zeros(net.N,2*net.N);
for i = 1:net.N
    net.DGs(i).K = [0,0];
end

switch methodKey
    case "current-stabilizing"
        [~,~,out] = net.design_MB_GSC(threshold,false);
        design.Problem = out.problem;
        design.Message = string(out.info);
        assert(out.problem==0,'DCMG:PaperMBGSCFailure', ...
            'Current model-based stabilizing design failed: %s',out.info);

    case "archived-dissipative"
        net.K = benchmark.ReferenceGlobalGain;
        for i = 1:net.N
            net.DGs(i).K = benchmark.ReferenceLocalGain(i,:);
        end
        [adjacency,~] = net.buildCommAdjFromK(threshold);
        assert(isequal(adjacency,benchmark.ReferenceCommunicationAdjacency), ...
            'DCMG:PaperReferenceTopology', ...
            'Archived paper gains no longer recover the recorded topology.');
        design.Message = "Archived accepted paper gains";

    case "current-dissipative"
        [~,~,out] = net.codesign_MB_DRC(threshold);
        design.Problem = out.problem;
        design.Message = string(out.info)+ ...
            "; current implementation uses the active softened global LMI";
        assert(out.problem==0,'DCMG:PaperMBDRCFailure', ...
            'Current model-based dissipative design failed: %s',out.info);

    otherwise
        error('DCMG:UnknownPaperMethod','Unknown controller method %s.',methodKey);
end

Kcommand = net.K+expandLocalGain(net);
assert(isequal(size(Kcommand),[net.N,2*net.N]) && ...
    all(isfinite(Kcommand(:))),'DCMG:InvalidPaperController', ...
    'Recovered paper-configuration controller is invalid.');
end


function [t,X,metrics] = simulateLinearPaperCase( ...
        net,benchmark,Kcommand,finalTime,outputStep)
Aphysical = net.A+net.E*net.YBar*net.D';
Aclosed = Aphysical+net.B*Kcommand;
poles = eig(Aclosed);

equilibriumResidual = net.A*net.x_s+ ...
    net.E*(net.YBar*net.D'*net.x_s)+net.B*net.u_s+net.E*net.wBar;

timeGrid = 0:outputStep:finalTime;
if timeGrid(end)<finalTime
    timeGrid(end+1) = finalTime;
end
initialError = benchmark.InitialState-net.x_s;
odeOptions = odeset('RelTol',1e-7,'AbsTol',1e-9);
[t,errorTrajectory] = ode45(@(~,errorState)Aclosed*errorState, ...
    timeGrid,initialError,odeOptions);
X = errorTrajectory+net.x_s';

stateScale = zeros(2*net.N,1);
for i = 1:net.N
    stateScale(2*i-1:2*i) = [net.DGs(i).Vrated;net.DGs(i).Irated];
end
normalizedError = errorTrajectory./stateScale';
errorNorm = vecnorm(normalizedError,2,2);

metrics = struct();
metrics.Poles = poles;
metrics.MaxRealPole = max(real(poles));
metrics.EquilibriumResidualInf = norm(equilibriumResidual,inf);
metrics.ErrorNorm = errorNorm;
metrics.FinalRatio = errorNorm(end)/max(errorNorm(1),eps);
metrics.VoltagePU = X(:,1:2:end)./net.x_s(1:2:end)';
metrics.CurrentPU = X(:,2:2:end)./net.x_s(2:2:end)';
metrics.IsFinite = all(isfinite(X(:))) && all(isfinite(errorNorm));
end


function Klocal = expandLocalGain(net)
Klocal = zeros(net.N,2*net.N);
for i = 1:net.N
    Klocal(i,2*i-1:2*i) = net.DGs(i).K;
end
end


function Klocal = collectLocalGain(net)
Klocal = zeros(net.N,2);
for i = 1:net.N
    Klocal(i,:) = net.DGs(i).K;
end
end


function distance = referenceGainDistance(net,benchmark,methodKey)
if methodKey == "current-dissipative" || methodKey == "archived-dissipative"
    globalDistance = norm(net.K-benchmark.ReferenceGlobalGain,'fro')/ ...
        max(norm(benchmark.ReferenceGlobalGain,'fro'),eps);
    localGain = collectLocalGain(net);
    localDistance = norm(localGain-benchmark.ReferenceLocalGain,'fro')/ ...
        max(norm(benchmark.ReferenceLocalGain,'fro'),eps);
    distance = hypot(globalDistance,localDistance);
else
    distance = NaN;
end
end


function matches = referenceTopologyMatch(net,benchmark,methodKey)
if methodKey == "current-stabilizing"
    matches = isequal(net.commAdj,ones(net.N));
else
    matches = isequal(net.commAdj,benchmark.ReferenceCommunicationAdjacency);
end
end


function status = classifyResult(designProblem,metrics,referenceTopologyMatches)
if designProblem~=0 || ~metrics.IsFinite || metrics.EquilibriumResidualInf>1e-6
    status = "FAIL";
elseif metrics.MaxRealPole>=0 || metrics.FinalRatio>=0.10 || ...
        ~referenceTopologyMatches
    status = "WARN";
else
    status = "PASS";
end
end


function summary = buildSummary(caseResults)
summary = table( ...
    string({caseResults.Method})',string({caseResults.Status})', ...
    [caseResults.DesignProblem]',[caseResults.EquilibriumResidualInf]', ...
    [caseResults.MaxRealPole]',[caseResults.FinalRatio]', ...
    [caseResults.OffDiagonalLinks]',[caseResults.ReferenceTopologyMatch]', ...
    [caseResults.ReferenceGainDistance]', ...
    'VariableNames',{'Method','Status','DesignProblem', ...
    'EquilibriumResidualInf','MaxRealPole','FinalRatio', ...
    'OffDiagonalDirectedLinks','ReferenceTopologyMatch', ...
    'ReferenceGainDistance'});
end


function fig = plotPaperTrajectories(caseResults)
valid = find(arrayfun(@(item)~isempty(item.Time),caseResults));
fig = figure('Color','w','Name','P6 Paper Configuration Comparison', ...
    'Position',[80,80,1250,650]);
layout = tiledlayout(fig,2,numel(valid),'TileSpacing','compact', ...
    'Padding','compact');
colors = lines(6);

for column = 1:numel(valid)
    item = caseResults(valid(column));
    voltageAxes = nexttile(layout,column); hold(voltageAxes,'on');
    currentAxes = nexttile(layout,numel(valid)+column); hold(currentAxes,'on');
    styleDiagnosticAxes(voltageAxes);
    styleDiagnosticAxes(currentAxes);
    for dg = 1:6
        plot(voltageAxes,item.Time,item.VoltagePU(:,dg),'LineWidth',1.1, ...
            'Color',colors(dg,:),'DisplayName',sprintf('DG %d',dg));
        plot(currentAxes,item.Time,item.CurrentPU(:,dg),'LineWidth',1.1, ...
            'Color',colors(dg,:),'DisplayName',sprintf('DG %d',dg));
    end
    yline(voltageAxes,1,'k--','HandleVisibility','off');
    yline(currentAxes,1,'k--','HandleVisibility','off');
    grid(voltageAxes,'on'); grid(currentAxes,'on');
    title(voltageAxes,sprintf('%s (%d directed links)', ...
        item.Method,item.OffDiagonalLinks),'Interpreter','none', ...
        'Color',[0.1,0.1,0.1]);
    ylabel(voltageAxes,'V_i [pu]');
    ylabel(currentAxes,'I_{ti} [pu]');
    xlabel(currentAxes,'Time [s]');
    ylim(voltageAxes,[-2,2]);
    ylim(currentAxes,[-6,4]);
    if column==1
        plotLegend = legend(voltageAxes,'Location','southoutside', ...
            'NumColumns',3);
        set(plotLegend,'Color','w','TextColor',[0.15,0.15,0.15], ...
            'EdgeColor',[0.75,0.75,0.75]);
    end
end
title(layout,'Recovered P6 Six-DG Model-Based Benchmark', ...
    'Color',[0.1,0.1,0.1]);
end


function fig = plotPaperNetworkStates( ...
        caseResults,caseNetworks,finalTime,outputStep)
valid = find(arrayfun(@(item)~isempty(item.Time),caseResults));
fig = figure('Color','w','Name','P6 Paper Configuration Network States', ...
    'Position',[60,60,1250,900]);
layout = tiledlayout(fig,2,2,'TileSpacing','compact', ...
    'Padding','compact');

[initialNetwork,~] = createPaperDCMicrogrid( ...
    FinalTime=finalTime,NoiseStep=outputStep);
initialAxes = nexttile(layout,1);
initialNetwork.draw(initialAxes,'Title','Initial paper configuration', ...
    'ShowCommunication',false,'ShowLineCurrents',false, ...
    'ShowLineLabels',false,'FontSize',7);
set(initialAxes,'Color','w');

for column = 1:numel(valid)
    itemIndex = valid(column);
    item = caseResults(itemIndex);
    net = caseNetworks{itemIndex};
    net.setStateVector(item.State(end,:)');
    ax = nexttile(layout,column+1);
    net.draw(ax,'Title',char(item.Method),'ShowCommunication',true, ...
        'ShowLineCurrents',false,'ShowLineLabels',false,'FontSize',7);
    set(ax,'Color','w');
end
end


function styleDiagnosticAxes(ax)
set(ax,'Color','w','XColor',[0.2,0.2,0.2],'YColor',[0.2,0.2,0.2], ...
    'GridColor',[0.45,0.45,0.45],'GridAlpha',0.22, ...
    'MinorGridAlpha',0.12);
end


function result = emptyCaseResult()
result = struct( ...
    'Method',"",'Status',"NOT RUN",'DesignProblem',NaN,'Message',"", ...
    'MaxRealPole',NaN,'InitialError',NaN,'FinalError',NaN, ...
    'FinalRatio',NaN,'EquilibriumResidualInf',NaN, ...
    'OffDiagonalLinks',NaN,'ReferenceTopologyMatch',false, ...
    'ReferenceGainDistance',NaN, ...
    'Time',[],'State',[],'VoltagePU',[],'CurrentPU',[], ...
    'ErrorNorm',[],'Poles',[],'GlobalGain',[],'LocalGain',[], ...
    'CommunicationAdjacency',[]);
end

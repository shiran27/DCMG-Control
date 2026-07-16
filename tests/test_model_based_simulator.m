function report = test_model_based_simulator(options)
%TEST_MODEL_BASED_SIMULATOR Exercise the current simulator with model-based controls.
%
%   report = test_model_based_simulator()
%   report = test_model_based_simulator(MakePlots=false,SaveArtifacts=false)
%
% The test reconstructs the same deterministic three-DG microgrid for each
% case so that handle-object state and controller mutations cannot leak
% between methods. It exercises:
%   1. Open loop (steady command only)
%   2. Existing model-based global stabilizing design
%   3. Centralized continuous-time LQR
%   4. Existing model-based dissipative local/global co-design
%
% This is a baseline diagnostic, not a replacement for the paper workflow.
% It does not write netFile.mat or Results/TrajResults.mat.

arguments
    options.MakePlots (1,1) logical = true
    options.SaveArtifacts (1,1) logical = true
    options.StopOnFailure (1,1) logical = false
    options.FinalTime (1,1) double {mustBePositive,mustBeFinite} = 0.08
    options.OutputStep (1,1) double {mustBePositive,mustBeFinite} = 1e-4
end

close all;
clc;
rng(17,'twister');

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(repoRoot);

methodKeys = ["open-loop","mb-stabilizing","lqr","mb-dissipative"];
methodLabels = ["Open loop","MB stabilizing","LQR","MB dissipative"];
cases = repmat(emptyCaseResult(),numel(methodKeys),1);
caseNetworks = cell(numel(methodKeys),1);

fprintf('\nMODEL-BASED DCMG SIMULATOR BASELINE\n');
fprintf('Repository: %s\n',repoRoot);
fprintf('Final time: %.4g s, output step: %.4g s\n\n', ...
    options.FinalTime,options.OutputStep);

for methodIndex = 1:numel(methodKeys)
    methodKey = methodKeys(methodIndex);
    methodLabel = methodLabels(methodIndex);
    fprintf('--- %s ---\n',methodLabel);

    result = emptyCaseResult();
    result.Method = methodLabel;

    try
        [net,x0,stateScale] = makeExampleDCMG(options.FinalTime,options.OutputStep);
        [design,Kcommand] = configureController(net,methodKey,stateScale);
        result.DesignProblem = design.Problem;

        Aphysical = net.A + net.E*net.YBar*net.D';
        Aclosed = Aphysical + net.B*Kcommand;
        poles = eig(Aclosed);

        net.setStateVector(net.x_s);
        equilibriumDerivative = net.dynamics(0,net.x_s,false);
        equilibriumResidual = norm(equilibriumDerivative,inf);

        timeGrid = makeTimeGrid(options.FinalTime,options.OutputStep);
        net.setStateVector(x0);
        odeOptions = odeset('RelTol',1e-7,'AbsTol',1e-9);
        [t,X] = ode45(@(time,state)net.dynamics(time,state,false), ...
            timeGrid,x0,odeOptions);

        metrics = evaluateTrajectory(net,t,X,Kcommand,stateScale);
        net.setStateVector(X(end,:)');

        result.Status = classifyResult(equilibriumResidual,poles,metrics);
        result.Message = design.Message;
        result.EquilibriumResidualInf = equilibriumResidual;
        result.MaxRealPole = max(real(poles));
        result.InitialError = metrics.ErrorNorm(1);
        result.FinalError = metrics.ErrorNorm(end);
        result.FinalRatio = result.FinalError/max(result.InitialError,eps);
        result.PeakError = max(metrics.ErrorNorm);
        result.MaxVoltageErrorV = max(metrics.VoltageMaxAbsError);
        result.MaxCurrentErrorA = max(metrics.CurrentMaxAbsError);
        result.MaxCommandDeviationV = max(metrics.CommandMaxAbsDeviation);
        result.CommandSaturationCount = metrics.CommandSaturationCount;
        result.CurrentClippingCount = metrics.CurrentClippingCount;
        result.Poles = poles;
        result.KCommand = Kcommand;
        result.Time = t;
        result.State = X;
        result.ErrorNorm = metrics.ErrorNorm;
        result.VoltageRmsError = metrics.VoltageRmsError;
        result.CommandRmsDeviation = metrics.CommandRmsDeviation;
        caseNetworks{methodIndex} = net;

        fprintf('status=%s, max Re(pole)=%.3e, final/error ratio=%.3e\n', ...
            result.Status,result.MaxRealPole,result.FinalRatio);
        fprintf('equilibrium residual=%.3e, command saturation=%d, line-current limiting=%d\n\n', ...
            result.EquilibriumResidualInf,result.CommandSaturationCount, ...
            result.CurrentClippingCount);
    catch exception
        result.Status = "ERROR";
        result.Message = string(getReport(exception,'extended','hyperlinks','off'));
        fprintf(2,'ERROR: %s\n\n',exception.message);
    end

    cases(methodIndex) = result;
end

summary = buildSummary(cases);
disp(summary);

report = struct();
report.GeneratedAt = datetime('now','TimeZone','local');
report.MatlabVersion = string(version);
report.Options = options;
report.Summary = summary;
report.Cases = cases;
report.Artifacts = struct( ...
    'MatFile',"", ...
    'PngFile',"", ...
    'FigFile',"", ...
    'NetworkPngFile',"", ...
    'NetworkFigFile',"");

if options.MakePlots
    comparisonFigure = plotComparison(cases);
    networkFigure = plotNetworkStates(cases,caseNetworks);
else
    comparisonFigure = gobjects(0);
    networkFigure = gobjects(0);
end

if options.SaveArtifacts
    resultsDirectory = fullfile(repoRoot,'Results');
    if ~isfolder(resultsDirectory)
        mkdir(resultsDirectory);
    end

    matFile = fullfile(resultsDirectory,'ModelBasedSimulatorBaseline.mat');
    report.Artifacts.MatFile = string(matFile);

    if options.MakePlots
        pngFile = fullfile(resultsDirectory,'ModelBasedSimulatorBaseline.png');
        figFile = fullfile(resultsDirectory,'ModelBasedSimulatorBaseline.fig');
        exportgraphics(comparisonFigure,pngFile,'Resolution',200, ...
            'BackgroundColor','white');
        savefig(comparisonFigure,figFile);
        report.Artifacts.PngFile = string(pngFile);
        report.Artifacts.FigFile = string(figFile);

        networkPngFile = fullfile(resultsDirectory, ...
            'ModelBasedSimulatorNetworkStates.png');
        networkFigFile = fullfile(resultsDirectory, ...
            'ModelBasedSimulatorNetworkStates.fig');
        exportgraphics(networkFigure,networkPngFile,'Resolution',200, ...
            'BackgroundColor','white');
        savefig(networkFigure,networkFigFile);
        report.Artifacts.NetworkPngFile = string(networkPngFile);
        report.Artifacts.NetworkFigFile = string(networkFigFile);
    end

    save(matFile,'report');
    fprintf('Saved baseline report: %s\n',matFile);
end

failed = summary.Status ~= "PASS";
if any(failed)
    failedMethods = strjoin(summary.Method(failed),', ');
    warning('DCMG:BaselineIncomplete', ...
        'One or more baseline cases did not pass: %s',failedMethods);
    if options.StopOnFailure
        error('DCMG:BaselineFailed','Model-based simulator baseline failed.');
    end
else
    fprintf('All model-based simulator baseline cases passed.\n');
end
end


function [net,x0,stateScale] = makeExampleDCMG(finalTime,outputStep)
% Deterministic, moderately heterogeneous 3-DG ring.

N = 3;
positions = [0,0;2.8,0;1.4,2.2];
L = 1e-3*[1.00,1.08,0.94];
C = 1e-3*[1.00,0.92,1.10];
Rf = [0.10,0.12,0.09];
RL = [15,18,16];
Ibar = [0.45,0.55,0.50];
Vrated = 48*ones(1,N);
Irated = 7*ones(1,N);
nodeColors = [ ...
    0.00,0.45,0.70; ...
    0.84,0.37,0.00; ...
    0.00,0.62,0.45];

DGs = DG.empty(0,N);
for i = 1:N
    parameters = struct( ...
        'L',L(i),'C',C(i),'Rf',Rf(i),'RL',RL(i),'Ibar',Ibar(i), ...
        'Vrated',Vrated(i),'Irated',Irated(i),'u_s',Vrated(i), ...
        'pos',positions(i,:),'color',nodeColors(i,:), ...
        'x',[Vrated(i);0],'K',[0,0]);
    DGs(i) = DG(i,parameters);
end

lineCells = { ...
    TransmissionLine(1,1,2,0.38), ...
    TransmissionLine(2,2,3,0.44), ...
    TransmissionLine(3,3,1,0.41)};
lines = [lineCells{:}];

net = DCMicrogrid(DGs,lines);
net.buildSystemMatrices();
steadyState = net.solveSteadyState();
assert(steadyState.problem==0,'DCMG:SteadyStateFailure', ...
    'The example microgrid steady-state problem is infeasible.');

net.K = zeros(N,2*N);
for i = 1:N
    net.DGs(i).K = [0,0];
end

noiseTime = 0:outputStep:finalTime;
net.setupNoise(noiseTime,outputStep,zeros(2));

% Keep the nominal verification inside the linear design region. Larger
% perturbations can be supplied in a separate nonlinear limiter stress test.
relativePerturbation = [0.025;-0.040;-0.020;0.035;0.0175;-0.030];
x0 = net.x_s.*(1+relativePerturbation);
net.setStateVector(x0);

stateScale = zeros(2*N,1);
for i = 1:N
    stateScale(2*i-1:2*i) = [net.DGs(i).Vrated;net.DGs(i).Irated];
end
end


function [design,Kcommand] = configureController(net,methodKey,stateScale)
N = net.N;
design = struct('Problem',0,'Message',"");

for i = 1:N
    net.DGs(i).K = [0,0];
end
net.K = zeros(N,2*N);

switch methodKey
    case "open-loop"
        design.Message = "Steady-state feedforward only";

    case "mb-stabilizing"
        [~,~,out] = net.design_MB_GSC(0,false);
        design.Problem = out.problem;
        design.Message = string(out.info);
        assert(out.problem==0,'DCMG:MBGSCFailure', ...
            'Existing model-based stabilizing design failed: %s',out.info);

    case "lqr"
        assert(exist('lqr','file')==2,'DCMG:LQRUnavailable', ...
            'The Control System Toolbox lqr function is unavailable.');
        Aphysical = net.A + net.E*net.YBar*net.D';
        Bphysical = net.B;

        qDiagonal = zeros(2*N,1);
        rDiagonal = zeros(N,1);
        for i = 1:N
            voltageScale = stateScale(2*i-1);
            currentScale = stateScale(2*i);
            qDiagonal(2*i-1:2*i) = [20/voltageScale^2;1/currentScale^2];
            rDiagonal(i) = 0.5/voltageScale^2;
        end
        [lqrGain,~,lqrPoles] = lqr(Aphysical,Bphysical, ...
            diag(qDiagonal),diag(rDiagonal));
        net.K = -lqrGain;
        design.Message = sprintf('Continuous-time LQR; max Re(pole)=%.3e', ...
            max(real(lqrPoles)));

    case "mb-dissipative"
        [~,~,out] = net.codesign_MB_DRC(0);
        design.Problem = out.problem;
        design.Message = string(out.info);
        assert(out.problem==0,'DCMG:MBDRCFailure', ...
            'Existing model-based dissipative co-design failed: %s',out.info);

    otherwise
        error('DCMG:UnknownMethod','Unknown controller method %s.',methodKey);
end

Kcommand = net.K + localCommandGain(net);
assert(isequal(size(Kcommand),[N,2*N]) && all(isfinite(Kcommand(:))), ...
    'DCMG:InvalidController','Recovered command gain is invalid.');
end


function Klocal = localCommandGain(net)
Klocal = zeros(net.N,2*net.N);
for i = 1:net.N
    Klocal(i,2*i-1:2*i) = net.DGs(i).K;
end
end


function metrics = evaluateTrajectory(net,t,X,Kcommand,stateScale)
stateError = X-net.x_s';
normalizedError = stateError./stateScale';
metrics.ErrorNorm = vecnorm(normalizedError,2,2);

sampleCount = numel(t);
voltageMaxAbsError = zeros(sampleCount,1);
currentMaxAbsError = zeros(sampleCount,1);
commandMaxAbsDeviation = zeros(sampleCount,1);
voltageRmsError = zeros(sampleCount,1);
commandRmsDeviation = zeros(sampleCount,1);
commandSaturationCount = 0;
currentClippingCount = 0;

steadyVoltage = net.D'*net.x_s;
for k = 1:sampleCount
    state = X(k,:)';
    error = state-net.x_s;
    voltage = net.D'*state;
    internalCurrent = net.DBar'*state;
    command = net.u_s+Kcommand*error;
    lineCurrent = net.YBar*voltage;

    voltageError = voltage-steadyVoltage;
    currentError = internalCurrent-net.DBar'*net.x_s;
    commandDeviation = command-net.u_s;

    voltageMaxAbsError(k) = max(abs(voltageError));
    currentMaxAbsError(k) = max(abs(currentError));
    commandMaxAbsDeviation(k) = max(abs(commandDeviation));
    voltageRmsError(k) = sqrt(mean(voltageError.^2));
    commandRmsDeviation(k) = sqrt(mean(commandDeviation.^2));

    for i = 1:net.N
        commandSaturationCount = commandSaturationCount + ...
            double(command(i)>2*net.DGs(i).Vrated || ...
            command(i)<-2*net.DGs(i).Vrated);
        currentClippingCount = currentClippingCount + ...
            double(abs(lineCurrent(i))>2*net.DGs(i).Irated);
    end
end

metrics.VoltageMaxAbsError = voltageMaxAbsError;
metrics.CurrentMaxAbsError = currentMaxAbsError;
metrics.CommandMaxAbsDeviation = commandMaxAbsDeviation;
metrics.VoltageRmsError = voltageRmsError;
metrics.CommandRmsDeviation = commandRmsDeviation;
metrics.CommandSaturationCount = commandSaturationCount;
metrics.CurrentClippingCount = currentClippingCount;
metrics.IsFinite = all(isfinite(X(:))) && all(isfinite(metrics.ErrorNorm));
end


function status = classifyResult(equilibriumResidual,poles,metrics)
finalRatio = metrics.ErrorNorm(end)/max(metrics.ErrorNorm(1),eps);
if ~metrics.IsFinite || equilibriumResidual>1e-6
    status = "FAIL";
elseif max(real(poles))>=1e-7 || finalRatio>=1 || ...
        metrics.CommandSaturationCount>0 || metrics.CurrentClippingCount>0
    status = "WARN";
else
    status = "PASS";
end
end


function summary = buildSummary(cases)
summary = table( ...
    string({cases.Method})', ...
    string({cases.Status})', ...
    [cases.DesignProblem]', ...
    [cases.EquilibriumResidualInf]', ...
    [cases.MaxRealPole]', ...
    [cases.FinalRatio]', ...
    [cases.MaxVoltageErrorV]', ...
    [cases.MaxCommandDeviationV]', ...
    [cases.CommandSaturationCount]', ...
    [cases.CurrentClippingCount]', ...
    'VariableNames',{'Method','Status','DesignProblem','EquilibriumResidualInf', ...
    'MaxRealPole','FinalRatio','MaxVoltageErrorV','MaxCommandDeviationV', ...
    'CommandSaturationCount','CurrentClippingCount'});
end


function comparisonFigure = plotComparison(cases)
valid = find(arrayfun(@(result)~isempty(result.Time),cases));
comparisonFigure = figure('Color','w','Name','Model-Based Simulator Baseline', ...
    'Position',[100,100,1000,700]);
layout = tiledlayout(comparisonFigure,2,2,'TileSpacing','compact', ...
    'Padding','compact');
colors = lines(max(numel(valid),1));

ax1 = nexttile(layout,1); hold(ax1,'on'); grid(ax1,'on');
ax2 = nexttile(layout,2); hold(ax2,'on'); grid(ax2,'on');
ax3 = nexttile(layout,3); hold(ax3,'on'); grid(ax3,'on');
for plotIndex = 1:numel(valid)
    result = cases(valid(plotIndex));
    semilogy(ax1,result.Time,max(result.ErrorNorm,1e-12), ...
        'LineWidth',1.4,'Color',colors(plotIndex,:), ...
        'DisplayName',result.Method);
    plot(ax2,result.Time,result.VoltageRmsError, ...
        'LineWidth',1.4,'Color',colors(plotIndex,:), ...
        'DisplayName',result.Method);
    plot(ax3,result.Time,result.CommandRmsDeviation, ...
        'LineWidth',1.4,'Color',colors(plotIndex,:), ...
        'DisplayName',result.Method);
end
xlabel(ax1,'Time [s]'); ylabel(ax1,'Normalized state-error norm');
xlabel(ax2,'Time [s]'); ylabel(ax2,'Voltage RMS error [V]');
xlabel(ax3,'Time [s]'); ylabel(ax3,'Command RMS deviation [V]');
legend(ax1,'Location','best');

ax4 = nexttile(layout,4); grid(ax4,'on');
methodNames = string({cases.Method});
maxRealPoles = [cases.MaxRealPole];
bar(ax4,categorical(methodNames,methodNames),maxRealPoles, ...
    'FaceColor',[0.12,0.48,0.55]);
yline(ax4,0,'k--','LineWidth',0.8);
ylabel(ax4,'Maximum real pole');
title(layout,'Deterministic 3-DG Model-Based Simulator Baseline', ...
    'Color','k');

allAxes = [ax1,ax2,ax3,ax4];
set(allAxes,'Color','w','XColor','k','YColor','k', ...
    'GridColor',[0.75,0.75,0.75],'GridAlpha',0.5);
for ax = allAxes
    ax.Toolbar.Visible = 'off';
end
legendHandle = legend(ax1);
set(legendHandle,'Color','w','TextColor','k','EdgeColor',[0.4,0.4,0.4]);
end


function networkFigure = plotNetworkStates(cases,caseNetworks)
valid = find(arrayfun(@(result)~isempty(result.Time),cases));
assert(~isempty(valid) && all(~cellfun(@isempty,caseNetworks(valid))), ...
    'DCMG:MissingNetworkSnapshot', ...
    'A network snapshot is required for every plotted case.');

networkFigure = figure('Color','w','Name','Model-Based Network States', ...
    'Position',[60,80,1450,820]);
layout = tiledlayout(networkFigure,2,3,'TileSpacing','compact', ...
    'Padding','compact');

initialIndex = valid(1);
initialNetwork = caseNetworks{initialIndex};
terminalState = initialNetwork.getStateVector();
initialNetwork.setStateVector(cases(initialIndex).State(1,:)');
initialAxes = nexttile(layout,1);
initialNetwork.draw(initialAxes, ...
    'Title','Initial network state', ...
    'ShowCommunication',false, ...
    'ShowLineCurrents',true, ...
    'FontSize',7.5);
initialNetwork.setStateVector(terminalState);

for plotIndex = 1:numel(valid)
    caseIndex = valid(plotIndex);
    network = caseNetworks{caseIndex};
    network.setStateVector(cases(caseIndex).State(end,:)');
    terminalAxes = nexttile(layout,plotIndex+1);
    terminalTitle = sprintf('%s terminal state (%.0f ms)', ...
        cases(caseIndex).Method,1e3*cases(caseIndex).Time(end));
    network.draw(terminalAxes, ...
        'Title',terminalTitle, ...
        'ShowCommunication',true, ...
        'ShowLineCurrents',true, ...
        'FontSize',7.5);
end

summaryAxes = nexttile(layout,6);
terminalErrors = [cases(valid).FinalError];
stem(summaryAxes,1:numel(valid),terminalErrors,'filled', ...
    'Color',[0.12,0.48,0.55],'LineWidth',1.4,'MarkerSize',5);
set(summaryAxes,'YScale','log','Color','w','XColor','k','YColor','k', ...
    'XTick',1:numel(valid),'XTickLabel',string({cases(valid).Method}), ...
    'XTickLabelRotation',20,'GridColor',[0.75,0.75,0.75], ...
    'GridAlpha',0.5);
grid(summaryAxes,'on');
xlim(summaryAxes,[0.5,numel(valid)+0.5]);
ylabel(summaryAxes,'Terminal normalized state error');
title(summaryAxes,'Terminal residual','Color',[0.08,0.18,0.28]);

allAxes = findall(networkFigure,'Type','axes');
for ax = allAxes'
    ax.Toolbar.Visible = 'off';
end
end


function timeGrid = makeTimeGrid(finalTime,outputStep)
timeGrid = 0:outputStep:finalTime;
if timeGrid(end)<finalTime
    timeGrid = [timeGrid,finalTime];
end
end


function result = emptyCaseResult()
result = struct( ...
    'Method',"", ...
    'Status',"NOT RUN", ...
    'Message',"", ...
    'DesignProblem',NaN, ...
    'EquilibriumResidualInf',NaN, ...
    'MaxRealPole',NaN, ...
    'InitialError',NaN, ...
    'FinalError',NaN, ...
    'FinalRatio',NaN, ...
    'PeakError',NaN, ...
    'MaxVoltageErrorV',NaN, ...
    'MaxCurrentErrorA',NaN, ...
    'MaxCommandDeviationV',NaN, ...
    'CommandSaturationCount',NaN, ...
    'CurrentClippingCount',NaN, ...
    'Poles',[], ...
    'KCommand',[], ...
    'Time',[], ...
    'State',[], ...
    'ErrorNorm',[], ...
    'VoltageRmsError',[], ...
    'CommandRmsDeviation',[]);
end

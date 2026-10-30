function checks = verifyLiteraturePI(report)
%VERIFYLITERATUREPI Cross-check the published absolute law and existing ODE.
[plant,~]=createPaperDCMicrogrid(LegacyEquilibrium=false);
index=report.LiteratureIndex; c=report.Controllers{index}; n=2*plant.N;
e0=report.Nominal{index}.State(1,:)'-plant.x_s;
y0=[plant.x_s+e0;zeros(plant.N,1)];
time=report.Nominal{index}.Time(1:101);
[~,Y]=ode45(@(t,y)augmentedODE(t,y,plant,c),time,y0, ...
    odeset('RelTol',1e-10,'AbsTol',1e-11));
checks.ODEResidual=max(abs(Y(:,1:n)-report.Nominal{index}.State(1:101,:)),[],'all');
checks.IntegratorResidual=max(abs(Y(:,n+1:end)-report.Nominal{index}.ControllerState(1:101,:)),[],'all');
x=plant.x_s+e0; h=sin((1:plant.N)');
z=c.AbsoluteIntegratorEquilibrium+c.TimeScale*h;
absoluteCommand=sum(c.PhysicalLocalGains.*[x(1:2:end),x(2:2:end),z],2);
errorCommand=plant.u_s+c.CommandGain*e0+c.IntegralGain*h;
checks.AbsoluteLawResidual=norm(absoluteCommand-errorCommand,inf);
A=report.Nominal{index}.ClosedLoopMatrix;
d=repmat([0.35;0],plant.N,1);
loadEquilibrium=-A\[plant.BBar*d;zeros(plant.N,1)];
checks.ConstantLoadVoltageResidual=norm(loadEquilibrium(1:2:n),inf);
checks.ConstantLoadCurrentResidual=norm(loadEquilibrium(2:2:n)-0.35,inf);
checks.MatchedDisturbances=true;
for scenario={'Disturbance','BoundedNoise','LoadStep'}
    runs=report.(scenario{1});
    for k=2:numel(runs)
        checks.MatchedDisturbances=checks.MatchedDisturbances && ...
            isequal(runs{1}.PrimitiveDisturbance,runs{k}.PrimitiveDisturbance) && ...
            isequal(runs{1}.StateJumps,runs{k}.StateJumps);
    end
end
checks.PeakBudgetSatisfied=report.Nominal{index}.Metrics.PeakCommandDeviation ...
    <=report.LiteratureTuning.PeakCommandBudget*(1+1e-10);
checks.Passed=c.Passivity.Accepted && checks.ODEResidual<1e-6 && ...
    checks.IntegratorResidual<1e-6 && checks.AbsoluteLawResidual<1e-10 && ...
    checks.ConstantLoadVoltageResidual<1e-10 && checks.ConstantLoadCurrentResidual<1e-10 ...
    && checks.MatchedDisturbances && checks.PeakBudgetSatisfied;
fprintf('Literature PI: ODE error %.3e, absolute-law error %.3e, checks passed=%d.\n', ...
    checks.ODEResidual,checks.AbsoluteLawResidual,checks.Passed);
end

function dy=augmentedODE(t,y,plant,c)
n=2*plant.N; x=y(1:n); h=y(n+1:end);
commandDeviation=c.CommandGain*(x-plant.x_s)+c.IntegralGain*h;
dx=plant.dynamics(t,x,false)+plant.B*commandDeviation;
dh=(plant.x_s(1:2:end)-x(1:2:end))/c.TimeScale;
dy=[dx;dh];
end

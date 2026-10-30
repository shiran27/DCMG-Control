function checks = verifyPaperImplementation(net,report)
%VERIFYPAPERIMPLEMENTATION Independent physics, implementation, and data gates.
checks=struct();
checks.EquilibriumResidual=net.ss.PhysicsResidual;
checks.LaplacianResidual=norm(net.YBar*ones(net.N,1),inf);
checks.LineCurrentResidual=0;
for j=1:net.M
    line=net.Lines(j); V=net.x_s(1:2:end);
    checks.LineCurrentResidual=max(checks.LineCurrentResidual, ...
        abs(line.current(V(line.i),V(line.j))-(V(line.i)-V(line.j))*line.g));
end
numberOfMethods=numel(report.Controllers);
checks.ODEResidual=zeros(1,numberOfMethods);
for k=1:numberOfMethods
    [plant,~]=createPaperDCMicrogrid(LegacyEquilibrium=false,FinalTime=0.12);
    c=report.Controllers{k};
    % Install the total command gain in one place. This is algebraically
    % equivalent to the separate local/global representation and also
    % covers the purely local conventional-droop baseline.
    plant.K=c.CommandGain;
    for i=1:plant.N, plant.DGs(i).K=[0,0]; end
    time=report.Nominal{k}.Time(1:51);
    initial=report.Nominal{k}.State(1,:)';
    [~,X]=ode45(@(t,x)plant.dynamics(t,x,false),time,initial, ...
        odeset('RelTol',1e-10,'AbsTol',1e-11));
    checks.ODEResidual(k)=max(abs(X-report.Nominal{k}.State(1:51,:)),[],'all');
end
assert(max(checks.ODEResidual)<1e-6,'DCMG:ODECrosscheck', ...
    'Existing simulator and exact linear propagation disagree.');
% A command-only hold must preserve the continuously coupled physical plant.
plant.Ts=1e-4; K=report.Controllers{3}.CommandGain;
plant.K=K;
for i=1:plant.N, plant.DGs(i).K=[0,0]; end
initial=report.Nominal{3}.State(1,:)';
[t,X]=plant.simulateStageZOH(0,0.001,initial,1e-4);
A=plant.A+plant.E*plant.YBar*plant.D'; n=size(A,1);
Z=expm(plant.Ts*[A,plant.B;zeros(plant.N,n+plant.N)]);
F=Z(1:n,1:n)+Z(1:n,n+1:end)*K;
pred=plant.x_s+F^10*(initial-plant.x_s);
checks.ZOHResidual=norm(X(end,:)'-pred,inf);
assert(abs(t(end)-0.001)<1e-12 && checks.ZOHResidual<1e-5, ...
    'DCMG:ZOHCrosscheck','Held-command propagation cross-check failed.');
% DD local synthesis must be invariant to changing the withheld truth A/B.
dg=net.DGs(1); previousA=dg.A; previousB=dg.BBar;
dg.DesignOptions.Rho=report.SelectedRho;
reference=designLocalIntegralDissipativity(dg,true);
dg.A=previousA+17*eye(2); dg.BBar=3*previousB;
blind=designLocalIntegralDissipativity(dg,true);
dg.A=previousA; dg.BBar=previousB;
assert(reference.problem==0 && blind.problem==0,'DCMG:BlindSynthesis');
checks.WithheldModelGainDifference=norm(reference.K-blind.K,inf);
assert(checks.WithheldModelGainDifference<1e-9,'DCMG:ModelLeakage');
% Reject rank-deficient input before accepting any controller.
original=dg.DataContract; dg.DataContract.Rank=3; caught=false;
try, designLocalIntegralDissipativity(dg,true);
catch e, caught=strcmp(e.identifier,'DCMG:DataRegularity'); end
dg.DataContract=original;
checks.RankFailureRejected=caught;
assert(caught,'DCMG:RankGuard');
% Published unknown-B lift obstruction: a positive QMI direction precludes
% lambda>0 because the final principal block is -lambda*N.
N=diag([1,-1]); checks.UnknownBPrincipalBlockEigenvalue=min(eig(-N));
checks.UnknownBLiftUsable=false;
checks.Passed=checks.EquilibriumResidual<1e-7 && checks.LaplacianResidual<1e-10 ...
    && checks.LineCurrentResidual<1e-10 && caught;
fprintf('Existing ODE vs exact propagation: %.3e; held-command check: %.3e\n', ...
    max(checks.ODEResidual),checks.ZOHResidual);
end

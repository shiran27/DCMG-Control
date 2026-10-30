function checks = test_conventional_droop_baseline()
%TEST_CONVENTIONALDROOPBASELINE Solver-free checks of the literature baseline.
root=fileparts(fileparts(mfilename('fullpath')));
addpath(root,fullfile(root,'support'));
[net,~]=createPaperDCMicrogrid(LegacyEquilibrium=false,FinalTime=0.12);
fractions=[0.01,0.03,0.05,0.10];
numberOfCases=numel(fractions);
maxRealPole=zeros(numberOfCases,1);
settlingMs=zeros(numberOfCases,1);
voltageRMS=zeros(numberOfCases,1);
sharingRMS=zeros(numberOfCases,1);
perturbation=repmat([1.5;0.35],net.N,1).* ...
    [1;-1;-1;1;1;1;-1;-1;1;-1;-1;1];

for k=1:numberOfCases
    controller=designConventionalDroop(net, ...
        VoltageDroopFraction=fractions(k));
    expected=fractions(k)*[net.DGs.Vrated]'./[net.DGs.Irated]';
    assert(norm(controller.DroopResistance-expected,inf)<1e-12);
    assert(controller.DirectedEdges==0 && ...
        nnz(controller.CommandGain(:,1:2:end))==0);
    assert(norm(controller.NoLoadVoltage- ...
        controller.DroopResistance.*net.x_s(2:2:end)-net.u_s,inf)<1e-12);
    assert(all(diag(controller.CommandGain(:,2:2:end))<0));

    nominal=simulatePaperCase(net,controller.CommandGain, ...
        InitialState=net.x_s+perturbation,Certificate=controller);
    disturbed=simulatePaperCase(net,controller.CommandGain, ...
        InitialState=net.x_s,Scenario="disturbance",Certificate=controller);
    assert(nominal.Accepted && disturbed.Accepted && ...
        nominal.Metrics.FinalRatio<1e-5);
    maxRealPole(k)=nominal.Metrics.MaxRealPole;
    settlingMs(k)=1e3*nominal.Metrics.SettlingTime2Percent;
    voltageRMS(k)=disturbed.Metrics.VoltageRMS;
    sharingRMS(k)=disturbed.Metrics.SharingRMS;
end

checks=struct();
checks.Passed=true;
checks.Formula='V_ti = V_0i - R_di I_ti';
checks.SolverIndependent=true;
checks.Table=table(100*fractions(:),maxRealPole,settlingMs,voltageRMS, ...
    sharingRMS,'VariableNames',{'DroopPercent','MaxRealPole','SettlingMs', ...
    'DisturbanceVoltageRMS','DisturbanceSharingRMS'});
disp(checks.Table);
fprintf('Conventional droop contract and all declared sensitivity cases passed.\n');
end

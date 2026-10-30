function checks = test_regenerated_results()
%TEST_REGENERATED_RESULTS Verify serialized certificates and failure rollback.
root=fileparts(fileparts(mfilename('fullpath')));
addpath(root,fullfile(root,'support'));
s=load(fullfile(root,'Results','Regenerated','paper_results.mat'),'report');
r=s.report;
assert(r.Accepted && r.ImplementationChecks.Passed);
[net,~]=createPaperDCMicrogrid(LegacyEquilibrium=false);
for k=2:3
    c=r.Controllers{k};
    aggregate=verifyGlobalStorage(net,c);
    assert(aggregate.Accepted && c.OriginalDiagnostics.Accepted);
    for i=1:6
        local=verifyLocalDissipativity(net.DGs(i),c.Local{i});
        assert(local.Accepted);
    end
    assert(c.PhysicalResidual<1e-9 && c.DirectedEdges==nnz(c.Adjacency));
    assert(r.BoundedNoise{k}.Metrics.EmpiricalGain<c.Gamma);
end
droopIndex=find(startsWith(r.Methods,"Conventional droop"),1);
assert(~isempty(droopIndex),'DCMG:MissingDroopBaseline', ...
    'The serialized report does not contain the droop baseline.');
droop=r.Controllers{droopIndex};
assert(droop.Accepted && droop.DirectedEdges==0 && ...
    nnz(droop.CommandGain(:,1:2:end))==0 && ...
    droop.EquilibriumCommandResidual<1e-12);
expected=r.Options.DroopFraction*[net.DGs.Vrated]'./[net.DGs.Irated]';
assert(norm(droop.DroopResistance-expected,inf)<1e-12, ...
    'DCMG:DroopCoefficientContract');
assert(all(r.DroopSensitivity.Table.Accepted) && ...
    height(r.DroopSensitivity.Table)==numel(r.Options.DroopSensitivity));
for i=1:6
    d=r.Data.Local{i};
    assert(d.Rank==4 && d.SlaterMargin>0 && d.TrueModelMargin>0 && d.BoundSlack>0);
    assert(d.EquationResidual<1e-9);
end
net.K=r.Controllers{1}.CommandGain;
old=captureControllerState(net);
net.DesignOptions.Rho=0.5;
[~,~,failed]=net.codesign_MB_DRC(0);
assert(failed.problem~=0,'DCMG:FailureCaseChanged','Revisit the intentionally rejected gain case.');
after=captureControllerState(net);
assert(isequaln(old,after),'DCMG:RollbackFailed','Failed design changed the installed controller.');
% No zero-gain network may be drawn as a fully connected graph.
net.K=zeros(net.N,2*net.N); [adj,~]=net.buildCommAdjFromK(0);
assert(nnz(adj)==0);
% Parallel conductances and mutable line resistance stay physically consistent.
l=TransmissionLine(99,1,2,2); l.R=4; assert(l.g==0.25);
checks=struct('Passed',true,'SerializedCertificates',true,'FailedDesignRollback',true, ...
    'ZeroGainTopology',true,'DependentConductance',true);
save(fullfile(root,'Results','Regenerated','regression_checks.mat'),'checks');
fprintf('Serialized certificates, data, failure rollback, and topology checks passed.\n');
end

function out = designLocalIntegralDissipativity(dg,useData)
%DESIGNLOCALINTEGRALDISSIPATIVITY Hard CT local Xi certificate, MB or integral DD.
% Variables use t/TimeScale and x/StateScale, u/InputScale. The DD branch
% receives only integral regression data and a declared disturbance bound.
% The S-procedure is sufficient: substitute [e;A'*e;B'*e] into its hard LMI.
opts=paperDesignOptions();
names=fieldnames(dg.DesignOptions);
for k=1:numel(names), opts.(names{k})=dg.DesignOptions.(names{k}); end
n=2; I=eye(n); O=zeros(n); eps0=opts.Margin;
P=sdpvar(n,n,'symmetric');
L=[zeros(1,n);sdpvar(1,n,'full')];
a=sdpvar(1); rho=opts.Rho; c=0.5;
cons=[P>=eps0*I,P<=100*I,a>=eps0,a<=100];
% Convex gain surrogate; the recovered gain bound is also checked.
cons=[cons,-opts.LocalGainLimit*I<=L+L',L+L'<=opts.LocalGainLimit*I];
if useData
    data=dg.DataContract;
    assert(isfield(data,'Kind') && strcmp(data.Kind,'continuous-integral'), ...
        'DCMG:DataContract','Use collectIntegralData before DD synthesis.');
    N=data.QMI;
    assert(data.Rank==4 && data.SlaterMargin>0,'DCMG:DataRegularity', ...
        'Integral data require full rank and strict QMI interior.');
    lambda=sdpvar(1);
    % Order [e; A''e; B''e; local disturbance; output Schur variable].
    M=[O,-P,-L',c*P-I,P; ...
       -P,O,O,O,O; -L,O,O,O,O; ...
       c*P-I,O,O,a*I,O; P,O,O,O,I/rho];
    M=M-lambda*blkdiag(N,zeros(4));
    cons=[cons,lambda>=0];
    % Fixed congruence whitens the observed regression directions. It does
    % not alter the model set or introduce model matrices into synthesis.
    T=blkdiag(data.Congruence,eye(4));
    W=T'*M*T;
else
    A=opts.TimeScale*dg.A;
    B=opts.TimeScale*dg.BBar*opts.InputScale/opts.StateScale;
    APBL=A*P+B*L;
    W=[-APBL-APBL',c*P-I,P; c*P-I,a*I,O; P,O,I/rho];
end
cons=[cons,W>=eps0*eye(size(W))];
ops=dcmgSolverOptions();
sol=optimize(cons,a+opts.LocalRegularization*(trace(P)+norm(L(:),1)),ops);
out=struct('problem',sol.problem,'info',sol.info,'Solver',ops.solver, ...
    'Options',opts,'DataDriven',useData);
if sol.problem~=0, return; end
out.P=value(P); out.L=value(L); out.a=value(a); out.rho=rho;
out.nu=-out.a; out.FullGain=out.L/out.P;
out.K=opts.InputScale/opts.StateScale*out.FullGain(2,:);
out.Storage=inv(out.P);
out.LMI=value(W); out.LMIDiagnostics=certificateDiagnostics(out.LMI);
out.Supply=[out.a*I,c*I;c*I,-rho*I];
out.GainRecoveryResidual=norm(out.FullGain*out.P-out.L,'fro');
out.Accepted=out.LMIDiagnostics.Accepted && all(isfinite(out.K)) ...
    && max(abs(out.K))<=opts.LocalGainLimit && out.a>0;
if useData
    out.Multiplier=value(lambda);
    out.DataSummary=rmfield(data,{'X','U','Delta','PrimitiveNoise'});
    out.Accepted=out.Accepted && out.Multiplier>=0;
end
if ~out.Accepted
    out.problem=4; out.info='Local certificate or recovered gain rejected';
    return
end
dg.K=out.K; dg.nu=out.nu; dg.rho=out.rho;
dg.DesignCertificate=out;
end

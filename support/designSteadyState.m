function ss = designSteadyState(net,opts)
%DESIGNSTEADYSTATE Bounded equilibrium with exact proportional sharing.
% The equilibrium is affine in one common current-utilization fraction.
% Eliminate all equality constraints before solving the scalar convex QP.
if nargin<2, opts=struct(); end
defaults=struct('VoltageTolerance',0.2,'CommandTolerance',0.2, ...
    'LineCurrentFraction',1,'CurrentUpper',0.98,'CurrentTarget',0.55, ...
    'VoltageWeight',1,'CurrentWeight',0.05,'LegacyObjective',false);
f=fieldnames(defaults);
for k=1:numel(f)
    if ~isfield(opts,f{k}), opts.(f{k})=defaults.(f{k}); end
end
N=net.N; Vr=[net.DGs.Vrated]'; Ir=[net.DGs.Irated]';
RL=[net.DGs.RL]'; Rf=[net.DGs.Rf]'; loads=[net.DGs.Ibar]';
G=net.YBar+diag(1./RL);
vSlope=G\Ir; vOffset=-G\loads;
uSlope=vSlope+Rf.*Ir; uOffset=vOffset;
lSlope=net.YBar*vSlope; lOffset=net.YBar*vOffset;
lo=0; hi=opts.CurrentUpper;
slopes=[vSlope;uSlope;lSlope]; offsets=[vOffset;uOffset;lOffset];
lower=[(1-opts.VoltageTolerance)*Vr;(1-opts.CommandTolerance)*Vr; ...
    -opts.LineCurrentFraction*Ir];
upper=[(1+opts.VoltageTolerance)*Vr;(1+opts.CommandTolerance)*Vr; ...
    opts.LineCurrentFraction*Ir];
for k=1:numel(slopes)
    if abs(slopes(k))<1e-12
        assert(offsets(k)>=lower(k) && offsets(k)<=upper(k),'DCMG:EquilibriumBounds');
    else
        ends=sort(([lower(k),upper(k)]-offsets(k))/slopes(k));
        lo=max(lo,ends(1)); hi=min(hi,ends(2));
    end
end
assert(lo<=hi,'DCMG:EquilibriumInfeasible','Equilibrium rating constraints are infeasible.');
a=vSlope./Vr; b=vOffset./Vr-1;
if opts.LegacyObjective
    alpha=-(2*(a'*b)+1)/(2*(a'*a));
else
    alpha=(opts.CurrentWeight*opts.CurrentTarget-opts.VoltageWeight*(a'*b))/ ...
        (opts.VoltageWeight*(a'*a)+opts.CurrentWeight);
end
alpha=min(hi,max(lo,alpha));
V=vSlope*alpha+vOffset; It=alpha*Ir; command=uSlope*alpha+uOffset;
x=zeros(2*N,1); x(1:2:end)=V; x(2:2:end)=It;
ss=struct('problem',0,'info','Exact scalar convex equilibrium solve', ...
    'x_ss',x,'u_s',command,'epsV',V./Vr,'epsI',alpha, ...
    'V_ss',V,'I_tss',It,'I_ss',net.YBar*V,'Options',opts, ...
    'FeasibleUtilizationInterval',[lo,hi]);
ss.PhysicsResidual=norm(net.A*x+net.E*ss.I_ss+net.B*command+net.E*loads,inf);
ss.LineResidual=norm(ss.I_ss-net.YBar*V,inf);
ss.BoundSlack=min([V-lower(1:N);upper(1:N)-V; ...
    command-lower(N+1:2*N);upper(N+1:2*N)-command; ...
    ss.I_ss-lower(2*N+1:end);upper(2*N+1:end)-ss.I_ss]);
assert(ss.PhysicsResidual<1e-7 && ss.LineResidual<1e-9 && ss.BoundSlack>=-1e-9, ...
    'DCMG:EquilibriumResidual','Equilibrium verification failed.');
net.x_s=x; net.u_s=command; net.ss=ss;
for i=1:N
    net.DGs(i).x_s=x(2*i-1:2*i); net.DGs(i).u_s=command(i);
    net.DGs(i).I_s=ss.I_ss(i);
end
end

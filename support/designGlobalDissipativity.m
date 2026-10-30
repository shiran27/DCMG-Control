function out = designGlobalDissipativity(net)
%DESIGNGLOBALDISSIPATIVITY Hard known-B global composition of local supplies.
% Physical line-current rows are fixed by the resistor network. Only the
% converter-command rows are optimized. No certificate slack or post-pruning.
opts=paperDesignOptions(); names=fieldnames(net.DesignOptions);
for k=1:numel(names), opts.(names{k})=net.DesignOptions.(names{k}); end
N=net.N; n=2*N; I=eye(n); O=zeros(n); e=opts.Margin;
p=sdpvar(N,1); H=sdpvar(N,n,'full'); gammaSq=sdpvar(1);
P=kron(diag(p),eye(2));
avec=zeros(N,1); rvec=zeros(N,1); localGain=zeros(N,n);
for i=1:N
    cert=net.DGs(i).DesignCertificate;
    assert(isfield(cert,'Accepted') && cert.Accepted, ...
        'DCMG:LocalCertificate','DG %d has no accepted local certificate.',i);
    avec(i)=cert.a; rvec(i)=cert.rho;
    localGain(i,2*i-1:2*i)=net.DGs(i).K;
end
A11=kron(diag(avec),eye(2));
Xp11=A11*P; Xp12=0.5*P; Xp22=-kron(diag(rvec),eye(2))*P;
B=opts.TimeScale*net.BBar*opts.InputScale/opts.StateScale;
physical=net.D*net.YBar*net.D'*opts.StateScale/opts.InputScale;
% Row scaling p_i commutes with the local block diagonal B matrix.
Kbar=P*physical+net.DBar*H;
BK=B*Kbar;
W=[Xp11,O,A11*BK,Xp11; ...
   O,I,I,O; ...
   BK'*A11,I,-0.5*(BK+BK')-Xp22,-Xp12; ...
   Xp11,O,-Xp12,gammaSq*I];
allowed=opts.AllowedEdges;
if isempty(allowed), allowed=true(N); end
assert(isequal(size(allowed),[N,N]));
weights=zeros(N,n);
for i=1:N
    for j=1:N
        ix=2*j-1:2*j;
        if i~=j, weights(i,ix)=norm(net.DGs(i).pos-net.DGs(j).pos); end
    end
end
weights=weights/max(1,max(weights(:)));
cost=norm(reshape(weights.*H,[],1),1);
cons=[p>=1e-4,p<=1e4,gammaSq>=e,gammaSq<=1e6, ...
    W>=e*eye(4*n),cost<=opts.CommunicationBudget];
for i=1:N
    cons=[cons,-opts.GlobalGainLimit*p(i)<=H(i,:),H(i,:)<=opts.GlobalGainLimit*p(i)];
    for j=1:N
        if ~allowed(i,j), cons=[cons,H(i,2*j-1:2*j)==0]; end
    end
end
ops=dcmgSolverOptions();
sol=optimize(cons,opts.PerformanceWeight*gammaSq+opts.EdgePenalty*cost ...
    +1e-5*sum(p),ops);
out=struct('problem',sol.problem,'info',sol.info,'Solver',ops.solver,'Options',opts);
if sol.problem~=0, return; end
out.Weights=value(p); out.TransformedGain=value(H);
out.K=diag(out.Weights)\out.TransformedGain;
% Explicit structural zeros remove only solver roundoff on imposed equalities.
for i=1:N
    for j=1:N
        if ~allowed(i,j), out.K(i,2*j-1:2*j)=0; end
    end
end
out.KFull=physical+net.DBar*out.K;
out.GammaSquared=value(gammaSq); out.Gamma=sqrt(out.GammaSquared);
out.WeightedCost=value(cost); out.LMI=value(W);
out.LMIDiagnostics=certificateDiagnostics(out.LMI);
Pv=kron(diag(out.Weights),eye(2));
Xp=[A11*Pv,0.5*Pv;0.5*Pv,-kron(diag(rvec),eye(2))*Pv];
M=[B*out.KFull,I;I,O];
out.OriginalResidual=blkdiag(-I,out.GammaSquared*I)-M'*Xp*M;
out.OriginalDiagnostics=certificateDiagnostics(out.OriginalResidual);
out.PhysicalResidual=norm(net.D'*out.KFull-net.YBar*net.D'* ...
    opts.StateScale/opts.InputScale,'fro');
out.CommandGain=localGain+opts.InputScale/opts.StateScale*out.K;
Acl=net.A+net.E*net.YBar*net.D'+net.B*out.CommandGain;
out.MaxRealPole=max(real(eig(Acl)));
out.Accepted=out.LMIDiagnostics.Accepted && out.OriginalDiagnostics.Accepted ...
    && out.PhysicalResidual<1e-9 && out.MaxRealPole<0;
out.Adjacency=false(N); out.GeometricCost=0;
for i=1:N
    for j=1:N
        block=out.K(i,2*j-1:2*j);
        out.Adjacency(i,j)=i~=j && norm(block)>1e-9;
        if i~=j, out.GeometricCost=out.GeometricCost+weights(i,2*j-1)*norm(block,1); end
    end
end
out.DirectedEdges=nnz(out.Adjacency);
out.UndirectedPairs=nnz(triu(out.Adjacency|out.Adjacency',1));
out.Scope='Known continuous input matrices; exact physical resistor interconnection';
if ~out.Accepted
    out.problem=4; out.info='Independent global certificate verification failed';
    return
end
net.K=opts.InputScale/opts.StateScale*out.K;
net.commAdj=out.Adjacency;
net.DesignCertificate=out;
end

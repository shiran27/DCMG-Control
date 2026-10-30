function data = collectIntegralData(net,options)
%COLLECTINTEGRALDATA Exact measured-integral regression on the existing plant.
% Delta=(x(t+h)-x(t))/h, X=average(x), U=average([Iline;Vt]).
% Thus Delta=A*X+B*U+W holds exactly even though line current is not held.
% Only the generator and independent validation use the simulation model.
arguments
    net
    options.Seed (1,1) double = 2718
    options.Samples (1,1) double = 400
    options.SampleTime (1,1) double = 1e-4
    options.Excitation (1,1) double = 0.25
    options.PrimitiveNoiseBound (1,1) double = 1e-4
    options.InputMatrixNormBound (1,1) double = 1.25
end
rng(options.Seed,'twister');
o=paperDesignOptions(); n=2*net.N; h=options.SampleTime/o.TimeScale;
A=o.TimeScale*(net.A+net.E*net.YBar*net.D');
B=o.TimeScale*net.B*o.InputScale/o.StateScale;
Bw=o.TimeScale*net.BBar*o.InputScale/o.StateScale;
% Augment with integral state qdot=x. Both command and primitive noise are
% held; the resistor network remains inside A throughout each interval.
G=[B,Bw]; m=size(G,2);
Z=expm(h*[A,zeros(n),G;eye(n),zeros(n,n+m);zeros(m,2*n+m)]);
F=Z(1:n,1:n); J=Z(1:n,2*n+1:end);
Fx=Z(n+1:2*n,1:n)/h; Jx=Z(n+1:2*n,2*n+1:end)/h;
T=options.Samples; states=zeros(n,T+1); averages=zeros(n,T);
commands=zeros(net.N,T); noise=zeros(n,T);
states(:,1)=0.05*(2*rand(n,1)-1);
for k=1:T
    % Fixed dissipative damping plus independent bounded PE, no gain corruption.
    commands(:,k)=-states(2:2:end,k)+options.Excitation*(2*rand(net.N,1)-1);
    noise(:,k)=options.PrimitiveNoiseBound/sqrt(2)*(2*rand(n,1)-1);
    v=[commands(:,k);noise(:,k)];
    states(:,k+1)=F*states(:,k)+J*v;
    averages(:,k)=Fx*states(:,k)+Jx*v;
end
delta=diff(states,1,2)/h;
lineAverages=net.YBar*net.D'*averages*o.StateScale/o.InputScale;
eta=T*(options.InputMatrixNormBound*options.PrimitiveNoiseBound)^2;
data=struct('Options',options,'Normalization',o,'States',states, ...
    'StateAverages',averages,'Command',commands,'PrimitiveNoise',noise, ...
    'Time',(0:T)*options.SampleTime,'Local',{{}},'Contract', ...
    'Exact continuous integral data; measured states/inputs, noise-free measurements');
for i=1:net.N
    ix=2*i-1:2*i; X=averages(ix,:); U=[lineAverages(i,:);commands(i,:)];
    D=delta(ix,:); R=[X;U]; sv=svd(R);
    regRank=sum(sv>sv(1)*1e-9);
    assert(regRank==4,'DCMG:DataRank','DG %d lacks excitation.',i);
    % Qw=diag(eta*I,-I_T) chosen BEFORE looking at realized disturbances.
    QMI=[eta*eye(2)-D*D',D*R';R*D',-R*R'];
    QMI=(QMI+QMI')/2;
    center=(R*R')\(R*D');
    slater=QMI(1:2,1:2)-QMI(1:2,3:end)*(QMI(3:end,3:end)\QMI(3:end,1:2));
    % Completing the square diagonalizes N; whitening avoids cancellations
    % between large regression terms and a small declared energy reserve.
    Tq=[eye(2),zeros(2,4);center,sqrt(eta)*inv(chol(R*R'))];
    local=struct('Kind','continuous-integral','X',X,'U',U,'Delta',D, ...
        'PrimitiveNoise',noise(ix,:),'QMI',QMI/eta,'QMIScale',eta,'Congruence',Tq, ...
        'EnergyBound',eta,'Rank',regRank,'SingularValues',sv, ...
        'SlaterMargin',min(eig((slater+slater')/2)), ...
        'OutputResidual',0,'BoundSource', ...
        'T*(declared input-matrix norm bound * primitive amplitude bound)^2');
    assert(local.SlaterMargin>0,'DCMG:NoSlater','No strict local QMI interior.');
    net.DGs(i).DataContract=local;
    % Truth is used only here for independent generator/bound verification.
    Ai=o.TimeScale*net.DGs(i).A;
    Bi=o.TimeScale*net.DGs(i).BBar*o.InputScale/o.StateScale;
    assert(norm(Bi,2)<=options.InputMatrixNormBound,'DCMG:BoundContract', ...
        'Declared input-matrix norm bound violated.');
    W=Bi*noise(ix,:);
    truth=[eye(2);Ai';Bi'];
    local.EquationResidual=norm(D-Ai*X-Bi*U-W,'fro');
    local.TrueModelMargin=min(eig((truth'*QMI*truth+truth'*QMI'*truth)/2));
    local.BoundSlack=min(eig(eta*eye(2)-W*W'));
    assert(local.EquationResidual<1e-9 && local.TrueModelMargin>0, ...
        'DCMG:IntegralDataValidation','Exact data equation or bound failed.');
    data.Local{i}=local;
end
end

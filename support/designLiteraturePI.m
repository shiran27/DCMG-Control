function c = designLiteraturePI(net,options)
%DESIGNLITERATUREPI Nahata et al., Automatica 113 (2020), 108770.
% Published multivariable PI law (18)-(19), gain region (39), storage (38).
% DOI: 10.1016/j.automatica.2019.108770. Specialization to the existing ZI
% plant and static resistive lines; no converter or line model is replaced.
% z integrates Vref-V in V*s. Internally h=(z-z_ss)/TimeScale is used.
arguments
    net
    options.VoltageGain (1,1) double {mustBePositive} = 1
    options.CurrentGain (1,1) double {mustBePositive} = 1
    options.IntegralFraction (1,1) double {mustBePositive} = 0.5
    options.TimeScale (1,1) double {mustBePositive} = 1e-3
    options.GainLimit (1,1) double {mustBePositive} = 20
end
assert(options.IntegralFraction<1,'DCMG:PIIntegralFraction');
N=net.N; n=2*N; t0=options.TimeScale;
Kx=zeros(N,n); Kz=zeros(N); P=zeros(3*N);
local=cell(1,N); physicalGains=zeros(N,3); margins=zeros(N,4);
for i=1:N
    dg=net.DGs(i); k1=-options.VoltageGain; k2=-options.CurrentGain;
    k3=options.IntegralFraction*(k1-1)*(k2-dg.Rf)/dg.L;
    alpha=(k1-1)/dg.L; beta=(k2-dg.Rf)/dg.L; gamma=k3/dg.L;
    omega=gamma-alpha*beta;
    Pphysical=[dg.C,0,0;0,beta/omega,gamma/omega; ...
        0,gamma/omega,alpha*gamma/omega];
    T=diag([1,1,t0]); Pi=T'*Pphysical*T;
    Fi=[-1/(dg.C*dg.RL),1/dg.C,0;alpha,beta,gamma*t0;-1/t0,0,0];
    Q=-(Pi*Fi+Fi'*Pi);
    portResidual=norm(Pi*[1/dg.C;0;0]-[1;0;0],inf);
    margins(i,:)=[1-k1,dg.Rf-k2,k3,(k1-1)*(k2-dg.Rf)/dg.L-k3];
    local{i}=struct('Storage',Pi,'Dissipation',Q,'PortResidual',portResidual, ...
        'StorageMinEig',min(eig(Pi)),'DissipationMinEig',min(eig((Q+Q')/2)));
    assert(all(margins(i,:)>0) && local{i}.StorageMinEig>0 && ...
        local{i}.DissipationMinEig>-1e-8 && portResidual<1e-10,'DCMG:PIPassivity');
    ix=[2*i-1,2*i,n+i]; P(ix,ix)=Pi;
    Kx(i,2*i-1:2*i)=[k1,k2]; Kz(i,i)=t0*k3;
    physicalGains(i,:)=[k1,k2,k3];
end
Aphysical=net.A+net.E*net.YBar*net.D';
A=[Aphysical+net.B*Kx,net.B*Kz;-net.D'/t0,zeros(N)];
Q=-(P*A+A'*P); Q=(Q+Q')/2;
gainNorm=max(vecnorm([Kx,Kz],2,2));
passivity=struct('Local',{local},'Storage',P,'Dissipation',Q, ...
    'StorageMinEig',min(eig(P)),'DissipationMinEig',min(eig(Q)), ...
    'MaxRealPole',max(real(eig(A))),'GainConditionMargins',margins);
passivity.Accepted=passivity.StorageMinEig>0 && passivity.DissipationMinEig>-1e-8 ...
    && passivity.MaxRealPole<0;
c=struct('Type','literature-pi','CommandGain',Kx,'IntegralGain',Kz, ...
    'PhysicalLocalGains',physicalGains,'TimeScale',t0,'Options',options, ...
    'Gamma',NaN,'DirectedEdges',0,'Adjacency',false(N),'Passivity',passivity, ...
    'GainNorm',gainNorm,'problem',0,'info','Verified published PI gain region and storage');
if ~passivity.Accepted || gainNorm>options.GainLimit
    c.problem=1; c.info='PI certificate or normalized gain limit rejected';
end
% Translating the published absolute-state law into error coordinates is exact.
c.AbsoluteIntegratorEquilibrium=(net.u_s-Kx*net.x_s)./physicalGains(:,3);
c.Citation=struct('Key','Nahata2020','DOI','10.1016/j.automatica.2019.108770', ...
    'URL','https://doi.org/10.1016/j.automatica.2019.108770', ...
    'Equations','(18)-(19), (38)-(39) in the author postprint', ...
    'Postprint','https://iris.unipv.it/retrieve/dac43a76-0142-4bad-867e-9972782278ac/A%20passivity%20based%20-%20postprint.pdf');
end

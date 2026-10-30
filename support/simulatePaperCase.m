function result = simulatePaperCase(net,K,options)
%SIMULATEPAPERCASE Exact propagation of the existing DCMicrogrid linear model.
% The resistor network is always continuous. Disturbances are piecewise held.
% Verification against DCMicrogrid.dynamics is performed by the regression test.
arguments
    net
    K
    options.FinalTime (1,1) double = 0.12
    options.OutputStep (1,1) double = 1e-4
    options.InitialState double = []
    options.Scenario (1,1) string = "nominal"
    options.Seed (1,1) double = 31415
    options.Certificate struct = struct()
end
rng(options.Seed,'twister'); n=2*net.N;
Aphysical=net.A+net.E*net.YBar*net.D';
if isstruct(K)
    assert(strcmp(K.Type,'literature-pi'),'DCMG:DynamicControllerType');
    Kx=K.CommandGain; Kz=K.IntegralGain;
    A=[Aphysical+net.B*Kx,net.B*Kz;-net.D'/K.TimeScale,zeros(net.N)];
    Bdist=[net.BBar;zeros(net.N,n)]; nz=net.N;
else
    Kx=K; Kz=zeros(net.N,0); A=Aphysical+net.B*Kx;
    Bdist=net.BBar; nz=0;
end
dt=options.OutputStep; steps=round(options.FinalTime/dt);
assert(abs(steps*dt-options.FinalTime)<1e-12,'DCMG:TimeGrid');
t=(0:steps)'*dt; x0=options.InitialState;
if isempty(x0), x0=net.x_s; end
augmentedError=zeros(steps+1,n+nz); augmentedError(1,1:n)=(x0-net.x_s)';
primitive=zeros(steps+1,n); kicks=zeros(steps+1,n);
if options.Scenario=="disturbance"
    % Same deterministic pulse and symmetric two-state jump for every method.
    pulse=t>=0.025 & t<0.035;
    amp=repmat([0.35,0.5],1,net.N);
    primitive(pulse,:)=repmat(amp.*sign(2*rand(1,n)-1),nnz(pulse),1);
    ix=round(0.065/dt)+1;
    kicks(ix,:)=repmat([1.0,0.25],1,net.N).*(2*rand(1,n)-1);
elseif options.Scenario=="bounded-noise"
    primitive=repmat([0.1,0.15],steps+1,net.N).*(2*rand(steps+1,n)-1);
elseif options.Scenario=="load-step"
    % A persistent +0.35 A load increment at every bus, from 25 ms.
    primitive(t>=0.025,1:2:end)=0.35;
end
na=n+nz; Z=expm(dt*[A,Bdist;zeros(n,na+n)]);
F=Z(1:na,1:na); G=Z(1:na,na+1:end);
for k=1:steps
    augmentedError(k+1,:)=(F*augmentedError(k,:)'+G*primitive(k,:)')';
    augmentedError(k+1,1:n)=augmentedError(k+1,1:n)+kicks(k+1,:);
end
error=augmentedError(:,1:n); controllerState=augmentedError(:,n+1:end);
X=error+net.x_s'; V=X(:,1:2:end); It=X(:,2:2:end);
U=error*Kx'+controllerState*Kz'+net.u_s'; lineCurrent=V*net.YBar';
Ir=[net.DGs.Irated]; Vr=[net.DGs.Vrated];
scaled=error/48; normError=vecnorm(scaled,2,2);
sharing=It./Ir-mean(It./Ir,2);
voltageError=V-net.x_s(1:2:end)';
ratio=normError(end)/max(normError(1),eps);
threshold=0.02*max(normError(1),eps);
lastOutside=find(normError>threshold,1,'last');
settle=NaN;
if isempty(lastOutside), settle=0;
elseif lastOutside<numel(t), settle=t(lastOutside+1); end
metrics=struct('MaxRealPole',max(real(eig(A))),'FinalRatio',ratio, ...
    'SettlingTime2Percent',settle,'VoltageRMS',sqrt(mean(voltageError(:).^2)), ...
    'MaxVoltageError',max(abs(voltageError),[],'all'), ...
    'SharingRMS',sqrt(mean(sharing(:).^2)), ...
    'TerminalSharingError',max(abs(sharing(end,:))), ...
    'PeakCommandDeviation',max(abs(U-net.u_s'),[],'all'), ...
    'CommandRMSDeviation',sqrt(mean((U-net.u_s').^2,'all')), ...
    'StateEnergy',trapz(t,sum(scaled.^2,2)), ...
    'CommandExceedances',nnz(abs(U)>2*Vr), ...
    'LineCurrentExceedances',nnz(abs(lineCurrent)>2*Ir));
tail=t>=max(0,options.FinalTime-0.02);
metrics.TailVoltageRMS=sqrt(mean(voltageError(tail,:).^2,'all'));
metrics.TailSharingRMS=sqrt(mean(sharing(tail,:).^2,'all'));
forcing=primitive*net.BBar'/48*1e-3;
metrics.DisturbanceEnergy=trapz(t,sum(forcing.^2,2));
metrics.EmpiricalGain=NaN;
if metrics.DisturbanceEnergy>0 && normError(1)<1e-12 && ~any(kicks(:))
    metrics.EmpiricalGain=sqrt(metrics.StateEnergy/metrics.DisturbanceEnergy);
end
result=struct('Time',t,'State',X,'StateError',error,'Command',U, ...
    'LineCurrent',lineCurrent,'PrimitiveDisturbance',primitive,'StateJumps',kicks, ...
    'SteadyState',net.ss,'CommandGain',Kx,'ControllerState',controllerState, ...
    'ClosedLoopMatrix',A,'Metrics',metrics,'Options',options, ...
    'Normalization',48,'LinearModel',true,'LimitersEnabled',false);
result.Certificate=options.Certificate;
result.Accepted=all(isfinite(X(:))) && metrics.MaxRealPole<0;
end

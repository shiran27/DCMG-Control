function controller = designConventionalDroop(net,options)
%DESIGNCONVENTIONALDROOP Conventional decentralized V-I droop baseline.
%   controller=designConventionalDroop(net) implements
%
%       V_ti = V_0i - R_di I_ti,
%       R_di = delta_V V_rated,i / I_rated,i.
%
%   The default delta_V=0.05 gives a 5% voltage reduction at rated current.
%   V_0i is chosen so that the controller shares the benchmark equilibrium
%   with every proposed method. In error coordinates this is simply
%   delta V_ti = -R_di delta I_ti. The baseline therefore has no remote
%   measurements, communication links, optimizer, or claimed certificate.
%
%   Literature contract:
%     J. M. Guerrero et al., IEEE TIE 58(1), 2011,
%       doi:10.1109/TIE.2010.2066534.
%     X. Lu et al., IEEE TPEL 29(4), 2014,
%       doi:10.1109/TPEL.2013.2266419.

arguments
    net
    options.VoltageDroopFraction (1,1) double {mustBePositive,mustBeFinite} = 0.05
end
assert(options.VoltageDroopFraction<1,'DCMG:DroopFraction', ...
    'VoltageDroopFraction must be strictly between zero and one.');
assert(~isempty(net.x_s) && ~isempty(net.u_s),'DCMG:DroopEquilibrium', ...
    'Solve and verify the common equilibrium before constructing droop control.');

N=net.N;
Vr=[net.DGs.Vrated]';
Ir=[net.DGs.Irated]';
Iss=net.x_s(2:2:end);
Rdroop=options.VoltageDroopFraction*Vr./Ir;
K=zeros(N,2*N);
for i=1:N
    K(i,2*i)=-Rdroop(i);
end

% V0 is the no-load intercept in the absolute droop law. This bias makes
% V0-Rd*Iss exactly equal to the already verified equilibrium command.
V0=net.u_s+Rdroop.*Iss;
equilibriumResidual=norm(V0-Rdroop.*Iss-net.u_s,inf);
Aphysical=net.A+net.E*net.YBar*net.D';
poles=eig(Aphysical+net.B*K);

controller=struct();
controller.Name=sprintf('Conventional droop (%g%%)', ...
    100*options.VoltageDroopFraction);
controller.Type='literature baseline';
controller.Formula='V_ti = V_0i - R_di I_ti';
controller.VoltageDroopFraction=options.VoltageDroopFraction;
controller.DroopResistance=Rdroop;
controller.NoLoadVoltage=V0;
controller.CommandGain=K;
controller.K=K;
controller.Local={};
controller.DirectedEdges=0;
controller.Adjacency=zeros(N);
controller.Gamma=NaN;
controller.Problem=0;
controller.problem=0;
controller.info='Fixed conventional V-I droop; no synthesis solver used';
controller.EquilibriumCommandResidual=equilibriumResidual;
controller.Poles=poles;
controller.MaxRealPole=max(real(poles));
controller.Accepted=equilibriumResidual<1e-12 && controller.MaxRealPole<0;
controller.CertificateClaim='None; empirical literature baseline only';
controller.CitationKeys={'Guerrero2011HierarchicalDroop', ...
    'Lu2014ImprovedDroop'};
assert(controller.Accepted,'DCMG:UnstableDroopBaseline', ...
    'The fixed droop baseline is not stable on the selected benchmark.');
end

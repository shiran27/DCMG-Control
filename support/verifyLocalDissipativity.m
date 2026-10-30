function diagnostics = verifyLocalDissipativity(dg,certificate)
%VERIFYLOCALDISSIPATIVITY Truth-only check of the original storage inequality.
o=certificate.Options; H=certificate.Storage;
A=o.TimeScale*dg.A;
B=o.TimeScale*dg.BBar*o.InputScale/o.StateScale;
Acl=A+B*certificate.FullGain;
R=[-certificate.rho*eye(2)-Acl'*H-H*Acl,0.5*eye(2)-H; ...
   0.5*eye(2)-H,certificate.a*eye(2)];
diagnostics=certificateDiagnostics(R);
diagnostics.Matrix=R;
diagnostics.MaxRealPole=max(real(eig(Acl)))/o.TimeScale;
diagnostics.Accepted=diagnostics.Accepted && diagnostics.MaxRealPole<0;
end

function out = verifyGlobalStorage(net,cert)
%VERIFYGLOBALSTORAGE Original aggregate CT bounded-real storage inequality.
H=[];
for i=1:net.N
    H=blkdiag(H,cert.Weights(i)*cert.Local{i}.Storage);
end
A=cert.Options.TimeScale*(net.A+net.E*net.YBar*net.D'+net.B*cert.CommandGain);
n=2*net.N;
R=[-A'*H-H*A-eye(n),-H;-H,cert.GammaSquared*eye(n)];
out=certificateDiagnostics(R); out.Matrix=R; out.Storage=H;
out.MaxRealPole=max(real(eig(A)))/cert.Options.TimeScale;
out.Accepted=out.Accepted && out.MaxRealPole<0;
end

function out = certificateDiagnostics(M)
%CERTIFICATEDIAGNOSTICS Scale-aware, eigenvalue-based certificate audit.
assert(all(isfinite(M(:))),'DCMG:NonfiniteCertificate','Nonfinite certificate.');
S=(M+M')/2;
out.MinEigenvalue=min(eig(S));
out.Scale=max(1,norm(S,2));
out.RelativeMargin=out.MinEigenvalue/out.Scale;
out.SymmetryResidual=norm(M-M','fro')/max(1,norm(M,'fro'));
out.Accepted=out.MinEigenvalue>=-1e-6 && out.RelativeMargin>=-1e-8 ...
    && out.SymmetryResidual<=1e-9;
end

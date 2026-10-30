function restoreControllerState(net,state)
%RESTORECONTROLLERSTATE Preserve the installed controller after a failed solve.
net.K=state.K; net.commAdj=state.Adjacency; net.DesignCertificate=state.Certificate;
for i=1:net.N
    dg=net.DGs(i); s=state.Local{i};
    dg.K=s.K; dg.nu=s.nu; dg.rho=s.rho;
    dg.DesignCertificate=s.Certificate; dg.DesignOptions=s.Options;
end
end

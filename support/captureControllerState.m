function state = captureControllerState(net)
%CAPTURECONTROLLERSTATE Snapshot gains before a transactional co-design.
state=struct('K',net.K,'Adjacency',net.commAdj,'Certificate',net.DesignCertificate, ...
    'Local',{cell(1,net.N)});
for i=1:net.N
    dg=net.DGs(i);
    state.Local{i}=struct('K',dg.K,'nu',dg.nu,'rho',dg.rho, ...
        'Certificate',dg.DesignCertificate,'Options',dg.DesignOptions);
end
end

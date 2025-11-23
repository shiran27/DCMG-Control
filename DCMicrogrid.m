classdef DCMicrogrid < handle
    % Holds arrays of DGs and TransmissionLines; simulates and draws
    properties
        DGs
        Lines
        N  (1,1) double
        M  (1,1) double
        G                 % NxN conductance Laplacian (diag=sum g_ij, offdiag=-g_ij)

        A           % (2N x 2N)
        BBar        % (2N x 2N)
        E           % (2N x  N)
        B           % (2N x  N)
        D           % (2N x 2N)   % picks vC per node
        DBar           % (2N x 2N)   % picks iL per node
        
       
        u_s         % (N x 1)
        wBar        % (N x 1)

        % --- Admittance matrices from line resistances ---
        Y       % (N x N)  off-diagonal admittances, diagonal = 0
        YBar    % (N x N)  Laplacian of admittances (diag=sum, offdiag = -Y)   

        x_s % steady-state
        ss
        % K_L % Strict local controller gains
        K           % (N  x 2N)
        commAdj % Adjacencuy matrix of the com topology


        % Integral action
        z           % (2N x 1) integrator state
        K_I         % (N x 4N) integral gain

        % Data 
        QBar_w
        E_perm

        % Kc          % continuous global controller gain 
        % Kd          % discrete-time global gain 
        Ts     % controller sampling period used in data-driven design
        u_G_Hold
        I_Hold

    end

    methods
        
        function obj = DCMicrogrid(DG_array, Line_array)

            obj.DGs   = DG_array(:).';    % row
            obj.Lines = Line_array(:).';

            obj.N = numel(obj.DGs); 
            obj.M = numel(obj.Lines);

            obj.buildConductance();

            for k = 1:obj.N
                obj.DGs(k).updateModel();
            end

        end

        % buildConductance: keep your Y, YBar, but don't push gsum into DGs anymore
        function buildConductance(obj)
            N = obj.N;  Y = zeros(N); YBar = zeros(N);
            for e = 1:obj.M
                i = obj.Lines(e).i; 
                j = obj.Lines(e).j; 
                g = obj.Lines(e).g;
                Y(i,j) = g; 
                Y(j,i) = g;
                
                YBar(i,j) = -g; 
                YBar(j,i) = -g;

                YBar(i,i) = YBar(i,i)+g; 
                YBar(j,j) = YBar(j,j)+g;
            end
            obj.Y = Y; 
            obj.YBar = YBar; 

            for i = 1:1:N
                obj.DGs(i).YBar = YBar(i,i);
            end
        end


        function X = getStateVector(obj)
            % Stack [vC1 iL1 vC2 iL2 ...]'
            X = zeros(2*obj.N,1);
            for k = 1:obj.N
                X(2*k-1:2*k) = obj.DGs(k).getState();
            end
        end

        function setStateVector(obj, X)
            for k = 1:obj.N
                obj.DGs(k).setState(X(2*k-1:2*k));
            end
        end

        
        function [t, X, U, Xs, Us] = simulate(obj, tspan, x0, useData, linkFiltThresh)

            ode_opts = odeset('RelTol',1e-6,'AbsTol',1e-6);
            
            dt = 1e-5;
            t_0 = tspan(1);
            t_f = tspan(end);

            t_1 = 0.05*t_f; % When the data-driven method is triggered 
            t_2 = 0.2*t_f; % When a high disturbance enters
            t_3 = 0.4*t_f; % When a state disturbance enters
            t_4 = 0.6*t_f; % When a parameter disturbance enters
            t_5 = 0.65*t_f; % When data-driven method reacts
            t_6 = 0.8*t_f; % When high disturbance and a state disturbance enters together
            
            t = [];
            X = [];
            U = [];
            Xs = [];
            Us = [];


            % Get ready for Stage 1: Only for Data-Driven Co-Designed DRC
            if useData 
                % A noisy version of Model-based design of a global stabilizing controller 
                obj.design_MB_GSC(linkFiltThresh,useData);
            end
            
            disp('Starting Stage 1...')
            % Stage 1
            x01 = x0;
            f = @(t,x) obj.dynamics(t, x, false);
            [t1, X1] = ode45(f, [t_0:dt:t_1], x01, ode_opts);
            [U1, ~, Utilde1, Wtilde1, Us1, Xs1] = obj.computeUWTrajectories(t1, X1, false);
            t  = [t;  t1];
            X  = [X;  X1];
            U  = [U;  U1];
            Xs = [Xs; Xs1];
            Us = [Us; Us1];
            disp('Finished Stage 1, getting ready for Stage 2...')
            
            
            % Get ready for Stage 2: Only for Data-Driven Co-Designed DRC
            if useData
                obj.loadDataMatrices(t1, X1, Xs1, Utilde1, Wtilde1);
                [AdjMat, KMat, out] = obj.codesign_DD_DRC(linkFiltThresh);
            end

            disp('Starting Stage 2...')
            % Stage 2
            x02 = X(end, :).';
            if useData
                [t2, X2, U2, Xs2, Us2] = obj.simulateStageZOH(t_1, t_2, x02, dt);
            else
                f = @(t,x) obj.dynamics(t, x, useData);
                [t2, X2] = ode45(f, [t_1:dt:t_2], x02, ode_opts);
                [U2, ~, ~, ~, Us2, Xs2] = obj.computeUWTrajectories(t2, X2, useData);
            end
            t = [t; t2(2:end)];
            X = [X; X2(2:end,:)];
            U = [U; U2(2:end,:)];
            Xs = [Xs; Xs2(2:end,:)];
            Us = [Us; Us2(2:end,:)];
            disp('Finished Stage 2, getting ready for Stage 3...')
            
            
            % Get ready for Stage 3: instantaneous State error
            mag = 0.2;
            for i = 1:obj.N
                if rand(1)<0.5
                    X(end, 2*i-1:2*i) = X(end, 2*i-1:2*i) + ...
                        X(end, 2*i-1:2*i).*(mag*(2*rand(1,2)-2));
                end
            end
            
            disp('Starting Stage 3...')
            % Stage 3
            x03 = X(end, :).';
            if useData
                [t3, X3, U3, Xs3, Us3] = obj.simulateStageZOH(t_2, t_3, x03, dt);
            else
                f = @(t,x) obj.dynamics(t, x, useData);
                [t3, X3] = ode45(f, [t_2:dt:t_3], x03, ode_opts);
                [U3, ~, ~, ~, Us3, Xs3] = obj.computeUWTrajectories(t3, X3, useData);
            end
            t = [t; t3(2:end)];
            X = [X; X3(2:end,:)];
            U = [U; U3(2:end,:)];
            Xs = [Xs; Xs3(2:end,:)];
            Us = [Us; Us3(2:end,:)];
            disp('Finished Stage 3, getting ready for Stage 4...')
            
            
            % Get ready for Stage 4: Heightened disturbance period 
            mag = 100;
            for i = 1:obj.N
                tIndices = (obj.DGs(i).noise.t > t_3 & obj.DGs(i).noise.t < (0.8*t_3 + 0.2*t_4));
                obj.DGs(i).noise.w(tIndices,:) = mag*obj.DGs(i).noise.w(tIndices,:);
            end

            disp('Starting Stage 4...')
            % Stage 4
            x04 = X(end, :).';
            if useData
                [t4, X4, U4, Xs4, Us4] = obj.simulateStageZOH(t_3, t_4, x04, dt);
            else
                f = @(t,x) obj.dynamics(t, x, useData);
                [t4, X4] = ode45(f, [t_3:dt:t_4], x04, ode_opts);
                [U4, ~, ~, ~, Us4, Xs4] = obj.computeUWTrajectories(t4, X4, useData);
            end
            t = [t; t4(2:end)];
            X = [X; X4(2:end,:)];
            U = [U; U4(2:end,:)];
            Xs = [Xs; Xs4(2:end,:)];
            Us = [Us; Us4(2:end,:)];
            disp('Finished Stage 4, getting ready for Stage 5...')
            
            
            % Get ready for Stage 5: Permenent model parameter variation
            mag = 0.2;
            % obj.A
            for i = 1:obj.N
                obj.DGs(i).A    = obj.DGs(i).A    + obj.DGs(i).A.*(mag*(2*rand(size(obj.DGs(i).A))-1));
                % obj.DGs(i).B    = obj.DGs(i).B    + obj.DGs(i).B.*(mag*(2*rand(size(obj.DGs(i).B))-1));
                % obj.DGs(i).B(1) = 0;
                % obj.DGs(i).E    = obj.DGs(i).E    + obj.DGs(i).E.*(mag*(2*rand(size(obj.DGs(i).E))-1));
                % obj.DGs(i).E(2) = 0;
                obj.DGs(i).Ibar = obj.DGs(i).Ibar + obj.DGs(i).Ibar.*(mag*(2*rand(size(obj.DGs(i).Ibar))-1));
            end
           
            obj.buildSystemMatrices();
            % obj.A
            % obj.u_s
            % res1 = obj.A*obj.x_s + obj.E*obj.ss.I_ss + obj.B*obj.u_s + obj.E*obj.wBar;
            % fprintf('||steady-state residual1||_2 = %.3e\n', norm(res1));
            % res2 = obj.ss.I_ss - obj.YBar* obj.D'*obj.x_s;
            % fprintf('||steady-state residual2||_2 = %.3e\n', norm(res2));
            obj.solveSteadyState();
            % obj.u_s
            % res1 = obj.A*obj.x_s + obj.E*obj.ss.I_ss + obj.B*obj.u_s + obj.E*obj.wBar;
            % fprintf('||steady-state residual||_2 = %.3e\n', norm(res1));
            % res2 = obj.ss.I_ss - obj.YBar* obj.D'*obj.x_s;
            % fprintf('||steady-state residual2||_2 = %.3e\n', norm(res2));
            
            disp('Starting Stage 5...')
            % Stage 5
            x05 = X(end, :).';
            f = @(t,x) obj.dynamics(t, x, false);
            [t5, X5] = ode45(f, [t_4:dt:t_5], x05, ode_opts);
            [U5, ~, Utilde5, Wtilde5, Us5, Xs5] = obj.computeUWTrajectories(t5, X5, false);
            t = [t; t5(2:end)];
            X = [X; X5(2:end,:)];
            U = [U; U5(2:end,:)];
            Xs = [Xs; Xs5(2:end,:)];
            Us = [Us; Us5(2:end,:)];
            disp('Finished Stage 5, getting ready for Stage 6...')
            
            
            % Get ready for Stage 6: Only Data-Driven Co-Designed DRC
            if useData
                obj.loadDataMatrices(t5, X5, Xs5, Utilde5, Wtilde5);
                [AdjMat, KMat, out] = obj.codesign_DD_DRC(linkFiltThresh);
            end

            disp('Starting Stage 6...')
            % Stage 6 
            x06 = X(end, :).';
            if useData
                [t6, X6, U6, Xs6, Us6] = obj.simulateStageZOH(t_5, t_6, x06, dt);
            else
                f = @(t,x) obj.dynamics(t, x, useData);
                [t6, X6] = ode45(f, [t_5:dt:t_6], x06, ode_opts);
                [U6, ~, ~, ~, Us6, Xs6] = obj.computeUWTrajectories(t6, X6, useData);
            end
            t = [t; t6(2:end)];
            X = [X; X6(2:end,:)];
            U = [U; U6(2:end,:)];
            Xs = [Xs; Xs6(2:end,:)];
            Us = [Us; Us6(2:end,:)];

            disp('Finished Stage 6, getting ready for Stage 7.')

            % Get ready for Stage 7: In this final stage, for the parameter
            % changed system, additional state and disturbance noise impacts.
            % State noise:
            mag = 0.2;
            for i = 1:obj.N
                if rand(1)<0.5
                    X(end, 2*i-1:2*i) = X(end, 2*i-1:2*i) + ...
                        X(end, 2*i-1:2*i).*(mag*(2*rand(1,2)-2));
                end
            end
            % Disturbance amplification
            mag = 100;
            for i = 1:obj.N
                tIndices = (obj.DGs(i).noise.t > t_6 & obj.DGs(i).noise.t < (0.8*t_6 + 0.2*t_f));
                obj.DGs(i).noise.w(tIndices,:) = mag*obj.DGs(i).noise.w(tIndices,:);
            end

            disp('Starting Stage 7...')
            % Stage 7
            x07 = X(end, :).';
            if useData
                [t7, X7, U7, Xs7, Us7] = obj.simulateStageZOH(t_6, t_f, x07, dt);
            else
                f = @(t,x) obj.dynamics(t, x, useData);
                [t7, X7] = ode45(f, [t_6:dt:t_f], x07, ode_opts);
                [U7, ~, ~, ~, Us7, Xs7] = obj.computeUWTrajectories(t7, X7, useData);
            end
            t = [t; t7(2:end)];
            X = [X; X7(2:end,:)];
            U = [U; U7(2:end,:)];
            Xs = [Xs; Xs7(2:end,:)];
            Us = [Us; Us7(2:end,:)];
            disp('Finished Stage 7, i.e., the Final.')
            

        end




        function [t_seg, X_seg, U_seg, Xs_seg, Us_seg] = ...
            simulateStageZOH(obj, t_start, t_end, x0, dt)
        
            % One-stage simulation with zero-order-held global controller:
            %   u_Gk = Kd * (x_k - x_s)
            %   xdot = dynamicsZOH(x, u_Gk)
            % Sample period = obj.Ts_ctrl
        
            ode_opts = odeset('RelTol',1e-6,'AbsTol',1e-6);
            Ts       = obj.Ts;
        
            t_seg  = [];
            X_seg  = [];
            U_seg  = [];
            Xs_seg = [];
            Us_seg = [];
        
            t0 = t_start;
            Xk = x0(:)';
        
            
            while t0 < t_end
                
                tk1 = min(t0 + Ts, t_end);
                
            
                % Error);
                XTilde = Xk' - obj.x_s;
            
                % Feedback
                u_G = obj.K*XTilde;

                obj.u_G_Hold = u_G;
                

                for i = 1:obj.N
                    
                    x_i = Xk(1,2*i-1:2*i);  % from setStateVector
                    xTilde_i = x_i' - obj.DGs(i).x_s;
                
                    obj.DGs(i).u_L_Hold = obj.DGs(i).K * xTilde_i;  % 2×1 local held input
                end

                % 2) integrate continuous plant with u_G(t) ≡ uGk on [t0, tk1]
                f = @(t_local, x_local) obj.dynamics(t_local, x_local, true);
                if length(t0:dt:tk1)<=1
                    [t, X] = ode45(f, [t0, tk1], Xk, ode_opts);
                else
                    [t, X] = ode45(f, t0:dt:tk1, Xk, ode_opts);
                end
        
                % 3) compute U, Xs, Us for this segment using ZOH uGk
                [U, ~, ~, ~, Us, Xs] = obj.computeUWTrajectories(t, X, true);
        
                % 4) append, avoid duplicate t
                if isempty(t_seg)
                    t_seg  = t;
                    X_seg  = X;
                    U_seg  = U;
                    Xs_seg = Xs;
                    Us_seg = Us;
                else
                    t_seg  = [t_seg;  t(2:end)];
                    X_seg  = [X_seg;  X(2:end,:)];
                    U_seg  = [U_seg;  U(2:end,:)];
                    Xs_seg = [Xs_seg; Xs(2:end,:)];
                    Us_seg = [Us_seg; Us(2:end,:)];
                end
        
                t0 = tk1;
                Xk = X(end,:);

            end
        end

        function out = loadDataMatrices(obj, t, X, Xs, Utilde, Wtilde)

            meanAbsDistVals = [];
            meanAbsDiscEVals = [];
            QBar_w = [];
            E1 = []; E2 = []; E3 = [];
            N = obj.N;

            % Data sampling
            Ts = 1e-4;
            obj.Ts = Ts;
            tVals = t(1):Ts:t(end);


            for k = 1:N
            
                obj.DGs(k).Ts = Ts;

                % Loading xTIlde and uTilde Data
                xTilde = X(:, 2*k-1:2*k) - Xs(2*k-1:2*k,1)';
                uTilde = Utilde(:, 2*k-1:2*k);
                wTilde = Wtilde(:, 2*k-1:2*k);
            
                xTildeSampled = interp1(t, xTilde, tVals, 'linear')';
                uTildeSampled = interp1(t, uTilde, tVals, 'previous', 'extrap')';
                wTildeSampled = Ts*interp1(t, wTilde, tVals, 'previous', 'extrap')';
            
                % Data loading
                obj.DGs(k).xTilde = xTildeSampled(:,1:end-1);
                obj.DGs(k).yTilde = obj.DGs(k).xTilde;
                obj.DGs(k).xTildeBar = xTildeSampled(:,2:end);
                obj.DGs(k).uTilde = uTildeSampled(:,1:end-1);
                obj.DGs(k).wTilde = wTildeSampled(:,1:end-1);
            
                % Investigating discretization Error and Disturbance impact
                meanAbsDistVal = mean(abs(obj.DGs(k).wTilde')); meanAbsDistVals = [meanAbsDistVals; meanAbsDistVal];
                discError = obj.DGs(k).xTildeBar - ((eye(2)+Ts*obj.DGs(k).A)*obj.DGs(k).xTilde ...
                                        + Ts*obj.DGs(k).BBar*obj.DGs(k).uTilde);
                meanAbsDiscEVal = mean(abs(discError')); meanAbsDiscEVals = [meanAbsDiscEVals; meanAbsDiscEVal];
                
                % Loading the discretization error as the disturbance 
                obj.DGs(k).wTilde = discError; %%%% check
                
                % Finding an upper bound for disturbance: Q_w 
                wwT = obj.DGs(k).wTilde*obj.DGs(k).wTilde';
                lambda = max(eig(wwT))*1.001; %%%% chech
                % eig(lambda*eye(2)-wwT);
                Q_ww = -eye(length(tVals)-1);
                Q_I = lambda*eye(2);
                Q_z = zeros(2,length(tVals)-1);
                Q_w = [Q_I, Q_z; Q_z', Q_ww];
                obj.DGs(k).Q_w = Q_w;
            
                % Loading QBar_w at DGs
                L = [eye(2), obj.DGs(k).xTildeBar; 
                    zeros(2,2), -obj.DGs(k).xTilde; 
                    zeros(2,2), -obj.DGs(k).uTilde];
                QBar_wk = L*Q_w*L';
                QBar_wk = 0.5*(QBar_wk + QBar_wk');
                obj.DGs(k).QBar_w = QBar_wk;
            
                % Loading for QBar_w at the DCMG
                QBar_w = blkdiag(QBar_w, QBar_wk);

                % Permutation matriz components
                E1 = blkdiag(E1, [eye(2), zeros(2), zeros(2)]);
                E2 = blkdiag(E2, [zeros(2), eye(2), zeros(2)]);
                E3 = blkdiag(E3, [zeros(2), zeros(2), eye(2)]);
            end
            obj.QBar_w = QBar_w;
            obj.E_perm = [E1; E2; E3];
            
            meanAbsDist = mean(meanAbsDistVals)
            meanAbsDiscErr = mean(meanAbsDiscEVals)

        end


        function [U, Uc, Utilde, Wtilde, Us, Xs] = computeUWTrajectories(obj, t, X, useData)
            % U: (nt x 2N)      -> [[I_i = I_iS + I_iG, V_ti = u_iT = u_iS + u_iL + u_iG]_i=1:N]
            % Uc: (nt x N x 5 x 1)      -> [[I_iS, I_iG, u_iS, u_iL, u_iG]_i=1:N]
            % Utilde: (nt x 2N)      -> [[I_iG, u_iL + u_iG]_i=1:N]
            % Wtilde: (nt x 2N)      -> [[(BBar_i*w_i)']_i=1:N]
        
            N = obj.N;
            nt = length(t);
        
            U = zeros(nt, 2*N);
            Uc = zeros(nt, N, 5, 1);
            UTilde = zeros(nt, 2*N);
            WTilde = zeros(nt, 2*N);
            Us = zeros(nt, 2*N);
            Xs = ones(nt,1)*obj.x_s';

            V_s = obj.D' * obj.x_s;
            I_S = obj.YBar * V_s;

            for k = 1:nt

                Xk = X(k, :).';           % current global state
                obj.setStateVector(Xk);   % pushes xk into each DG.x internally
        
                % Network quantities
                Vk = obj.D' * Xk;
                Ik = obj.YBar * Vk;
                
                xTildek = Xk - obj.x_s;
                if useData
                    u_Gk = obj.u_G_Hold;        % global/secondary input (size N x 1)
                else
                    u_Gk = obj.K * xTildek;        % global/secondary input (size N x 1)
                end

                % DG-by-DG local and total inputs
                for i = 1:N

                    DG   = obj.DGs(i);
        
                    % DG.x was set by setStateVector
                    x_i  = DG.x;                       % local state of DG i

                    ITilde_ki = Ik(i) - I_S(i);

                    u_Si = DG.u_s;
                    if useData
                        u_Li = DG.u_L_Hold;
                    else
                        u_Li = DG.K * (x_i - DG.x_s);      % local input
                    end
                    u_Gi = u_Gk(i);
                    u_i  = u_Si + u_Li + u_Gi;     % total input
        
                    uTilde_ki = u_i - u_Si;
                    
                    tVal = t(k);
                    w_i = interp1(DG.noise.t, DG.noise.w, tVal, 'previous', 'extrap')';
                    Wtilde_ki = (DG.BBar*w_i);

                    U(k,2*i-1:2*i) = [Ik(i), u_i];
                    Uc(k,i,:,:) = [I_S(i); ITilde_ki; u_Si; u_Li; u_Gi];
                    Utilde(k,2*i-1:2*i) = [ITilde_ki, uTilde_ki];
                    Wtilde(k,2*i-1:2*i) = Wtilde_ki; 
                    Us(k,2*i-1:2*i) = [I_S(i), u_Si];
                end
            end

        end

        



        function dX = dynamics(obj, t, X, useData)
            N  = obj.N;
            
            obj.setStateVector(X);
            % Electrical injections
            V = obj.D'*X;
            Iline = obj.YBar * V;
        
            if useData
                u_G = obj.u_G_Hold;
            else
                % Error);
                XTilde = X - obj.x_s;
            
                % Feedback
                u_G = obj.K*XTilde;
            end
        
            % Plant derivatives DG-by-DG with controller-provided u_i
            dX = zeros(2*N,1);
            for i = 1:N
                dxi = obj.DGs(i).dynamics(t, Iline(i), u_G(i), useData);
                dX(2*i-1:2*i) = dxi;
            end

        end


       


        function obj = setupNoise(obj, tspan, dt_noise, sigma)
            t0 = tspan(1);
            tf = tspan(end);
        
            t_noise = (t0:dt_noise:tf).';
            % One scalar disturbance per DG; adjust size as needed

            for i = 1:obj.N
                w = randn(length(t_noise), 2)*sigma;
                obj.DGs(i).noise.t = t_noise;
                obj.DGs(i).noise.w = w;
            end
        end


        


        function draw(obj, ax)

            if nargin < 2 || isempty(ax), ax = gca; end
            cla(ax); hold(ax, 'on'); axis(ax, 'equal');

            obj.drawComm(ax);

            % Lines first
            for e = 1:obj.M
                i = obj.Lines(e).i; j = obj.Lines(e).j;
                obj.Lines(e).draw(ax, obj.DGs(i).pos, obj.DGs(j).pos);
            end
            % DGs on top
            for k = 1:obj.N
                obj.DGs(k).draw(ax);
            end
            grid(ax, 'on'); xlabel(ax,'x'); ylabel(ax,'y');
            title(ax, 'DC Microgrid Topology');
            set(ax,'Visible','off') 

        end

        function drawComm(obj, ax, varargin)
            % Draw dashed lines for nonzero 1x2 blocks in K (thresholded).
            if nargin < 2 || isempty(ax), ax = gca; end
            if isempty(obj.K), return; end
            N = obj.N; K = obj.K; tol = 1e-12;
            hold(ax, 'on');
            for i = 1:N
                for j = 1:N
                    blk = K(i, 2*j-1:2*j);
                    if norm(blk,2) > tol && i ~= j
                         
                        pi = obj.DGs(i).pos; 
                        pj = obj.DGs(j).pos;
                        
                        % Perpendicular midpoint bulge
                        dline = (pj - pi)/norm(pj - pi);
                        dperp = [dline(2) - dline(1)];
                        
                        bulge = 0.4;
                        pm = (0.5*pi + 0.5*pj) + bulge * dperp;
                        
                        % ----- Quadratic Bézier -----
                        num_points = 200;
                        t = linspace(0, 1, num_points)';
                        B = (1 - t).^2 .* pi + 2*(1 - t).*t .* pm + t.^2 .* pj;
                        
                        x_curve = B(:,1);
                        y_curve = B(:,2);
                        
                        % Plot curve
                        plot(ax, x_curve, y_curve, '--', ...
                            'Color', [0.1 0.5 1.0], 'LineWidth', 0.8);
                        hold(ax, 'on');
                        
                        % ----- Arrowhead at MIDDLE of curve -----
                        tm = 0.5;   % midpoint along Bezier
                        
                        % Position on curve at t = 0.5
                        P = (1 - tm)^2 * pi + 2*(1 - tm)*tm * pm + tm^2 * pj;
                        
                        % Tangent direction on quadratic Bezier:
                        
                        % Arrowhead size (scale with edge length)
                        arrow_len = 0.1;
                        arrow_wid = 0.05;
                        
                        % Triangle vertices (tip at P)
                        tip     = P;
                        leftpt  = P - arrow_len*dline + arrow_wid*dperp;
                        rightpt = P - arrow_len*dline - arrow_wid*dperp;
                        
                        % Draw triangle arrowhead
                        patch(ax, [tip(1) leftpt(1) rightpt(1)], ...
                                  [tip(2) leftpt(2) rightpt(2)], ...
                                  [0.1 0.5 1.0], ...
                                  'EdgeColor', 'none');


                        
                    end
                end
            end
        end



        function buildSystemMatrices(obj)
            N = obj.N;

            Ablk   = cell(1,N);
            Eblk   = cell(1,N);
            Bblk   = cell(1,N);
            BBarblk   = cell(1,N);
            Dblk   = cell(1,N);
            DBarblk   = cell(1,N);

            obj.wBar = zeros(N,1);

            for i = 1:N
                Ai = obj.DGs(i).A;
                Ei = obj.DGs(i).E;   % (2x1)
                Bi = obj.DGs(i).B;   % (2x1)

                Ablk{i}   = Ai;
                Eblk{i}   = Ei;
                Bblk{i}   = Bi;
                BBarblk{i}   = [Ei, Bi];
                Dblk{i}   = [1; 0];  % vC selector
                DBarblk{i}   = [0; 1];  % iL selector

                obj.wBar(i) = obj.DGs(i).Ibar;
            end

            obj.A    = blkdiag(Ablk{:});
            obj.E    = blkdiag(Eblk{:});
            obj.B    = blkdiag(Bblk{:});
            obj.BBar = blkdiag(BBarblk{:});
            obj.D    = blkdiag(Dblk{:});
            obj.DBar    = blkdiag(DBarblk{:});

        end


        function ss = solveSteadyState(obj, opts)
            % Solve steady state including epsilon matrices:
            %   A x_ss + E I_ss + B u_ss + E wBar = 0
            %   V_ss = diag(epsV) * V_rated
            %   I_ss = epsI * I_rated             % equal current sharing (scalar epsI)
            %
            % Minimizes  wV*||epsV-1||^2 + wI*(epsI-1)^2 + wu*||u-us_prev||^2
            % with optional bounds on epsV and epsI.
            %
            % Requires YALMIP + MOSEK.
            
            % arguments
            %     obj
            %     opts.wV (1,1) double = 1.0
            %     opts.wI (1,1) double = 1.0
            %     opts.wu (1,1) double = 1e-3
            %     opts.deltaV (1,1) double = 0.10   % |epsV_i-1| <= deltaV
            %     opts.deltaI (1,1) double = 0.10   % |epsI-1| <= deltaI
            %     opts.set_state (1,1) logical = true
            %     opts.set_us    (1,1) logical = true
            % end
        
            
            N    = obj.N;    
            nX   = 2*N;
        
            % Gather ratings from DGs
            Vr = zeros(N,1); 
            Ir = zeros(N,1);
            for i = 1:N
                Vr(i) = obj.DGs(i).Vrated;
                Ir(i) = obj.DGs(i).Irated;
            end
        
            % Plant
            Am = obj.A;
            Bm = obj.B;   
            Em = obj.E;

            wBar = obj.wBar; 
            D = obj.D; 
            DBar = obj.DBar;
            YBar = obj.YBar;
        
            % Decision vars
            x = sdpvar(nX,1);

            u = sdpvar(N,1); % Second component of u_E
            I = sdpvar(N,1); % First compoenent of u_E
            epsV = sdpvar(N,1);           % diagonal entries of Sigma_V
            epsI = sdpvar(1,1);           % equal-sharing scalar
        
            % Selections
            V = D'*x; 
            I_t = DBar'*x; 
            
            % Relations to ratings
            constr = [];
            constr = [constr, V == diag(epsV) * Vr];   % voltage regulation
            constr = [constr, I_t == epsI * Ir];       % equal current sharing
        
            % Physics at steady state
            constr = [constr, Am*x + [Em Bm]*[I; u] + Em*wBar == 0];
            constr = [constr, I == YBar*V];
        
            % Epsilon bounds
            constr = [constr, epsV >= 0.8, epsV <= 1.2];
            constr = [constr, epsI >= 0, epsI <= 0.98];
        
            % Mild input regularization
            % us_prev = mats.u_s; if isempty(us_prev), us_prev = zeros(N,1); end
        
            objFun = norm(epsV - 1, 2)^2 + epsI;
        
            ops = sdpsettings('solver','mosek','verbose',0);
            info = optimize(constr, objFun, ops);
        
            ss = struct();
            ss.problem = info.problem;
            ss.info    = yalmiperror(info.problem);
            ss.x_ss    = value(x);
            ss.u_s     = value(u);
            ss.epsV    = value(epsV);
            ss.epsI    = value(epsI);
            ss.V_ss    = value(V);
            ss.I_tss    = value(I_t);
            ss.I_ss    = value(I);
        
            if info.problem ~= 0
                warning('Steady-state design failed: %s', ss.info);
                return;
            end
            
            disp('Steady-state design SUCCESS!');

            % Push solution to plant (optional)
            
            obj.u_s = ss.u_s;
            obj.x_s = ss.x_ss;
            obj.z   = zeros(obj.N,1);   % start integrators at 0 error
            for k = 1:N
                obj.DGs(k).u_s = ss.u_s(k);
                obj.DGs(k).x_s = ss.x_ss(2*k-1:2*k,1);
                obj.DGs(k).I_s = ss.I_ss(k);
            end
            

            fprintf('SS current usage (fraction) at each DG: %s \n', ss.epsI);
            fprintf('SS voltage levels (Volt) at DGs: \n');    
            ss.epsV'

            obj.ss = ss;

            % % % Verify physics at solution (should be ~zero)
            % Am = obj.A; Em = obj.E; Bm = obj.B; YB = obj.YBar;
            % xss = ss.x_ss; uss = ss.u_s; Vss = ss.V_ss; Iliness = YB*Vss;
            % res = Am*xss + Em*Iliness + Bm*uss + Em*obj.wBar;
            % fprintf('||steady-state residual||_2 = %.3e\n', norm(res));
            
            % res1 = obj.A*obj.x_s + obj.E*obj.ss.I_ss + obj.B*obj.u_s + obj.E*obj.wBar;
            % fprintf('||steady-state residual||_2 = %.3e\n', norm(res1));
            % res2 = obj.ss.I_ss - obj.YBar* obj.D'*obj.x_s;
            % fprintf('||steady-state residual2||_2 = %.3e\n', norm(res2));

        end

        function [AdjMat, KMat, out] = design_MB_GSC(obj, linkFiltThresh, useData)
        % Design a dense K (N x 2N) s.t. A + B K is Hurwitz (continuous-time).
        % No sparsity constraints; we infer comm graph from K afterward.
            
            A = obj.A;  
            BBar = obj.BBar;
            YBar = obj.YBar;
            DBar = obj.DBar;
            D = obj.D;

            N = obj.N;  
            nX = 2*N;
        
            % Decision vars
            P = sdpvar(nX,nX,'symmetric');
            L = sdpvar(2*N, 2*N,'full');        % Y = K*P
            epsilon = sdpvar(1,1,'full');

            eps = 1e-3;
            cons = [P >= epsilon*eye(nX), epsilon >= eps];
            W = - (A*P + BBar*L)' - (A*P + BBar*L);
            cons = [cons, W >= epsilon*eye(size(W))]; % -2*opts.alpha*P
            % cons = [cons, W >= 2*eps*P]; % 2*opts.alpha*P
            cons = [cons, 1*(D'*L - YBar*D'*P)==zeros(size(D'*L))];
        
            costFun = epsilon + trace(P);

            ops = sdpsettings('solver','mosek','verbose',0);
            sol = optimize(cons, costFun, ops);
        
            out.problem = sol.problem;
            out.info    = yalmiperror(sol.problem);
            
            if sol.problem ~= 0
                warning('Model-based GSC design failed: %s', out.info);
                obj.K = [];
                AdjMat = []; 
                KMat = [];
                return;
            end

            disp('Model-based GSC design SUCCESS!');

            K = value(L) / value(P);           % K = Y * P^{-1}
            % obj.K_L = zeros(N,2*N);   % strictly local component off for this test
            for i = 1:1:N
                obj.DGs(i).K = [0, 0];
            end

            K_G = DBar'*K;             % store full K
            if useData
                mag = 2;
                K_G = K_G + K_G.*(mag*(2*rand(size(K_G))-1));
            end
            obj.K = K_G;

            % eig(A + BBar*K)
            % phyError = norm(value(D'*L - YBar*D'*P))

            [AdjMat, KMat] = obj.buildCommAdjFromK(linkFiltThresh);   % threshold for nonzero blocks
            
        end


                
        function [Adj, K] = buildCommAdjFromK(obj, thr)
        % Build an undirected adjacency by thresholding 1x2 block norms of K.
            N = obj.N; 
            Adj = zeros(N);
            K = obj.K;
            maxK = max(max(abs(K)));
            for i = 1:N
                for j = 1:N
                    K_ij = K(i, 2*j-1:2*j);
                    if norm(K_ij,2) >= thr*maxK
                        Adj(i,j) = 1;
                    else
                        K(i, 2*j-1:2*j) = [0, 0];
                    end
                end
            end

            obj.K = K;
            obj.commAdj = Adj;
        end



        % function out = designLocalK(obj, useData)
        % % Run local dissipativity on each DG;
        %     out = [];
        %     for i = 1:obj.N
        % 
        %         if ~useData
        %             oi = obj.DGs(i).designLocalXiDissipative();
        %         else
        %             oi = obj.DGs(i).designLocalXiDissipative_DataDriven();
        %         end
        % 
        %         if oi.problem~=0
        %             warning('Local design failed at DG %d: %s', i, oi.info);
        %         else
        %             out = [out; oi.K];
        %         end
        %     end
        % end
        

        function [AdjMat, KMat, out] = codesign_MB_DRC(obj, linkFiltThresh) 
            
            % Design local controllers
            for i = 1:obj.N
                oi = obj.DGs(i).designLocalXiDissipative();
                if oi.problem~=0
                    warning('Local design failed at DG %d: %s', i, oi.info);
                end
            end


            YBar = obj.YBar;
            % YBar = YBar - 0*diag(diag(YBar));
            BBar = obj.BBar;
            N = obj.N; 
            D = obj.D;
            DBar = obj.DBar;
        
            % Whether to use a soft or hard graph constraint
            isSoft = 1;
            normType = 1;
            maxCostVal = 0.001;
            % minCostVal = 0.001;
            
            % Set up the LMI problem
            I = eye(2*N);
            I_n = eye(2);
            O = zeros(2*N);

            % Variables
            K = sdpvar(2*N, 2*N,'full'); 
            P = sdpvar(N, N, 'diagonal');
            gammaSq = sdpvar(1, 1, 'full');
            epsilon = sdpvar(1, 1);    


            X_p_11 = [];
            X_11 = [];
            X_p_12 = [];
            X_12 = [];
            X_p_22 = [];

            for i = 1:1:N
                nu_i = obj.DGs(i).nu;
                rho_i = obj.DGs(i).rho;
            
                X_p_11 = blkdiag(X_p_11, -nu_i*P(i,i)*I_n);
                X_11 = blkdiag(X_11, -nu_i*I_n);
                X_p_12 = blkdiag(X_p_12, 0.5*P(i,i)*I_n);
                X_12 = blkdiag(X_12, (-1/(2*nu_i))*I_n);
                X_p_22 = blkdiag(X_p_22, -rho_i*P(i,i)*I_n);
            end
            X_p_21 = X_p_12';
            X_21 = X_12';
            
            costMat = [];
            for i = 1:1:N
                costMatRow = [];
                for j = 1:1:N
                    dist_ij = norm(obj.DGs(i).pos-obj.DGs(j).pos);
                    if dist_ij==0
                        dist_ij = 0.001;  %%%% Key point check
                    end
                    costMatRow = [costMatRow, dist_ij*ones(2)];
                end
                costMat = [costMat; costMatRow];
            end
            costMat = costMat;

            % Objective Function
            % costFun = 1*norm(Q.*costMatBlock,normType);
            KMat1 = DBar'*(K.*costMat)*D;
            KMat2 = DBar'*(K.*costMat)*DBar;
            KMat = KMat1 + KMat2;
            % costFun00 = sum(sum(KMat));
            costFun0 = norm(KMat,normType);

            % Minimum Budget Constraints
            con0 = [];
            % con0 = [con0, costFun00 >= minCostVal];
            con0 = [con0, costFun0 <= maxCostVal, gammaSq >= 0.001];
           
                        
            % Basic Constraints
            con1 = [P >= 0.001*eye(N)];

            %% Since: 
            % KBar = D*K + DBar*YBar*D';
            % L_uy = X_p_11*(BBar * KBar);
            %% We have with KHat = X_p_11*BBar*D*K = (struct of DK)
            L_uy = X_11*BBar*K;

            DMat = [X_p_11, O; O, I];
            MMat = [L_uy, X_p_11; I, O];
            ThetaMat = [- X_21*L_uy - L_uy'*X_12 - X_p_22, - X_p_21; - X_p_12, gammaSq*I];
            W = [DMat, MMat; MMat', ThetaMat];
            con2 = W >= epsilon*eye(size(W)); % The real one

            con3 = 1*(D'*K - P*YBar*D')==zeros(size(D'*K));
                        
            % Total Cost and Constraints
            if isSoft
                cons = [con0, con1, con2, con3]; % Without the hard graph constraint con7
                costFun = 1*costFun0 + 1*gammaSq + 1*trace(P) - 1e12*epsilon; % soft 
            else
                cons = [con0, con1, con2, con3]; % With the hard graph constraint con7
                costFun = 1*costFun0 + 1*gammaSq + 1*trace(P); % hard (same as soft)
            end
            
            % Solving
            ops = sdpsettings('solver','mosek','verbose',0);
            sol = optimize(cons, costFun, ops);
        
            out.problem = sol.problem;
            out.info    = yalmiperror(sol.problem);
        
            if sol.problem ~= 0
                warning('Model-based DRC co-design failed: %s', out.info);
                obj.K = [];
                AdjMat = []; 
                KMat = [];
                return;
            end
        
            disp('Model-based DRC co-design SUCCESS!');
            [k_crit, minor_crit, minors] = obj.criticalLeadingMinor(value(W))

            PVal        = value(P);
            KVal        = value(K);
            costFun0Val = value(costFun0);
            gammaSqVal  = value(gammaSq);
        
            fprintf('epsilon    = %.4e\n', value(epsilon));
            fprintf('||K||-weighted    = %.4e\n', costFun0Val);
            fprintf('gamma^2 (global)  = %.4e\n', gammaSqVal);
            fprintf('trace(P)          = %.4e\n', trace(PVal));

            
            % Consistency of DᵀK = P Ȳ Dᵀ
            phyErr = norm(value(D'*K - P*YBar*D')) 
            fprintf('||D''K - PȲD''|| = %.4e\n', phyErr);
            
            KBarVal = PVal \ DBar'*KVal
            obj.K = KBarVal;

            [AdjMat, KMat] = obj.buildCommAdjFromK(linkFiltThresh);   % threshold for nonzero blocks

        end

        
        function [AdjMat, KMat, out] = codesign_DD_DRC(obj, linkFiltThresh)
            
            % Local controller design
            for i = 1:1:obj.N
                oi = obj.DGs(i).designLocalXiDissipative_DataDriven();
                if oi.problem~=0
                    warning('Local design failed at DG %d: %s', i, oi.info);
                end
            end
            
            %--------------------------------------------------------------
            % Build global system matrices (same as model-based version)
            %--------------------------------------------------------------
            
            YBar = obj.YBar;
            % BBar = obj.BBar;
            N    = obj.N;
            D    = obj.D;
            DBar = obj.DBar;
        
            % Aggregate data-driven noise information
            QBar_w = obj.QBar_w;
            E = obj.E_perm;

            scaleQ = max(1, max(abs(QBar_w(:))))   % e.g. ~1e4 for your case
            QBar_w = QBar_w / scaleQ;


            %--------------------------------------------------------------
            % Graph / communication cost parameters
            %--------------------------------------------------------------
            isSoft     = 1;
            normType   = 1;
            maxCostVal = 0.001; %%%% check
        
            %--------------------------------------------------------------
            % Basic identities
            %--------------------------------------------------------------
            nN = 2*N;       % total state dim (each DG has 2 states)
            I_nN    = eye(nN);
            O_nN    = zeros(nN);
            I_n = eye(2);
        
            %--------------------------------------------------------------
            % Decision variables
            %--------------------------------------------------------------
            K       = sdpvar(nN, nN, 'full');    % pre-scaled global K
            P       = sdpvar(N, N, 'diagonal');          % p_i > 0 (scaling)
            gammaSq = sdpvar(1, 1, 'full');             % global Y-dissipativity gain
            epsilon = sdpvar(1, 1);                     % small slack
            lambda  = sdpvar(1, 1);                     % data-driven robust multiplier
        
            %--------------------------------------------------------------
            % Build X_p^{kl} and X^{kl} from local (nu_i, rho_i)
            %   X_i^{11} = -nu_i I_2,  X_i^{12} = 0.5 I_2,  X_i^{22} = -rho_i I_2
            %   X_p^{kl} = diag(p_i X_i^{kl})
            %--------------------------------------------------------------
            X_p_11 = [];
            X_11   = [];
            X_p_12 = [];
            X_12   = [];
            X_p_22 = [];
        
            for i = 1:N
                nu_i  = obj.DGs(i).nu;
                rho_i = obj.DGs(i).rho;
        
                X_p_11 = blkdiag(X_p_11, -nu_i*P(i,i)*I_n);
                X_11   = blkdiag(X_11,   -nu_i*I_n);
        
                X_p_12 = blkdiag(X_p_12,  0.5*P(i,i)*I_n);
                X_12   = blkdiag(X_12,  (-1/(2*nu_i))*I_n);
        
                X_p_22 = blkdiag(X_p_22, -rho_i*P(i,i)*I_n);
            end
            X_p_21 = X_p_12';
            X_21   = X_12';

            XBar_11 = eye(nN)/X_11;
            XBar_p_11 = XBar_11*X_p_11*XBar_11;     

            Y_11 = gammaSq*I_nN;
            Y_12 = O_nN;
            Y_21 = O_nN;
            Y_22 = -I_nN;

            %--------------------------------------------------------------
            % Communication cost matrix 
            %--------------------------------------------------------------
            costMat = [];
            for i = 1:N
                row = [];
                for j = 1:N
                    dist_ij = norm(obj.DGs(i).pos - obj.DGs(j).pos);
                    if dist_ij == 0
                        dist_ij = 0.001;   % avoid zero weight on self-links
                    end
                    row = [row, dist_ij*ones(2)];
                end
                costMat = [costMat; row];
            end
        
            % weighted controller matrix used in objective
            KMat1   = DBar'*(K.*costMat)*D;
            KMat2   = DBar'*(K.*costMat)*DBar;
            KMat    = KMat1 + KMat2;
            costFun0 = norm(KMat, normType);
        

            %--------------------------------------------------------------
            % Global LMI (data-driven robust version of Prop. 4)
            %--------------------------------------------------------------
            % Closed-loop term L_uy = 𝓧^{11} B K̄
            Q = [Y_11, O_nN; O_nN, -Y_22];

            S_11 = [X_p_11*XBar_11, O_nN, O_nN];
            S_12 = [-X_p_12+Y_12, O_nN, O_nN];
            S_21 = [O_nN, O_nN, O_nN];
            S_22 = [-Y_22, O_nN, O_nN];
            O_S = [O_nN, O_nN, O_nN];
            S = [S_11, S_12, O_S; 
                 S_21, S_22, O_S];%2nNx9nN
            
            R_11 = [XBar_p_11, O_nN, O_nN;
                       O_nN,   O_nN, O_nN;
                       O_nN,   O_nN, O_nN];
            R_12 = [O_nN,   O_nN,   O_nN;
                    O_nN,   O_nN,   O_nN;
                    K,   O_nN,   O_nN];
            R_21 = R_12';
            R_22 = [X_p_22, O_nN, O_nN;
                     O_nN,   O_nN, O_nN;
                     O_nN,   O_nN, O_nN];
            O_R = zeros(3*nN);

            R1 = [R_11, R_12, O_R;
                 R_21, -R_22, -R_21;
                 O_R,  -R_12, O_R];
            
            EQ_wMat = [E*QBar_w*E', O_R,    O_R;
                        O_R,        O_R,    O_R;
                        O_R,        O_R,    E*QBar_w*E'];

            % EQ_wMat = [E*QBar_w*E', O_R,    O_R;
            %             O_R,        O_R,    0.5*E*QBar_w*E';
            %             O_R,        0.5*E*QBar_w*E',    O_R];
            
            R = R1 - lambda*EQ_wMat;
            

            % Full LMI matrix
            W = [Q,  S;
                 S', R];
        

            %--------------------------------------------------------------
            % Constraints
            %--------------------------------------------------------------
            cons = [];
        
            % Communication budget
            cons = [cons, costFun0 <= maxCostVal, gammaSq >= 0.001]; %%%% check
        
            % Positivity / slacks 
            cons = [cons, P >= 0.001*eye(N), lambda >= 0];
        
            % Main Y-dissipativity LMI (data-driven)
            cons = [cons, W >= epsilon*eye(size(W))];
        
            % Structural coupling constraint: DᵀK = P Ȳ Dᵀ
            cons = [cons, 1*(D'*K - P*YBar*D')==zeros(size(D'*K))]; %%%% check
        
            %--------------------------------------------------------------
            % Objective: sparse K, small gamma, small λ and P
            %--------------------------------------------------------------
            if isSoft
                costFun = 1*costFun0 + 1*gammaSq + 1*trace(P) - 1e12*epsilon; %%%% check
            else
                % cons = [cons, ...]
                % costFun = 1*costFun0 + 1*gammaSq + 1e3*lambda + 1*trace(P);
            end
        
            %--------------------------------------------------------------
            % Solve
            %--------------------------------------------------------------
            ops = sdpsettings('solver','mosek','verbose',0);
            sol = optimize(cons, costFun, ops);
        
            out.problem = sol.problem;
            out.info    = yalmiperror(sol.problem);
        
            if sol.problem ~= 0
                warning('Data-driven DRC co-design failed: %s', out.info);
                obj.K = [];
                AdjMat = []; 
                KMat = [];
                return;
            end
        
            disp('Data-driven DRC co-design SUCCESS!');
        
            [k_crit, minor_crit, minors] = obj.criticalLeadingMinor(value(W))

            %--------------------------------------------------------------
            % Extract and store results
            %--------------------------------------------------------------
            PVal        = value(P);
            KVal        = value(K);
            costFun0Val = value(costFun0);
            gammaSqVal  = value(gammaSq);
            lambdaVal   = value(lambda);
        
                    
            fprintf('epsilon    = %.4e\n', value(epsilon));
            fprintf('||K||-weighted    = %.4e\n', costFun0Val);
            fprintf('gamma^2 (global)  = %.4e\n', gammaSqVal);
            fprintf('trace(P)          = %.4e\n', trace(PVal));
            fprintf('lambda (robust)   = %.4e\n', lambdaVal);

            % Consistency of DᵀK = P Ȳ Dᵀ
            phyErr = norm(value(D'*K - P*YBar*D'));
            fprintf('||D''K - PȲD''|| = %.4e\n', phyErr);
        
            % Recover K̄ as in model-based code
            KBarVal   = PVal \ (DBar'*KVal)
            obj.K     = KBarVal;

            [AdjMat, KMat] = obj.buildCommAdjFromK(linkFiltThresh);   % threshold for nonzero blocks
        
            
        end


        function [k_crit, minor_crit, minors] = criticalLeadingMinor(obj, M)
            
            M;
            n = size(M,1);
            minors = zeros(n,1);
        
            for k = 1:n
                Mk = M(1:k, 1:k);
                kVal = k
                detVal = det(Mk)
                minors(k) = detVal;
            end
        
            [minor_crit, k_crit] = min(minors);

        end


    end
end

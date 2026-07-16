classdef DCMicrogrid < handle
    % DCMicrogrid Network assembly, simulation, and controller co-design.
    %
    % Runtime role: key system component used by main.mlx. It assembles DG
    % and line models, computes operating points, runs staged simulations,
    % records data matrices, implements model/data-driven global designs,
    % and draws the physical and communication topologies.
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

                % REVIEW PROPOSAL SYS-01 (ADD): validate endpoints and line
                % conductance before indexing the network matrices.
                % assert(i>=1 && i<=N && j>=1 && j<=N && i~=j, ...
                %     'DCMicrogrid:InvalidLineEndpoints','Invalid line endpoints.');
                % assert(isfinite(g) && g>0, 'DCMicrogrid:InvalidConductance', ...
                %     'Every physical line must have positive finite conductance.');

                Y(i,j) = g; 
                Y(j,i) = g;
                
                YBar(i,j) = -g; 
                YBar(j,i) = -g;

                % REVIEW PROPOSAL SYS-02 (REPLACE): if parallel physical
                % lines are allowed, replace the four assignments above by
                % accumulation so off-diagonal and diagonal terms agree.
                % Y(i,j) = Y(i,j)+g;       Y(j,i) = Y(j,i)+g;
                % YBar(i,j) = YBar(i,j)-g; YBar(j,i) = YBar(j,i)-g;

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

            % REVIEW PROPOSAL SIM-01 (ADD): uncomment this block to make the
            % event schedule correct for a tspan that does not start at zero.
            % duration = t_f-t_0;
            % t_1=t_0+0.05*duration; t_2=t_0+0.20*duration;
            % t_3=t_0+0.40*duration; t_4=t_0+0.60*duration;
            % t_5=t_0+0.65*duration; t_6=t_0+0.80*duration;
            
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
            mag = 0.1;
            for i = 1:obj.N
                if rand(1)<0.5
                    X(end, 2*i-1:2*i) = X(end, 2*i-1:2*i) + ...
                        X(end, 2*i-1:2*i).*(mag*(2*[rand(1,1), 1] - 2));

                    % REVIEW PROPOSAL SIM-02 (REPLACE): the active expression
                    % only decreases voltage and applies exactly zero current
                    % disturbance. Replace it with a symmetric two-state draw.
                    % X(end,2*i-1:2*i) = X(end,2*i-1:2*i).*( ...
                    %     1 + mag*(2*rand(1,2)-1));
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
            mag = 500;
            for i = 1:obj.N
                tIndices = (obj.DGs(i).noise.t > t_3 & obj.DGs(i).noise.t < (0.98*t_3 + 0.02*t_4));
                % REVIEW PROPOSAL SIM-03 (ADD): uncomment to make the
                % heightened disturbance span all of Stage 4 instead of 2%.
                % tIndices = obj.DGs(i).noise.t > t_3 & obj.DGs(i).noise.t <= t_4;
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
            mag = 0.4;
            for i = 1:obj.N
                obj.DGs(i).A    = obj.DGs(i).A + obj.DGs(i).A.*(mag*(2*rand(size(obj.DGs(i).A))-1));
                obj.DGs(i).Ibar = obj.DGs(i).Ibar + obj.DGs(i).Ibar.*(mag*(2*rand(size(obj.DGs(i).Ibar))-1));

                % REVIEW PROPOSAL SIM-04 (ADD): uncomment this structured
                % perturbation block to overwrite the arbitrary A perturbation
                % using physically meaningful component changes instead.
                % obj.DGs(i).L  = obj.DGs(i).L *(1+mag*(2*rand()-1));
                % obj.DGs(i).C  = obj.DGs(i).C *(1+mag*(2*rand()-1));
                % obj.DGs(i).Rf = obj.DGs(i).Rf*(1+mag*(2*rand()-1));
                % obj.DGs(i).RL = obj.DGs(i).RL*(1+mag*(2*rand()-1));
                % obj.DGs(i).updateModel();
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
                [AdjMat, KMat, out] = obj.codesign_DD_DRC(linkFiltThresh); %%%% Check
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
            mag = 0.1;
            for i = 1:obj.N
                if rand(1)<0.5
                    X(end, 2*i-1:2*i) = X(end, 2*i-1:2*i) + ...
                        X(end, 2*i-1:2*i).*(mag*(2*[rand(1,1), 1] - 2));

                    % REVIEW PROPOSAL SIM-05 (REPLACE): same correction as
                    % SIM-02 for the final-stage state disturbance.
                    % X(end,2*i-1:2*i) = X(end,2*i-1:2*i).*( ...
                    %     1 + mag*(2*rand(1,2)-1));
                end
            end
            % Disturbance amplification
            mag = 250;
            for i = 1:obj.N
                tIndices = (obj.DGs(i).noise.t > t_6 & obj.DGs(i).noise.t < (0.98*t_6 + 0.02*t_f));
                % REVIEW PROPOSAL SIM-06 (ADD): uncomment for a full Stage-7
                % disturbance interval rather than the first 2% only.
                % tIndices = obj.DGs(i).noise.t > t_6 & obj.DGs(i).noise.t <= t_f;
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

            % meanAbsDistVals = [];
            % meanAbsDiscEVals = [];
            QBar_w = [];
            E1 = []; E2 = []; E3 = [];
            N = obj.N;

            % Data sampling
            Ts = 1e-5;
            obj.Ts = Ts;
            tVals = t(1):Ts:t(end);

            for k = 1:N
            
                obj.DGs(k).Ts = Ts;

                % Loading xTIlde and uTilde Data
                x_k = X(:, 2*k-1:2*k);
                x_sk = Xs(:, 2*k-1:2*k);
                xTilde = x_k - x_sk; %Xs(2*k-1:2*k,1)'
                uTilde = Utilde(:, 2*k-1:2*k);
                % wTilde = Wtilde(:, 2*k-1:2*k);
            
                xTildeSampled = interp1(t, xTilde, tVals, 'linear')';
                uTildeSampled = interp1(t, uTilde, tVals, 'previous', 'extrap')';
                % wTildeSampled = Ts*interp1(t, wTilde, tVals, 'previous', 'extrap')';
            
                % Data loading
                obj.DGs(k).xTilde = xTildeSampled(:,1:end-1);
                obj.DGs(k).yTilde = obj.DGs(k).xTilde;
                obj.DGs(k).xTildeBar = xTildeSampled(:,2:end);
                obj.DGs(k).uTilde = uTildeSampled(:,1:end-1);
                % obj.DGs(k).wTilde = wTildeSampled(:,1:end-1);

                % REVIEW PROPOSAL DATA-01 (ADD): verify Assumption 5 before
                % constructing either Proposition 7 or Proposition 8.
                % dataStack = [obj.DGs(k).uTilde;obj.DGs(k).xTilde];
                % sData = svd(dataStack);
                % rankTol = max(size(dataStack))*eps(max(sData));
                % assert(sum(sData>rankTol)==size(dataStack,1), ...
                %     'DCMicrogrid:DataRank','DG %d violates Assumption 5.',k);
                
                % Disturbance impact
                % meanAbsDistVal = mean(abs(obj.DGs(k).wTilde')); 
                % meanAbsDistVals = [meanAbsDistVals; meanAbsDistVal];
                
                % discretization Error Impact
                discError = obj.DGs(k).xTildeBar - ( (eye(2) + Ts*obj.DGs(k).A)*obj.DGs(k).xTilde ...
                                        + Ts*obj.DGs(k).BBar*obj.DGs(k).uTilde);
                % meanAbsDiscEVal = mean(abs(discError')); 
                % meanAbsDiscEVals = [meanAbsDiscEVals; meanAbsDiscEVal];
                
                % Loading the discretization error as the disturbance 
                obj.DGs(k).wTilde = discError; %%%% check

                % REVIEW PROPOSAL DATA-02 (REPLACE): the active residual uses
                % the supposedly unknown A and BBar. For a simulation-only
                % check, use the logged disturbance impact instead. The Ts
                % scaling below is forward-Euler only; exact integration is
                % preferable when the disturbance is held over the interval.
                % wImpact = Wtilde(:,2*k-1:2*k);
                % wImpactSampled = Ts*interp1(t,wImpact,tVals,'previous','extrap')';
                % obj.DGs(k).wTilde = wImpactSampled(:,1:end-1);

                % % Approach 1: Finding an upper bound for disturbance: Q_w 
                % wwT = obj.DGs(k).wTilde*obj.DGs(k).wTilde';
                % lambda = max(eig(wwT))*1.001; %%%% chech
                % % eig(lambda*eye(2)-wwT);
                % Q_ww = -eye(length(tVals)-1);
                % Q_I = lambda*eye(2);
                % Q_z = zeros(2,length(tVals)-1);
                % Q_w = [Q_I, Q_z; Q_z', Q_ww];
                % obj.DGs(k).Q_w = Q_w;

                % % Approach 2: Compute Q_w
                Q_w = obj.compute_Qw_from_wtilde(discError);
                % REVIEW PROPOSAL DATA-03 (ADD): when DATA-02 is enabled,
                % uncomment this override. For the paper claim, replace both
                % fitted bounds by an independently specified a-priori Q_w.
                % Q_w = obj.compute_Qw_from_wtilde(obj.DGs(k).wTilde);
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
            
            % meanAbsDist = mean(meanAbsDistVals);
            % meanAbsDiscErr = mean(meanAbsDiscEVals);

        end

        function Qw_val = compute_Qw_from_wtilde(obj, W)

            % Given: W (nW x T)
            [nW, T] = size(W);
            
            % Identity over time
            I_n = eye(nW);
            
            % Decision variable Q_w (fully general symmetric for now)
            
            Q_I = sdpvar(nW, nW, 'symmetric');
            Q_z = zeros(nW,T); %sdpvar(nW, 1, 'full')*ones(1, T);
            Q_ww = -eye(T);
            Q_w = [Q_I, Q_z; Q_z', Q_ww];
            epsilon = sdpvar(1, 1, 'full');

            % Construct the stacked data matrix
            Phi = [I_n; W'];    % size: (T + nW*T) x T
            
            % QMI constraint
            M = Phi' * Q_w * Phi;    % size: T x T
            
            Constraints = [M >= epsilon*eye(size(M)), epsilon >= 1e-3, Q_I >= epsilon*eye(size(Q_I))];
            
            % To avoid the trivial Q_w = 0, add some objective.
            % Example: minimize Frobenius norm of Q_w, or add structure.
            Objective = trace(Q_I) + epsilon;   % or something else
            
            opts = sdpsettings('solver','mosek','verbose',0);
            sol  = optimize(Constraints, Objective, opts);
            
            if sol.problem ~= 0
                error('SDP infeasible or solver failed: %s', sol.info);
            end
            
            Qw_val = value(Q_w);

            % REVIEW PROPOSAL DATA-04 (ADD): a bound fitted to the exact
            % observed realization has no reserve. Inflate its state block or,
            % preferably, supply a physical bound selected before collecting data.
            % disturbanceMargin = 1.20;
            % Qw_val(1:nW,1:nW) = disturbanceMargin^2*Qw_val(1:nW,1:nW);

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
            Utilde = zeros(nt, 2*N);
            Wtilde = zeros(nt, 2*N);
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

                    DG = obj.DGs(i);
        
                    % DG.x was set by setStateVector
                    x_i  = DG.x;                       % local state of DG i

                    Iline_i = Ik(i);
                    % Saturation
                    Imax = 2*DG.Irated; Imin = -2*DG.Irated;
                    Iline_i = min(Imax, max(Imin, Iline_i));
                    % REVIEW PROPOSAL DATA-05 (ADD): keep recorded input data
                    % consistent with the network current used by the theory.
                    % Iline_i = Ik(i);
            

                    u_Si = DG.u_s;
                    if useData
                        u_Li = DG.u_L_Hold;
                    else
                        u_Li = DG.K * (x_i - DG.x_s);      % local input
                    end
                    u_Gi = u_Gk(i);
                    u_i  = u_Si + u_Li + u_Gi;     % total input
                    % Saturation
                    Vmax = 2*DG.Vrated; Vmin = -2*DG.Vrated;
                    u_i = min(Vmax, max(Vmin, u_i));
                    % REVIEW PROPOSAL DATA-06 (ADD): disable command clipping
                    % for the certified linear run, or model/log saturation.
                    % u_i = u_Si + u_Li + u_Gi;
        
                    ITilde_ki = Iline_i - I_S(i);
                    uTilde_ki = u_i - u_Si;
                    
                    tVal = t(k);
                    w_i = interp1(DG.noise.t, DG.noise.w, tVal, 'previous', 'extrap')';
                    Wtilde_ki = (DG.BBar*w_i);

                    U(k,2*i-1:2*i) = [Iline_i, u_i];
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
            % sigma is currently used as a standard-deviation/factor matrix,
            % not as a covariance matrix or a scalar disturbance variance.

            % REVIEW PROPOSAL DATA-07 (ADD): validate the actual two-channel
            % disturbance factor expected by randn(...,2)*sigma.
            % validateattributes(dt_noise,{'numeric'},{'scalar','positive','finite'});
            % assert(isequal(size(sigma),[2 2]) && all(isfinite(sigma(:))), ...
            %     'DCMicrogrid:NoiseScale','sigma must be a finite 2-by-2 factor.');

            for i = 1:obj.N
                w = randn(length(t_noise), 2)*sigma;
                % REVIEW PROPOSAL DATA-08 (REPLACE): if the caller supplies a
                % covariance matrix instead, use its Cholesky factor.
                % noiseFactor = chol(sigma,'lower');
                % w = randn(length(t_noise),2)*noiseFactor';
                obj.DGs(i).noise.t = t_noise;
                obj.DGs(i).noise.w = w;
            end

        end


        


        function draw(obj, ax, varargin)

            if nargin < 2 || isempty(ax), ax = gca; end

            parser = inputParser;
            parser.FunctionName = 'DCMicrogrid.draw';
            addParameter(parser,'Title','DC Microgrid Topology', ...
                @(value)ischar(value) || (isstring(value) && isscalar(value)));
            addParameter(parser,'ShowCommunication',true, ...
                @(value)islogical(value) && isscalar(value));
            addParameter(parser,'ShowLineCurrents',false, ...
                @(value)islogical(value) && isscalar(value));
            addParameter(parser,'ShowLineLabels',true, ...
                @(value)islogical(value) && isscalar(value));
            addParameter(parser,'ShowNodeLabels',true, ...
                @(value)islogical(value) && isscalar(value));
            addParameter(parser,'FontSize',8, ...
                @(value)isnumeric(value) && isscalar(value) && value>0);
            parse(parser,varargin{:});
            options = parser.Results;

            cla(ax); hold(ax,'on'); axis(ax,'equal');
            set(ax,'Color','w','XColor',[0.25 0.25 0.25], ...
                'YColor',[0.25 0.25 0.25]);

            positions = vertcat(obj.DGs.pos);
            networkCenter = mean(positions,1);
            xSpan = max(max(positions(:,1))-min(positions(:,1)),1);
            ySpan = max(max(positions(:,2))-min(positions(:,2)),1);

            voltage = zeros(obj.N,1);
            for k = 1:obj.N
                voltage(k) = obj.DGs(k).x(1);
            end
            networkCurrent = obj.YBar*voltage;

            % Draw physical lines first so communication links remain visible.
            for e = 1:obj.M
                i = obj.Lines(e).i; j = obj.Lines(e).j;
                lineCurrent = obj.Lines(e).current(voltage(i),voltage(j));
                obj.Lines(e).draw(ax,obj.DGs(i).pos,obj.DGs(j).pos, ...
                    'ShowLabel',options.ShowLineLabels, ...
                    'ShowCurrent',options.ShowLineCurrents, ...
                    'LineCurrent',lineCurrent,'FontSize',options.FontSize);
            end

            if options.ShowCommunication
                obj.drawComm(ax);
            end

            for k = 1:obj.N
                if positions(k,2) >= networkCenter(2)+0.25*ySpan
                    labelPosition = positions(k,:)+[0,0.2];
                    horizontalAlignment = 'center';
                    verticalAlignment = 'bottom';
                elseif positions(k,1) <= networkCenter(1)
                    labelPosition = positions(k,:)+[-0.16,-0.14];
                    horizontalAlignment = 'right';
                    verticalAlignment = 'top';
                else
                    labelPosition = positions(k,:)+[0.16,-0.14];
                    horizontalAlignment = 'left';
                    verticalAlignment = 'top';
                end
                obj.DGs(k).draw(ax,'NetworkCurrent',networkCurrent(k), ...
                    'ShowLabels',options.ShowNodeLabels, ...
                    'FontSize',options.FontSize, ...
                    'LabelPosition',labelPosition, ...
                    'HorizontalAlignment',horizontalAlignment, ...
                    'VerticalAlignment',verticalAlignment);
            end

            xlim(ax,[min(positions(:,1))-max(1.6,0.45*xSpan), ...
                max(positions(:,1))+max(1.6,0.45*xSpan)]);
            ylim(ax,[min(positions(:,2))-max(0.9,0.3*ySpan), ...
                max(positions(:,2))+max(0.9,0.3*ySpan)]);

            title(ax,char(options.Title),'Color',[0.08 0.18 0.28], ...
                'FontWeight','bold');
            set(ax,'Visible','off');
            ax.Title.Visible = 'on';

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

                        % REVIEW PROPOSAL DRAW-01 (ADD): protect the
                        % normalization below if two DG positions coincide.
                        % if norm(pj-pi) <= eps, continue; end
                        
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

            % REVIEW PROPOSAL SS-01 (ADD): Eq. (45) bounds the equilibrium
            % exported current and VSC voltage command as well as the state.
            % deltaI = 1.0; deltaVt = 0.20;
            % constr = [constr, -deltaI.*Ir <= I, I <= deltaI.*Ir];
            % constr = [constr, (1-deltaVt).*Vr <= u, u <= (1+deltaVt).*Vr];
        
            % Mild input regularization
            % us_prev = mats.u_s; if isempty(us_prev), us_prev = zeros(N,1); end
        
            objFun = norm(epsV - 1, 2)^2 + epsI;
            % REVIEW PROPOSAL SS-02 (ADD): uncomment to override the active
            % objective, which currently drives epsI toward zero.
            % objFun = norm(epsV-1,2)^2 + (epsI-1)^2 + ...
            %          1e-4*norm((u-Vr)./Vr,2)^2;
        
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

            % REVIEW PROPOSAL SS-03 (ADD): verify Eqs. (46)-(47) numerically.
            % ssResidual = Am*ss.x_ss + [Em Bm]*[ss.I_ss;ss.u_s] + Em*wBar;
            % lineResidual = ss.I_ss-YBar*ss.V_ss;
            % assert(norm(ssResidual,inf)<=1e-7 && norm(lineResidual,inf)<=1e-7, ...
            %     'DCMicrogrid:SteadyStateResidual','Steady-state equations failed.');

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
            % REVIEW PROPOSAL GSC-01 (ADD): uncomment to disable gain
            % corruption. Collect persistently exciting additive input data
            % instead of perturbing a stabilizing gain by as much as 200%.
            % K_G = DBar'*K;
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
            % REVIEW PROPOSAL GSC-02 (ADD): avoid classifying every zero block
            % as an active edge when the entire recovered gain is zero.
            % if maxK==0, Adj=zeros(N); obj.commAdj=Adj; return; end
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

            % REVIEW PROPOSAL GSC-03 (DESIGN): thresholding changes the gain
            % after optimization, so the LMI certificate no longer applies.
            % Re-evaluate stability/dissipativity with this pruned K, or impose
            % the desired block sparsity inside the optimization and re-solve.

            obj.K = K;
            obj.commAdj = Adj;
        end



       
        

        function [AdjMat, KMat, out] = codesign_MB_DRC(obj, linkFiltThresh) 
            
            % Design local controllers
            for i = 1:obj.N
                oi = obj.DGs(i).designLocalXiDissipative();
                if oi.problem~=0
                    warning('Local design failed at DG %d: %s', i, oi.info);
                end
                % REVIEW PROPOSAL MB-G01 (ADD): do not continue to divisions
                % by nu_i or to a global certificate after a local failure.
                % assert(oi.problem==0, 'DCMicrogrid:LocalMBFailure', ...
                %     'Model-based local design failed at DG %d.',i);
            end

            YBar = obj.YBar
            BBar = obj.BBar;
            N = obj.N; 
            D = obj.D;
            DBar = obj.DBar;
        
            % Whether to use a soft or hard graph constraint
            isSoft = 1;
            normType = 1;
            maxCostVal = 1;
            
            % Set up the LMI problem
            I = eye(2*N);
            I_n = eye(2);
            O = zeros(2*N);

            % Variables
            P = sdpvar(N, N, 'diagonal');
            gammaSq = sdpvar(1, 1, 'full');
            epsilon = sdpvar(1, 1); 

            KHat = sdpvar(N, 2*N,'full');
            K = [];
            for i = 1:1:N
                K_i = [];
                for j = 1:1:N
                    K_ij = [P(i,i)*YBar(i,j), 0; 
                            KHat(i,2*j-1:2*j)];
                    K_i = [K_i, K_ij];
                end
                K = [K; K_i];
            end

            % P\D'*K*D = YBar iff YHat = P*YBar
            % D'*K = P*YBar*D'
            % KBar = P\DBar'*KVal iff KBar = P\KHat;
            % K = sdpvar(2*N, 2*N,'full');

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
                    dist_ij = 1*norm(obj.DGs(i).pos-obj.DGs(j).pos);
                    if dist_ij==0
                        dist_ij = 0; % 1e9; %Key point check
                    end
                    costMatRow = [costMatRow, dist_ij*ones(1,2)];
                end
                costMat = [costMat; costMatRow];
            end
            costMat = costMat

            % Objective Function
            KMat = KHat.*costMat;
            costFun0 = norm(KMat(:),normType);

            % Minimum Budget Constraints
            con0 = [];
            con0 = [con0, costFun0 <= maxCostVal];
            con0 = [con0, epsilon == 1e-6, gammaSq <= 1e6, gammaSq >= 0];
                        
            % Basic Constraints
            con1 = [P >= epsilon*eye(N)];
            limitVal = 10;
            con1 = [con1, -P*limitVal*ones(size(KHat)) <= KHat, KHat <= P*limitVal*ones(size(KHat))];
            % REVIEW PROPOSAL MB-G02 (ADD): if using -alpha*trace(P) as in
            % Remark 6, add an upper bound to prevent unbounded scaling.
            % pMax = 1e3; con1 = [con1, P <= pMax*eye(N)];

            % Main LMI
            L_uy = X_11*BBar*K;
            DMat = [X_p_11, O; O, I];
            MMat = [L_uy, X_p_11; I, O];
            ThetaMat = [- X_21*L_uy - L_uy'*X_12 - X_p_22, - X_p_21; - X_p_12, gammaSq*I];
            W = [DMat, MMat; MMat', ThetaMat];
            
            W_dim = size(W,1);
            epsilonSlack = sdpvar(W_dim, W_dim, 'diagonal'); 

            con2 = (W - epsilonSlack) >= epsilon*eye(W_dim); % The real one
            con2 = [con2, epsilonSlack >= -1*eye(W_dim), epsilonSlack <= 1*eye(W_dim)];
            % REVIEW PROPOSAL MB-G03 (ADD, recommended): this hard constraint
            % restores the Proposition-4 certificate. Negative epsilonSlack
            % currently permits an infeasible W to be reported as success.
            % con2 = [con2, W >= epsilon*eye(W_dim)];
                      
            % Total Cost and Constraints
            if isSoft
                cons = [con0, con1, con2]; % Without the hard graph constraint con7
                costFun = 1*costFun0 + 0*gammaSq + 1*trace(P) + 0*epsilon + 1*norm(epsilonSlack,1); 
                % REVIEW PROPOSAL MB-G04 (REPLACE): align the objective with
                % Prop. 4/Remark 6 and optimize the performance certificate.
                % alphaP=1e-3; betaGamma=1;
                % costFun = costFun0-alphaP*trace(P)+betaGamma*gammaSq + ...
                %           norm(diag(epsilonSlack),1);
                % epsilonslackcoef is a Key point check
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
            if abs(minor_crit)>1e-6
                disp('Error in Minors!')
                k_crit
                minor_crit 
                minors
            end

            PVal        = value(P)
            KVal        = value(K);
            KHatVal = value(KHat)
            costFun0Val = value(costFun0);
            gammaSqVal  = value(gammaSq);
            epsilonSlackVal = diag(value(epsilonSlack))';
            [minEps,minEpsIdx] = min(epsilonSlackVal)
            [maxEps,maxEpsIdx] = max(epsilonSlackVal)
        
            fprintf('epsilon    = %.4e\n', value(epsilon));
            fprintf('epsilonSlackmin    = %.4e\n', minEps);
            fprintf('epsilonSlackmax    = %.4e\n', maxEps);
            fprintf('||K||-weighted    = %.4e\n', costFun0Val);
            fprintf('gamma^2 (global)  = %.4e\n', gammaSqVal);
            fprintf('trace(P)          = %.4e\n', trace(PVal));

            
            % Consistency of DᵀK = P Ȳ Dᵀ
            % phyErr = norm(value(D'*K - P*YBar*D')); 
            % fprintf('||D''K - PȲD''|| = %.4e\n', phyErr);
            
            KBarVal = PVal \ KHatVal;
            obj.K = KBarVal;

            [AdjMat, KMat] = obj.buildCommAdjFromK(linkFiltThresh)   % threshold for nonzero blocks

        end

        
        function [AdjMat, KMat, out] = codesign_DD_DRC(obj, linkFiltThresh)
            
            % Local controller design
            for i = 1:1:obj.N
                oi = obj.DGs(i).designLocalXiDissipative_DataDriven();
                if oi.problem~=0
                    warning('Local design failed at DG %d: %s', i, oi.info);
                end
                % REVIEW PROPOSAL DD-G01 (ADD): Proposition 8 requires every
                % local Proposition-7 certificate to succeed.
                % assert(oi.problem==0, 'DCMicrogrid:LocalDDFailure', ...
                %     'Data-driven local design failed at DG %d.',i);
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

            scaleQ = max(1, max(abs(QBar_w(:))));   % e.g. ~1e4 for your case
            QBar_w = QBar_w / scaleQ;
            % REVIEW PROPOSAL DD-G02 (DESIGN): assemble exactly the normalized
            % per-DG QBar_w blocks used by the local designs. A new aggregate
            % scaling is not equivalent when one global lambda multiplies all
            % subsystem blocks with different original magnitudes.


            %--------------------------------------------------------------
            % Graph / communication cost parameters
            %--------------------------------------------------------------
            isSoft     = 1;
            normType   = 1;
            maxCostVal = 1e-3; %%%% check
            % REVIEW PROPOSAL DD-G03 (DESIGN): the model-based budget is 1,
            % while the data-driven budget is 1e-3. Use the same dimensionless
            % communication budget before attributing a difference to data.
        
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
                costMatRow = [];
                for j = 1:N
                    dist_ij = norm(obj.DGs(i).pos - obj.DGs(j).pos);
                    if dist_ij == 0
                        dist_ij = 0;   % avoid zero weight on self-links
                    end
                    costMatRow = [costMatRow, dist_ij*ones(1,2)];
                end
                costMat = [costMat; costMatRow];
            end
            costMat = costMat;

            % weighted controller matrix used in objective
            % KMat1   = DBar'*(K.*costMat)*D;
            % KMat2   = DBar'*(K.*costMat)*DBar;
            % KMat    = KMat1 + KMat2;
            KMat = (DBar'*K).*costMat;
            costFun0 = norm(KMat(:),normType);

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
            
            XBar_11_disp = XBar_11;

            EQ_wMat = [E*QBar_w*E', O_R,    O_R;
                        O_R,        O_R,    O_R;
                        O_R,        O_R,    E*QBar_w*E'];


            QBar_w_disp = QBar_w;
            EQBarET_neg = -E*QBar_w*E';
            

            [k_crit, minor_crit, minors] = obj.criticalLeadingMinor(EQBarET_neg);
            if abs(minor_crit)>1e-6
                disp('Error in Minors 0!')
                if k_crit(1)>0
                    col_crit = EQBarET_neg(1:k_crit(1),k_crit(1))
                end
            end
            

            % EQ_wMat = [E*QBar_w*E', O_R,    O_R;
            %             O_R,        O_R,    0.5*E*QBar_w*E';
            %             O_R,        0.5*E*QBar_w*E',    O_R];
            
            R = R1 - lambda*EQ_wMat;
            

            % Full LMI matrix
            W = [Q,  S;
                 S', R];

            W_dim = size(W,1);
            epsilonSlack = sdpvar(W_dim, W_dim, 'diagonal'); 

            %--------------------------------------------------------------
            % Constraints
            %--------------------------------------------------------------
            cons = [];
        
            % Communication budget
            cons = [cons, costFun0 <= maxCostVal, gammaSq >= epsilon, epsilon >= 0.001, gammaSq<=1e-3]; %%%% check
            % REVIEW PROPOSAL DD-G04 (REPLACE): the active bounds force both
            % gammaSq and epsilon to exactly 1e-3. Decouple feasibility margin
            % from performance and normalize before selecting tolerances.
            % cons = [costFun0<=maxCostVal, gammaSq>=1e-8, ...
            %         epsilon>=1e-7, epsilon<=1e-3];
        
            % Positivity / slacks 
            cons = [cons, lambda >= 0, epsilonSlack >= -1*eye(W_dim), epsilonSlack <= 1*eye(W_dim)];
        
            % Main Y-dissipativity LMI (data-driven)
            cons = [cons, P >= epsilon*eye(N), W - epsilonSlack >= epsilon*eye(W_dim)];
            % REVIEW PROPOSAL DD-G05 (ADD, recommended): require the actual
            % Proposition-8 LMI, not only its softened version.
            % cons = [cons, W >= epsilon*eye(W_dim)];
        
            % Structural coupling constraint: DᵀK = P Ȳ Dᵀ
            cons = [cons, 1e6*(D'*K - P*YBar*D') == zeros(size(D'*K))]; %%%% check
            % REVIEW PROPOSAL DD-G06 (REPLACE): remove the 1e6 equality
            % multiplier. It is algebraically redundant and harms scaling.
            % cons = [cons, D'*K == P*YBar*D'];
        
            %--------------------------------------------------------------
            % Objective: sparse K, small gamma, small λ and P
            %--------------------------------------------------------------
            if isSoft
                costFun = 1*costFun0 + 1*gammaSq + 1*trace(P) + 1*epsilon + 1*norm(epsilonSlack,1); %%%% check
                % REVIEW PROPOSAL DD-G07 (REPLACE): after adding a P upper
                % bound, use the Prop. 8 objective direction and maximize a
                % bounded margin rather than minimizing epsilon.
                % alphaP=1e-3; betaGamma=1; betaMargin=1e-3;
                % costFun = costFun0-alphaP*trace(P)+betaGamma*gammaSq ...
                %     -betaMargin*epsilon+norm(diag(epsilonSlack),1);
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
        
            [k_crit, minor_crit, minors] = obj.criticalLeadingMinor(value(W));
            if abs(minor_crit)>1e-6
                disp('Error in Minors!')
                k_crit
                minor_crit 
                minors
            end

            %--------------------------------------------------------------
            % Extract and store results
            %--------------------------------------------------------------
            PVal        = value(P);
            KVal        = value(K);
            costFun0Val = value(costFun0);
            gammaSqVal  = value(gammaSq);
            lambdaVal   = value(lambda);
            epsilonSlackVal = diag(value(epsilonSlack))';
            [minEps,minEpsIdx] = min(epsilonSlackVal)
            [maxEps,maxEpsIdx] = max(epsilonSlackVal)

                    
            fprintf('epsilon    = %.4e\n', value(epsilon));
            fprintf('epsilonSlackmin    = %.4e\n', minEps);
            fprintf('epsilonSlackmax    = %.4e\n', maxEps);
            fprintf('||K||-weighted    = %.4e\n', costFun0Val);
            fprintf('gamma^2 (global)  = %.4e\n', gammaSqVal);
            fprintf('trace(P)          = %.4e\n', trace(PVal));
            

            % Consistency of DᵀK = P Ȳ Dᵀ
            phyErr = norm(value(D'*K - P*YBar*D'));
            fprintf('||D''K - PȲD''|| = %.4e\n', phyErr);
            fprintf('lambda (robust)   = %.4e\n', lambdaVal);
        
            % Recover K̄ as in model-based code
            KBarVal   = PVal \ (DBar'*KVal);
            obj.K     = KBarVal;

            [AdjMat, KMat] = obj.buildCommAdjFromK(linkFiltThresh)   % threshold for nonzero blocks
        
            
        end


        function [k_crit, minor_crit, minors] = criticalLeadingMinor(obj, M)

            % REVIEW PROPOSAL NUM-01 (REPLACE): determinants of leading
            % minors are poorly scaled and the current callers incorrectly
            % treat a positive minor far from zero as an error. Prefer:
            % minEigenvalue = min(eig((value(M)+value(M)')/2));
            % assert(minEigenvalue >= -1e-7, ...
            %     'DCMicrogrid:LMIViolation','LMI minimum eigenvalue is negative.');
            
            M;
            n = size(M,1);
            minors = zeros(1,n);
            firstFound = 0;
            k_crit1 = 0;

            for k = 1:n
                Mk = M(1:k, 1:k);
                kVal = k;
                detVal = det(Mk);
                minors(k) = detVal;
                if detVal < 0 && ~firstFound
                    k_crit1 = k;
                    firstFound = 1;
                end
            end

            [minor_crit, k_crit2] = min(minors);
            k_crit = [k_crit1, k_crit2];

        end


    end
end

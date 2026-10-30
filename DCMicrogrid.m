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
        DesignOptions = struct()
        DesignCertificate = struct()
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

                assert(i>=1 && i<=N && j>=1 && j<=N && i~=j);
                assert(isfinite(g) && g>0);
                Y(i,j) = Y(i,j)+g;
                Y(j,i) = Y(j,i)+g;
                
                YBar(i,j) = YBar(i,j)-g;
                YBar(j,i) = YBar(j,i)-g;

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

            assert(~useData,'DCMG:LegacyAdaptiveContract', ...
                'For certified integral DD results use run_paper_results; legacy adaptive stages use a different sampling contract.');
            ode_opts = odeset('RelTol',1e-6,'AbsTol',1e-6);
            
            dt = 1e-5;
            t_0 = tspan(1);
            t_f = tspan(end);

            t_1 = t_0+0.05*(t_f-t_0); % When the data-driven method is triggered
            t_2 = t_0+0.2*(t_f-t_0); % When a high disturbance enters
            t_3 = t_0+0.4*(t_f-t_0); % When a state disturbance enters
            t_4 = t_0+0.6*(t_f-t_0); % When a parameter disturbance enters
            t_5 = t_0+0.65*(t_f-t_0); % When data-driven method reacts
            t_6 = t_0+0.8*(t_f-t_0); % When high disturbance and a state disturbance enters together

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
                        X(end, 2*i-1:2*i).*(mag*(2*rand(1,2)-1));

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
                obj.DGs(i).L=obj.DGs(i).L*(1+mag*(2*rand()-1));
                obj.DGs(i).C=obj.DGs(i).C*(1+mag*(2*rand()-1));
                obj.DGs(i).Rf=obj.DGs(i).Rf*(1+mag*(2*rand()-1));
                obj.DGs(i).RL=obj.DGs(i).RL*(1+mag*(2*rand()-1));
                obj.DGs(i).updateModel();
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
                        X(end, 2*i-1:2*i).*(mag*(2*rand(1,2)-1));

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
        
            ode_opts = odeset('RelTol',1e-9,'AbsTol',1e-10);
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
                    grid=unique([t0:dt:tk1,tk1]);
                    [t, X] = ode45(f, grid, Xk, ode_opts);
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

        function out = loadDataMatrices(obj,data,varargin)
            % Ingest the declared exact integral-data contract. Historical
            % interpolated forward-Euler samples do not meet this contract.
            assert(isstruct(data) && isfield(data,'Local') && isempty(varargin), ...
                'DCMG:LegacyDataContract', ...
                'Use data=collectIntegralData(net); net.loadDataMatrices(data).');
            assert(numel(data.Local)==obj.N,'DCMG:DataDimensions');
            for i=1:obj.N
                entry=data.Local{i};
                % Validation truth fields are deliberately withheld from synthesis.
                fields={'EquationResidual','TrueModelMargin','BoundSlack'};
                for j=1:numel(fields)
                    if isfield(entry,fields{j}), entry=rmfield(entry,fields{j}); end
                end
                obj.DGs(i).DataContract=entry;
            end
            out=struct('Kind','continuous-integral','Samples',data.Options.Samples, ...
                'SampleTime',data.Options.SampleTime,'Accepted',true);
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
            Phi = [I_n; W'];    % size: (nW + T) x nW
            
            % QMI constraint
            M = Phi' * Q_w * Phi;    % size: nW x nW
            
            Constraints = [M >= epsilon*eye(size(M)), epsilon >= 1e-3, Q_I >= epsilon*eye(size(Q_I))];
            
            % To avoid the trivial Q_w = 0, add some objective.
            % Example: minimize Frobenius norm of Q_w, or add structure.
            Objective = trace(Q_I) + epsilon;   % or something else
            
            opts = dcmgSolverOptions();
            sol  = optimize(Constraints, Objective, opts);
            
            if sol.problem ~= 0
                error('SDP infeasible or solver failed: %s', sol.info);
            end
            
            Qw_val = value(Q_w);

            % REVIEW PROPOSAL DATA-04 (ADD): this fit has only a fixed 1e-3
            % post-hoc reserve. Prefer a physical bound selected before data.
            % If retained, scale its reserve to normalized disturbance units.
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
                    if DG.EnableLimiters, Iline_i = min(Imax, max(Imin, Iline_i)); end
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
                    if DG.EnableLimiters, u_i = min(Vmax, max(Vmin, u_i)); end
                    % REVIEW PROPOSAL DATA-06 (ADD): disable command clipping
                    % for the certified linear run, or model/log saturation.
                    % u_i = u_Si + u_Li + u_Gi;
        
                    ITilde_ki = Iline_i - I_S(i);
                    uTilde_ki = u_i - u_Si;
                    
                    tVal = t(k);
                    if isempty(DG.noise) || ~isfield(DG.noise,'t') || isempty(DG.noise.t)
                        w_i=zeros(2,1);
                    else
                        w_i=interp1(DG.noise.t,DG.noise.w,tVal,'previous','extrap')';
                    end
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


        function ss = solveSteadyState(obj,opts)
            if nargin<2, opts=struct(); end
            ss=designSteadyState(obj,opts);
        end

        function [AdjMat,KMat,out] = design_MB_GSC(obj,linkFiltThresh,useData)
            % Existing dense Lyapunov design, with physical rows eliminated.
            % Additive excitation is handled by collectIntegralData.
            if nargin<2, linkFiltThresh=0; end
            if nargin<3, useData=false; end %#ok<NASGU>
            n=2*obj.N; tau=1e-3;
            A=tau*(obj.A+obj.E*obj.YBar*obj.D'); B=tau*obj.B;
            P=sdpvar(n,n,'symmetric'); H=sdpvar(obj.N,n,'full');
            W=-(A*P+B*H)-(A*P+B*H)';
            cons=[P>=1e-3*eye(n),W>=1e-3*eye(n)];
            ops=dcmgSolverOptions();
            sol=optimize(cons,trace(P)+1e-3*norm(H(:),1),ops);
            out=struct('problem',sol.problem,'info',sol.info,'Solver',ops.solver);
            AdjMat=[]; KMat=[];
            if sol.problem~=0, return; end
            out.P=value(P); out.H=value(H); out.LMI=value(W);
            gain=out.H/out.P; out.UnprunedGain=gain;
            out.LMIDiagnostics=certificateDiagnostics(out.LMI);
            candidate=gain;
            for i=1:obj.N
                for j=1:obj.N
                    ix=2*j-1:2*j;
                    if norm(candidate(i,ix))<linkFiltThresh*max(abs(gain(:)))
                        candidate(i,ix)=0;
                    end
                end
            end
            % Retain a proposed pruning only when the original storage checks.
            storage=inv(out.P); Ac=A+B*candidate;
            check=certificateDiagnostics(-Ac'*storage-storage*Ac);
            if ~check.Accepted, candidate=gain; Ac=A+B*candidate; end
            out.OriginalResidual=-Ac'*storage-storage*Ac;
            out.OriginalDiagnostics=certificateDiagnostics(out.OriginalResidual);
            out.MaxRealPole=max(real(eig(Ac)))/tau;
            out.Accepted=out.LMIDiagnostics.Accepted && ...
                out.OriginalDiagnostics.Accepted && out.MaxRealPole<0;
            if ~out.Accepted
                out.problem=4; out.info='Stabilizing certificate rejected'; return
            end
            obj.K=candidate;
            for i=1:obj.N, obj.DGs(i).K=[0,0]; end
            [AdjMat,KMat]=obj.buildCommAdjFromK(0);
            out.CommandGain=KMat;
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
                    if norm(K_ij,2) > max(thr*maxK,1e-12)
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



       
        

        function [AdjMat,KMat,out] = codesign_MB_DRC(obj,linkFiltThresh)
            if nargin<2, linkFiltThresh=0; end
            [AdjMat,KMat,out]=obj.codesignCertified(false,linkFiltThresh);
        end

        function [AdjMat,KMat,out] = codesign_DD_DRC(obj,linkFiltThresh)
            % Integral-data local synthesis followed by the known-B global LMI.
            if nargin<2, linkFiltThresh=0; end
            [AdjMat,KMat,out]=obj.codesignCertified(true,linkFiltThresh);
        end

        function [AdjMat,KMat,out] = codesignCertified(obj,useData,linkFiltThresh)
            previous=captureControllerState(obj);
            AdjMat=[]; KMat=[]; local=cell(obj.N,1);
            try
            for i=1:obj.N
                obj.DGs(i).DesignOptions=obj.DesignOptions;
                if useData
                    oi=obj.DGs(i).designLocalXiDissipative_DataDriven();
                else
                    oi=obj.DGs(i).designLocalXiDissipative();
                end
                if oi.problem~=0
                    out=oi; out.FailedDG=i; restoreControllerState(obj,previous); return
                end
                oi.TruthVerification=verifyLocalDissipativity(obj.DGs(i),oi);
                assert(oi.TruthVerification.Accepted,'DCMG:LocalVerification', ...
                    'Original local storage inequality failed at DG %d.',i);
                local{i}=oi;
            end
            out=designGlobalDissipativity(obj);
            out.Local=local;
            if out.problem~=0, restoreControllerState(obj,previous); return; end
            % Threshold proposes a graph only. A second hard solve certifies it.
            if linkFiltThresh>0
                mask=true(obj.N);
                for i=1:obj.N
                    for j=1:obj.N
                        if i~=j
                            mask(i,j)=norm(out.K(i,2*j-1:2*j))>= ...
                                linkFiltThresh*max(abs(out.K(:)));
                        end
                    end
                end
                previousOptions=obj.DesignOptions;
                obj.DesignOptions.AllowedEdges=mask;
                sparse=designGlobalDissipativity(obj);
                obj.DesignOptions=previousOptions;
                out.SparsityAttempt=sparse;
                if sparse.problem==0
                    dense=out; out=sparse; out.DenseDesign=dense; out.Local=local;
                end
            end
            AdjMat=out.Adjacency; KMat=obj.K;
            obj.DesignCertificate=out;
            catch exception
                restoreControllerState(obj,previous);
                rethrow(exception);
            end
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

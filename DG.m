classdef DG < handle
    % DG  Inverter-based DC source with L-C filter and local load (current sink)
    % State order: x = [vC; iL]
    % Model: xdot = A*x + B*u + E*w + sum_j Ytilde_ij * x_j
    properties
        
        id              (1,1) double

        % Physical params
        L               (1,1) double = 1e-3   % H
        C               (1,1) double = 1e-3   % F
        Rf              (1,1) double = 0.1    % Ohm
        YBar 

        % Operating/loads
        RL      (1,1) double = 12.0;   % constant-impedance load [Ohm]
        Ibar    (1,1) double = 0.0;    % constant-current part [A]  (== wbar_i)
        iload_fun             = @(t) 0; % small time-varying perturbation (A)
        
        % Control parameters
        u_s            (1,1) double = 48     % steady input command (V)
        K               (1,2) double = [0 0]  % neighbor-state feedback row (placeholder)
        
        % State
        x               (2,1) double = [48; 0]  % [vC; iL]
        noise

        % Drawing
        pos             (1,2) double = [0 0]
        color           (1,3) double = [0 0.45 0.74]
        
        % Model matrices (filled by updateModel)
        A   (2,2) double
        B   (2,1) double
        E   (2,1) double
        BBar
        % gsum (1,1) double = 0.0 % sum_j g_ij (set by network)

        % Ratings
        Vrated  (1,1) double = 48;     % [V]
        Irated  (1,1) double = 3.5;    % [A] (we'll set better below)

        %steady state 
        x_s = [48; 3];
        I_s = 0.5;

        nu
        rho

        xTilde
        xTildeBar
        yTilde
        uTilde
        wTilde
        Q_w
        QBar_w

        Kc          % continuous global controller gain 
        Kd          % discrete-time global gain
        Ts     % controller sampling period used in data-driven design
        u_L_Hold

    end

    methods

        function obj = DG(id, params)
            if nargin >= 1, obj.id = id; else, obj.id = 1; end
            if nargin >= 2 && ~isempty(params)
                f = fieldnames(params);
                for k = 1:numel(f)
                    obj.(f{k}) = params.(f{k});
                end
            end
            obj.updateModel(); % initialize A,B,E
        end

        function updateModel(obj)

            % A,B,E for state order [vC; iL]
            obj.A = [ -(1/(obj.RL*obj.C))       1/obj.C ; ...
                        -1/obj.L             -obj.Rf/obj.L ];
            obj.B = [0; 1/obj.L];
            obj.E = [-1/obj.C; 0];   % multiplies any current sink (Ibar, Iline, perturbation)
            obj.BBar = [obj.E, obj.B];

        end

        function setState(obj, xnew), obj.x = xnew(:); end
        
        function x = getState(obj),   x = obj.x;       end

        function dx = dynamics(obj, t, Iline_i, u_G, useData)
            
            % Piecewise-constant noise (sample & hold)
            w_i = interp1(obj.noise.t, obj.noise.w, t, 'previous', 'extrap')';
        
            if useData
                u_L = obj.u_L_Hold;
            else
                u_L = obj.K*(obj.x - obj.x_s);
            end
            
            u_i = obj.u_s + u_L + u_G;                     % controller-computed input

            % saturated_x = min(upper_bound, max(lower_bound, x));
            % Saturation
            Imax = 2*obj.Irated; Imin = -2*obj.Irated; 
            Iline_i = min(Imax, max(Imin, Iline_i));
            Vmax = 2*obj.Vrated; Vmin = -2*obj.Vrated;
            u_i = min(Vmax, max(Vmin, u_i));
            
        
            dx  = obj.A*obj.x + obj.BBar*[Iline_i; u_i] + obj.E*obj.Ibar + obj.BBar*w_i;

        end


        function draw(obj, ax)

            if nargin < 2 || isempty(ax), ax = gca; end

            p = obj.pos; hold(ax, 'on');

            % Node marker
            plot(ax, p(1), p(2), 'o', 'MarkerSize', 8, ...
                'MarkerFaceColor', obj.color, 'MarkerEdgeColor','k');
            % Little pictogram and labels
            % line(ax, [p(1)-0.5 p(1)], [p(2) p(2)], 'Color', obj.color, 'LineWidth', 1.5);
            text(p(1)+0.1, p(2)+0.2, sprintf('DG %d', obj.id), 'Color', obj.color);
            text(p(1)+0.1, p(2)+0, sprintf('V=%.1fV (%.1f pu)', obj.x(1), obj.x(1)/obj.Vrated), ...
                 'Color',[0.25 0.25 0.25], 'FontSize',8);
            text(p(1)+0.1, p(2)-0.2, sprintf('It=%.1fA (%.1f%%)', obj.x(2), 100*obj.x(2)/obj.Irated), ...
                 'Color',[0.25 0.25 0.25], 'FontSize',8);
            text(p(1)+0.1, p(2)-0.4, sprintf('I=%.1fA', obj.I_s(1)), ...
                 'Color',[0.25 0.25 0.25], 'FontSize',8);

        end

        % function C = outputSelector(obj)
        %     % Local regulated output y_i; default = bus voltage error
        %     C = [1 0];
        % end
        
        function out = designLocalXiDissipative(obj)
            
            A = obj.A;  
            BBar = obj.BBar;
            YBar = obj.YBar;
            D = [1; 0];
            DBar = [0; 1];
        
            P   = sdpvar(2,2,'symmetric');
            L   = [zeros(1,2); sdpvar(1,2,'full')];       % Y = K_iL * P
            xBar11 = sdpvar(1,1);
            x22 = sdpvar(1,1);
            x12 = 0.5;
            x21 = x12';

            epsilon = sdpvar(1,1);

            I = eye(2);
            O = zeros(2);
            O_12 = zeros(1,2);

            % For Necessary Conditions
            KTilde = sdpvar(2,2,'full');
            KHat = sdpvar(2,2,'full');
            yBar11 = sdpvar(1,1,'full');
            yBar12 = 0;
            yBar21 = 0;
            yBar22 = sdpvar(1,1,'full');
        
            % Main LMI
            AP_BL = A*P + BBar*L;                % Acl_i * P
            LMI1 = [I,    P,                  O;
                   P,  -AP_BL'-AP_BL,      x22*I + x21*P;
                   O,   x22*I + x12*P       xBar11*I];

            BK = BBar*KHat;
            LMI2 = [xBar11*I,           O,              BBar*KTilde,      xBar11*I;
                    O,                  -yBar22*I     -yBar22*I,          O;
                    KTilde'*BBar'    -yBar22*I       -BK-BK'+I,          -x21*I;
                    xBar11*I,           O,              -x12*I,             yBar11*I];
            
            scalarCons = [x22 <= -epsilon, yBar22 <= -epsilon, epsilon >= 0.001];
            
            % phyCons = [1e12*D'*L == O_12];
            
            
            mainCons = [P >= epsilon*eye(2), LMI1 >= epsilon*eye(size(LMI1))];
            mainCons = [mainCons, LMI2 >= epsilon*eye(size(LMI2))];
            
            cons = [scalarCons, mainCons];

            % BMI Cons
            xBar11ValGuess = 0.002; % revise this
            BMICons = [xBar11ValGuess*KHat == x21*KTilde];
            cons = [cons, BMICons];

            
            % Optimize
            ops  = sdpsettings('solver','mosek','verbose',0);
            
            objective = 1*yBar11 - 1*yBar22 + trace(P) + epsilon;

            sol  = optimize(cons, objective, ops);
        
            out.problem = sol.problem;
            out.info    = yalmiperror(sol.problem);
            out.P       = value(P);
            out.L       = value(L);
            
            LVal = value(L);
            PVal = value(P);

            xBar11Val = value(xBar11)
            if abs(xBar11Val-xBar11ValGuess) > 1e-6
                disp(['Revise xBar11ValGuess towards: ',num2str(xBar11Val), ' from ', num2str(xBar11ValGuess),'.'])
            end
            x22Val = value(x22);
            x11Val = (-x22Val)\xBar11Val;

            out.nu      = -x11Val;
            out.rho     = -x22Val;

            [k_crit, minor_crit, ~] = obj.criticalLeadingMinor(value(LMI1));
            [k_crit2, minor_crit2, ~] = obj.criticalLeadingMinor(value(LMI2));
            if abs(minor_crit)>1e-3 || abs(minor_crit2)>1e-3
                disp('Error in Minors!')
            end

            if sol.problem==0
                
                K = out.L / out.P;
                out.K = DBar'*K;

                eigs = eig(obj.A + obj.BBar*K)

                % Store results in the DG
                obj.K       = out.K;
                obj.nu      = out.nu;
                obj.rho     = out.rho;

                disp(['Local Design Success at DG ',num2str(obj.id),'!']);
                disp(['(nu,rho)=(',num2str(out.nu),',',num2str(out.rho),'), ','K=[',num2str(out.K(1)),',',num2str(out.K(2)),']']);

            else
                warning('Local Design Fail: %s', sol.info);
                out.K = [0 0];

                obj.K       = out.K;
                obj.nu      = 0;
                obj.rho     = 0;
            end

            
        end


        function out = designLocalXiDissipative_DataDriven(obj)
            % Data–driven local Xi–dissipative controller design (Proposition 7)
            
            % Basic dimensions from the DG model
            n = size(obj.A,1);        % state dimension
            m = size(obj.BBar,2);     % input dimension
            r = n;                    % output dimension (in this case n = m = r)
            T = size(obj.xTilde,2);   % number of columns / samples

            % eps0 = 1e-6;
            % Sx = diag( 1./(max(abs(obj.xTilde),[],2) + eps0) )
            % Su = diag( 1./(max(abs(obj.uTilde),[],2) + eps0) )

            % Data matrices
            X = obj.xTilde;
            % U = obj.uTilde;
            Y = obj.yTilde;

            QBar_w  = 0.5*(obj.QBar_w + obj.QBar_w');
            scaleQ = max(1, max(abs(QBar_w(:))));   % e.g. ~1e4
            QBar_w = QBar_w / scaleQ;

            % Basic matrices
            I_n  = eye(n);
            I_m  = eye(m);
            I_r = eye(r);
            O_n = zeros(n); O_m = zeros(m); O_r = zeros(r);
            O_nm = zeros(n,m); O_mn = zeros(m,n); 
            O_rn = zeros(r,n); O_rm = zeros(r,m); 
            O_nr = zeros(n,r); O_mr = zeros(m,r);
            D = [1; 0];
            DBar = [0; 1];
      
            % ---------------------------------------------------------------------
            % Decision variables of Proposition 7
            % ---------------------------------------------------------------------
            % Main controller matrices
            
            KTilde = [zeros(1,2); sdpvar(1,2,'full')]; % \tilde K_i
            KHat   = sdpvar(T,n,'full');   % \hat K_i
        
            % Scalar multipliers for the QMI
            lambda1 = sdpvar(1,1);
            lambda2 = sdpvar(1,1);
        
            % Xi-dissipativity scalars {x̄_i^{11}, x_i^{12}, x_i^{21}, x_i^{22}}
            xBar11  = sdpvar(1,1); %=-nu*rho
            x12     = 0.5;
            x21     = 0.5;
            x22     = sdpvar(1,1);

            % For the "necessary condition" LMI (same set as in Prop. 5 code)
            KTildeii = sdpvar(m,n,'full');
            KHatii   = sdpvar(m,n,'full');
            yBar11 = sdpvar(1,1,'full');
            yBar12 = 0;
            yBar21 = 0;
            yBar22 = sdpvar(1,1,'full');
        
            % Tuning / regularization parameter
            epsilon = sdpvar(1,1);
        
            % ---------------------------------------------------------------------
            % LMI #1:  [Qi1  Si1;  Si1'  Ri1] >= 0  (data–driven robust LMI)
            % ---------------------------------------------------------------------
            
            % Q_{i1}
            XKhatSym = 0.5*(X*KHat + (X*KHat)');
            Qi1 = [ I_r,            Y*KHat,                     O_r;
                    (Y*KHat)',     XKhatSym,                x21*KHat'*Y';
                    O_r,          (x21*KHat'*Y')',              xBar11*I_r ];
        
            % S_{i1}
            Si1 = [ O_rn,            O_rn,           O_rm;
                    O_n,           XKhatSym,          KTilde';
                    (-x22)*I_r,     O_rn,           O_rm ];    % Here the fact that I_r = I_n helps
        
            % R_{i1}
            % First build the "structured" part, then subtract lambda1*Qw_bar
            R_struct_i1 = [ XKhatSym,     O_n,     O_nm;
                            O_n,        O_n,     O_nm;
                            O_mn,       O_mn,    O_m];
        
            Ri1 = R_struct_i1 - lambda1*QBar_w;
        

            LMI1 = [Qi1,  Si1;
                    Si1', Ri1];
        
            % ---------------------------------------------------------------------
            % LMI #2:  [Qi2  Si2;  Si2'  Ri2] >= 0  (necessary conditions, data–driven)
            % ---------------------------------------------------------------------
            
            % Q_{i2}
            Qi2 = [ -yBar22*I_n,      O_n;
                    O_n,               yBar11*I_n ];
        
            % S_{i2}
            % (block structure follows directly from the screenshot; dimensions are
            % multiples of n; adapt if your Xi-parameterization differs)
            Si2 = [ O_n,           O_n,  O_nm, -yBar22*I_n,          O_n,        O_nm;
                    xBar11*I_n,    O_n,  O_nm, (-x21 + yBar21)*I_n,    O_n,        O_nm];
        
            % R_{i2} – large block matrix then minus lambda2*blkdiag(Qw_bar,Qw_bar)
            R_struct_i2 = [ xBar11*I_n,  O_n,        O_nm,        O_n,            O_n,       O_nm;
                            O_n,         O_n,        O_nm,        O_n,            O_n,       O_nm;
                            O_mn,        O_mn,        O_m,        KTildeii,        O_mn,      O_m;
                            O_n,         O_n,        KTildeii',   I_n,            O_n,      -KHatii';
                            O_n,         O_n,        O_nm,        O_n,            O_n,       O_nm;
                            O_mn,         O_mn,      O_m,         -KHatii,       O_mn,       O_m ];
        
            Ri2 = R_struct_i2 - lambda2*blkdiag(QBar_w, QBar_w);
        
            LMI2 = [ Qi2,  Si2;
                     Si2', Ri2 ];
        
            % ---------------------------------------------------------------------
            % Collect constraints
            % ---------------------------------------------------------------------
            % Sign constraints from Proposition 7
            
            scalarCons = [x22 <= -epsilon, yBar22 <= -epsilon,...
                          lambda1 >= 0, lambda2 >= 0, ...
                          epsilon >= 1e-9]; %%%% check

            % Physical Constraints:
            % O_12 = zeros(1,2);
            % phyCons = [1e12*D'*KTilde == O_12]; %%%% check

            % Symmetry contraitns
            symCons = [X*KHat == (X*KHat)'];

            % Main LMIs
            mainCons = [XKhatSym >= epsilon*eye(size(XKhatSym))];
            mainCons = [mainCons, LMI1 >= epsilon*eye(size(LMI1))];
            mainCons = [mainCons, LMI2 >= epsilon*eye(size(LMI2))];
        
            cons = [scalarCons, symCons, mainCons];


            % BMI Cons
            xBar11ValGuess = 1e-8; % revise this
            BMICons = [xBar11ValGuess*KHatii == x21*KTildeii];
            % cons = [cons, BMICons];


            % cons = [cons, -0>=x22, x22>=-0.5, 2>=xBar11, xBar11>=0.1];
            % ---------------------------------------------------------------------
            % Solve the SDP
            % ---------------------------------------------------------------------
            ops = sdpsettings('solver','mosek','verbose',0,'showprogress',0);
        
            % Simple objective: favor "small" Xi, y-variables and multipliers
            objective = 1*yBar11 - 1*yBar22 + 1*trace(XKhatSym) + 1*epsilon; %%%% check
            % tempCost  = (x22-(-0.296))^2+(xBar11-0.9584)^2;
            sol = optimize(cons, objective, ops);
        
            % ---------------------------------------------------------------------
            % Pack results
            % ---------------------------------------------------------------------
            out.problem = sol.problem;
            out.info    = yalmiperror(sol.problem);
            
            out.P = value(XKhatSym);
            out.L = value(KTilde);
            
            LVal = value(KTilde);
            PVal = value(XKhatSym);

            % XBar = obj.xTildeBar;
            % XXT = XBar*XBar'
            % W = obj.wTilde;
            % WWT = W*W'
            Q_w_disp = obj.Q_w(1:5,1:5);
            QBar_w_disp = QBar_w(1:5,1:5);  

            lam1 = value(lambda1);
            lam2 = value(lambda2);
            
            Qi1Val = value(Qi1);
            Si1Val = value(Si1);
            R_struct_i1Val = value(R_struct_i1);
            Ri1Val = value(Ri1);

            % KTildeVal = value(KTilde)
            % LVal_2 = DBar'*LVal
            % Y_0 = D'*KTilde
            yBar11Val = value(yBar11);
            yBar22Val = value(yBar22);
            % KTYT = value(x21*KHat'*Y')
            % value(X*KHat)/value(XKhatSym)

            xBar11Val = value(xBar11)
            if abs(xBar11Val-xBar11ValGuess) > 1e-6
                disp(['Revise xBar11ValGuess towards: ',num2str(xBar11Val), ' from ', num2str(xBar11ValGuess),'.'])
            end
            x22Val = value(x22);
            x11Val = (-x22Val)\xBar11Val;
            out.nu = -x11Val;
            out.rho = -x22Val;
        
            [k_crit, minor_crit, ~] = obj.criticalLeadingMinor(value(LMI1));
            [k_crit2, minor_crit2, ~] = obj.criticalLeadingMinor(value(LMI2));
            if abs(minor_crit)>1e-3 || abs(minor_crit2)>1e-3
                disp('Error in Minors!')
            end

            if sol.problem == 0
                
                K = out.L / out.P;
                out.K = DBar'*K;

                A = (eye(2) + obj.Ts*obj.A);
                B = obj.Ts*obj.BBar;
                eigs = eig(A + B*K);

                eigs = abs(eigs)
                % K = place(A,-B,real(eigs))
                % eigs = eig(A + B*K)
                % out.K = DBar'*K;

                % Store results in the DG
                obj.K       = out.K;
                obj.nu      = out.nu;
                obj.rho     = out.rho;

                disp(['Local Design Success at DG ',num2str(obj.id),'!']);
                disp(['(nu,rho)=(',num2str(out.nu),',',num2str(out.rho),'), ','K=[',num2str(out.K(1)),',',num2str(out.K(2)),']']);

            else

                warning('Local Design Fail: %s', sol.info);
                out.K = [0 0]; 

                obj.K       = out.K;
                obj.nu      = 0;
                obj.rho     = 0;

            end

        end

        
        function [k_crit, minor_crit, minors] = criticalLeadingMinor(obj, M)
            
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

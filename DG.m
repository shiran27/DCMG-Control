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
            obj.A = [ -(1/(obj.RL*obj.C))   1/obj.C ; ...
                  -1/obj.L             -obj.Rf/obj.L ];
            obj.B = [0; 1/obj.L];
            obj.E = [-1/obj.C; 0];   % multiplies any current sink (Ibar, Iline, perturbation)
            obj.BBar = [obj.E, obj.B];

        end

        function setState(obj, xnew), obj.x = xnew(:); end
        function x = getState(obj),   x = obj.x;       end

        function dx = f(obj, t, Iline, u_G)
            
            iload_dist = 1*sin(1*t);

            u_L = obj.K*(obj.x - obj.x_s);
            
            u_i = obj.u_s + u_L + u_G;                     % controller-computed input
        
            dx  = obj.A*obj.x + obj.BBar*[Iline; u_i] + obj.E*(obj.Ibar + iload_dist);

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
        
        function out = designLocalXiDissipative(obj, opts)
            

            A = obj.A;  
            BBar = obj.BBar;
            YBar = obj.YBar;
            D = [1; 0];
            DBar = [0; 1];
        
            P   = sdpvar(2,2,'symmetric');
            L   = sdpvar(2,2,'full');       % Y = K_iL * P
            xBar11 = sdpvar(1,1);
            x22 = sdpvar(1,1);
            x12 = 0.5;
            x21 = x12';

            epsilon = sdpvar(1,1);

            I = eye(2);
            O = zeros(2);

            % For Necessary Conditions
            KTilde11 = sdpvar(2,2,'full');
            KHat21 = sdpvar(2,2,'full');
            yTilde11 = sdpvar(1,1,'full');
            yHat22 = sdpvar(1,1,'full');
            yCheck22 = sdpvar(1,1,'full');
            yBar12 = 0;
            yBar21 = 0;
        
            % Main LMI
            AP_BL = A*P + BBar*L;                % Acl_i * P
            LMI = [I,    P,                  O;
                   P,  -AP_BL'-AP_BL,      x22*I + P*x21;
                   O,   x22*I + P*x12       xBar11*I];
        
            cons = [P >= epsilon*eye(2), epsilon >= 0.001, LMI >= epsilon*eye(size(LMI)), D'*L == 0*YBar*D'*P, x22 <= -epsilon];

           BK = BBar*KHat21;
            LMI2 = [xBar11*I,           O,              BBar*KTilde11,      xBar11*I;
                    O,                  -yCheck22*I     -yHat22*I,          O;
                    KTilde11'*BBar'    -yHat22*I       -BK-BK'+I,          -x21*I;
                    xBar11*I,           O,              -x21*I,             yTilde11*I];
            
            cons = [cons, LMI2 >= epsilon*eye(size(LMI2)), yHat22 <= -epsilon];
            
            % Optimize


            ops  = sdpsettings('solver','mosek','verbose',0);
            sol  = optimize(cons, 1*yTilde11 + trace(P) + epsilon, ops);
        
            out.problem = sol.problem;
            out.info    = yalmiperror(sol.problem);
            out.P       = value(P);
            out.L       = value(L);
            xBar11Val = value(xBar11); 
            x22Val = value(x22); 
            x11Val = (-x22Val)\xBar11Val;
            out.nu      = -x11Val;
            out.rho     = -x22Val;

            if sol.problem==0
                
                K = out.L / out.P;
                out.K = DBar'*K;
                disp(['Local Design Success! (nu,rho)=(',num2str(out.nu),',',num2str(out.rho),')','K=[',num2str(K(1)),',',num2str(K(2)),']'])

            else
                warning('Local Design Fail: %s', sol.info);
                out.K = [0 0];
            end

            % Store results in the DG
            obj.K       = out.K;
            obj.nu      = out.nu;
            obj.rho     = out.rho;
        end

    end
end

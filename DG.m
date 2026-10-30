classdef DG < handle
    % DG Inverter-based DC source, local dynamics, and local design routines.
    %
    % Runtime role: key simulation component used by main.mlx through
    % DCMicrogrid. The state order is x = [vC; iL]. The network supplies
    % the net line current, while the DG applies its local/global feedback
    % command and additive disturbance channels.
    %
    % Design role: implements the model-based and data-driven local
    % dissipativity controller designs used by the distributed co-design.
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
        EnableLimiters (1,1) logical = false % certified linear mode by default
        DesignOptions = struct()
        DesignCertificate = struct()
        DataContract = struct()

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

            % REVIEW PROPOSAL DG-01 (ADD): fail early on nonphysical data.
            % validateattributes(obj.L,  {'numeric'}, {'scalar','real','positive','finite'});
            % validateattributes(obj.C,  {'numeric'}, {'scalar','real','positive','finite'});
            % validateattributes(obj.Rf, {'numeric'}, {'scalar','real','nonnegative','finite'});
            % validateattributes(obj.RL, {'numeric'}, {'scalar','real','positive','finite'});
            % validateattributes(obj.Ibar, {'numeric'}, {'scalar','real','finite'});
            % assert(isequal(size(obj.K), [1 2]), 'DG:InvalidGain', ...
            %     'The implemented local voltage-command gain must be 1-by-2.');

            validateattributes(obj.L,{'numeric'},{'scalar','positive','finite'});
            validateattributes(obj.C,{'numeric'},{'scalar','positive','finite'});
            validateattributes(obj.Rf,{'numeric'},{'scalar','nonnegative','finite'});
            validateattributes(obj.RL,{'numeric'},{'scalar','positive','finite'});
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
            % REVIEW PROPOSAL DG-02 (REPLACE): use this guarded block in
            % place of the next interp1 line if zero-noise runs are needed.
            % if isempty(obj.noise) || ~isfield(obj.noise,'t') || isempty(obj.noise.t)
            %     w_i = zeros(2,1);
            % else
            %     w_i = interp1(obj.noise.t,obj.noise.w,t,'previous','extrap')';
            % end
            if isempty(obj.noise) || ~isfield(obj.noise,'t') || isempty(obj.noise.t)
                w_i=zeros(2,1);
            else
                w_i = interp1(obj.noise.t,obj.noise.w,t,'previous','extrap')';
            end
        
            if useData
                u_L = obj.u_L_Hold;
            else
                u_L = obj.K*(obj.x - obj.x_s);
            end
            
            u_i = obj.u_s + u_L + u_G;                     % controller-computed input

            % saturated_x = min(upper_bound, max(lower_bound, x));
            % Saturation
            % REVIEW PROPOSAL DG-03 (ADD, part 1 of 2): save the network
            % current before the active clipping line.
            % Iline_i_network = Iline_i;
            Imax = 2*obj.Irated; Imin = -2*obj.Irated; 
            if obj.EnableLimiters, Iline_i = min(Imax, max(Imin, Iline_i)); end
            % REVIEW PROPOSAL DG-03 (ADD, part 2 of 2): restore the network
            % current after the active clipping line. Recommended for the
            % linear/theory-consistent run; a real limiter needs its own model.
            % Iline_i = Iline_i_network;

            Vmax = 2*obj.Vrated; Vmin = -2*obj.Vrated;
            if obj.EnableLimiters, u_i = min(Vmax, max(Vmin, u_i)); end
            % REVIEW PROPOSAL DG-04 (DESIGN): command saturation is outside
            % the LMI certificate. Log every activation, or disable it for
            % certificate-validation runs and report it as a nonlinear stress test.
            
        
            dx  = obj.A*obj.x + obj.BBar*[Iline_i; u_i] + obj.E*obj.Ibar + obj.BBar*w_i;
            % REVIEW PROPOSAL DG-05 (REPLACE): include the configured load
            % ripple if iload_fun is intended to be an active disturbance.
            % dx = obj.A*obj.x + obj.BBar*[Iline_i;u_i] + ...
            %      obj.E*(obj.Ibar + obj.iload_fun(t)) + obj.BBar*w_i;

        end


        function draw(obj, ax, varargin)

            if nargin < 2 || isempty(ax), ax = gca; end

            parser = inputParser;
            parser.FunctionName = 'DG.draw';
            addParameter(parser,'NetworkCurrent',obj.I_s(1), ...
                @(value)isnumeric(value) && isscalar(value));
            addParameter(parser,'ShowLabels',true, ...
                @(value)islogical(value) && isscalar(value));
            addParameter(parser,'FontSize',8, ...
                @(value)isnumeric(value) && isscalar(value) && value>0);
            addParameter(parser,'MarkerSize',9, ...
                @(value)isnumeric(value) && isscalar(value) && value>0);
            addParameter(parser,'LabelPosition',obj.pos+[0.1 0.05], ...
                @(value)isnumeric(value) && isequal(size(value),[1 2]));
            addParameter(parser,'HorizontalAlignment','left', ...
                @(value)ischar(value) || (isstring(value) && isscalar(value)));
            addParameter(parser,'VerticalAlignment','middle', ...
                @(value)ischar(value) || (isstring(value) && isscalar(value)));
            parse(parser,varargin{:});
            options = parser.Results;

            p = obj.pos; hold(ax,'on');

            plot(ax,p(1),p(2),'o','MarkerSize',options.MarkerSize, ...
                'MarkerFaceColor',obj.color,'MarkerEdgeColor',[0.1 0.1 0.1], ...
                'LineWidth',0.9);

            if ~options.ShowLabels
                return;
            end

            labelColor = [0.18 0.21 0.24];
            label = sprintf(['DG %d\nV=%.2f V (%.3f pu)\n', ...
                'It=%.2f A (%.1f%%)\nInet=%+.2f A'], ...
                obj.id,obj.x(1),obj.x(1)/obj.Vrated,obj.x(2), ...
                100*obj.x(2)/obj.Irated,options.NetworkCurrent);
            text(ax,options.LabelPosition(1),options.LabelPosition(2),label, ...
                'Color',labelColor,'FontSize',options.FontSize, ...
                'FontWeight','normal','Clipping','off', ...
                'HorizontalAlignment',char(options.HorizontalAlignment), ...
                'VerticalAlignment',char(options.VerticalAlignment), ...
                'BackgroundColor','w','Margin',1,'Interpreter','none');

        end

        % function C = outputSelector(obj)
        %     % Local regulated output y_i; default = bus voltage error
        %     C = [1 0];
        % end
        
        function out = designLocalXiDissipative(obj)
            out=designLocalIntegralDissipativity(obj,false);
        end

        function out = designLocalXiDissipative_DataDriven(obj)
            out=designLocalIntegralDissipativity(obj,true);
        end

        function [k_crit, minor_crit, minors] = criticalLeadingMinor(obj, M)

            % REVIEW PROPOSAL NUM-02 (REPLACE): use the minimum eigenvalue of
            % the symmetrized LMI instead of determinants of leading minors.
            % minEigenvalue = min(eig((value(M)+value(M)')/2));
            % assert(minEigenvalue >= -1e-7, 'DG:LMIViolation', ...
            %     'Recovered local LMI violates positive semidefiniteness.');
            
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

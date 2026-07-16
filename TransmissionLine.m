classdef TransmissionLine < handle
    % TransmissionLine Resistive physical connection between two DG buses.
    %
    % Runtime role: key simulation component used by DCMicrogrid to build
    % the nodal conductance Laplacian. Positive current is defined from
    % endpoint i to endpoint j as (vCi-vCj)/R.
    properties
        id      (1,1) double
        i       (1,1) double    % endpoint i
        j       (1,1) double    % endpoint j
        R       (1,1) double = 1.0   % Ohm
        g       (1,1) double = 1.0   % S (filled in ctor)
        color   (1,3) double = [0.3 0.3 0.3];
        width   (1,1) double = 1.6
    end

    methods
        function obj = TransmissionLine(id, i, j, R)
            % REVIEW PROPOSAL LINE-01 (ADD): validate the physical line data.
            % validateattributes(id, {'numeric'}, {'scalar','integer','positive','finite'});
            % validateattributes(i,  {'numeric'}, {'scalar','integer','positive','finite'});
            % validateattributes(j,  {'numeric'}, {'scalar','integer','positive','finite'});
            % validateattributes(R,  {'numeric'}, {'scalar','real','positive','finite'});
            % assert(i ~= j, 'TransmissionLine:SelfLoop', ...
            %     'A physical transmission line must connect two distinct DGs.');

            obj.id = id; 
            obj.i = i; 
            obj.j = j;
            obj.R = R;
            obj.g = 1/R;

            % REVIEW PROPOSAL LINE-02 (DESIGN): R and g are both public and
            % can become inconsistent after construction. Prefer making g a
            % Dependent property with get.g = 1/obj.R, or update both values
            % together whenever R changes.
        end

        

        function i_ij = current(obj, vCi, vCj)
            i_ij = (vCi - vCj)/obj.R;
        end

        function draw(obj, ax, pos_i, pos_j, varargin)
            if nargin < 2 || isempty(ax), ax = gca; end

            parser = inputParser;
            parser.FunctionName = 'TransmissionLine.draw';
            addParameter(parser,'ShowLabel',true, ...
                @(value)islogical(value) && isscalar(value));
            addParameter(parser,'ShowCurrent',false, ...
                @(value)islogical(value) && isscalar(value));
            addParameter(parser,'LineCurrent',NaN, ...
                @(value)isnumeric(value) && isscalar(value));
            addParameter(parser,'FontSize',8, ...
                @(value)isnumeric(value) && isscalar(value) && value>0);
            parse(parser,varargin{:});
            options = parser.Results;

            plot(ax,[pos_i(1) pos_j(1)],[pos_i(2) pos_j(2)],'-', ...
                'Color',obj.color,'LineWidth',obj.width);
            if ~options.ShowLabel
                return;
            end

            mid = (pos_i + pos_j)/2;
            if options.ShowCurrent && isfinite(options.LineCurrent)
                label = sprintf('R=%.2f Ohm\nI%d->%d=%+.2f A', ...
                    obj.R,obj.i,obj.j,options.LineCurrent);
            else
                label = sprintf('R=%.2f Ohm',obj.R);
            end
            text(ax,mid(1),mid(2),label,'HorizontalAlignment','center', ...
                'VerticalAlignment','middle','Color',[0.18 0.21 0.24], ...
                'FontSize',options.FontSize,'BackgroundColor','w', ...
                'Margin',1,'Clipping','off','Interpreter','none');
        end
    end
end

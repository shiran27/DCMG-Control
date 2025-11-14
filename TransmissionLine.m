classdef TransmissionLine < handle
    % Resistive line between DGs (i,j)
    % g = 1/R. Directional Ytilde_ij depends on the RECEIVER'S capacitor C_i:
    %   Ytilde_ij = [ g/C_i  0; 0  0]
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
            obj.id = id; 
            obj.i = i; 
            obj.j = j;
            obj.R = R;
            obj.g = 1/R;
        end

        

        function i_ij = current(obj, vCi, vCj)
            i_ij = (vCi - vCj)/obj.R;
        end

        function draw(obj, ax, pos_i, pos_j)
            if nargin < 2 || isempty(ax), ax = gca; end
            plot(ax, [pos_i(1) pos_j(1)], [pos_i(2) pos_j(2)], '-', ...
                'Color', obj.color, 'LineWidth', obj.width);
            mid = (pos_i + pos_j)/2;
            text(mid(1), mid(2), sprintf('R=%.2fΩ', obj.R), ...
                'HorizontalAlignment','center', 'Color', [0.2 0.2 0.2], 'FontSize',8);
        end
    end
end

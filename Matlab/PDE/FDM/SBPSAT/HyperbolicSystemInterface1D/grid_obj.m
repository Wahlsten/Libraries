classdef grid_obj
    % Grid object
    properties
        minv;
        maxv;
        N;
        d;
        grid;
        gridT;
        dim;
        t;
        Trans; % Transformation parameters
        geometry; % Geometry: 'Square', 'Circle Sector'
    end
    methods
        function plotGrid(obj)
            
            plot(obj.grid, zeros(obj.N(1) + 1, 1), 'b')
            axis equal
            
            axis([obj.minv(1) - .1 * (obj.maxv(1) - obj.minv(1)), ...
                obj.maxv(1) + .1 * (obj.maxv(1) - obj.minv(1)), - .1 * ...
                (obj.maxv(1) - obj.minv(1)), .1 * (obj.maxv(1) - obj.minv(1))])
            axis equal
            
        end
        function obj = createGrid(obj)
            
            % Spatial grid
            
            obj.grid{1} = linspace(obj.minv(1), obj.maxv(1), obj.N(1)+1);
            obj.grid{2} = linspace(obj.minv(2), obj.maxv(2), obj.N(2)+1);
            obj.d(1) = (obj.maxv(1)-obj.minv(1))/(obj.N(1));
            obj.d(2) = (obj.maxv(2)-obj.minv(2))/(obj.N(2));
            obj.dim = 1;
            
            % Time grid
            
            obj.d(end + 1) = (obj.maxv(end) - obj.minv(end)) / (obj.N(end) - 1);
            obj.t          = linspace(obj.minv(end), obj.maxv(end), obj.N(end) + 1);
            
            % Transformation
            
            if ~strcmp(obj.geometry, 'Unit')

                obj = transformGrid(obj);
                
            end
            
        end
        function obj = transformGrid(obj)
            
            s = 1;
            
            obj.grid = obj.grid * (1 + s);
                
        end
    end
end
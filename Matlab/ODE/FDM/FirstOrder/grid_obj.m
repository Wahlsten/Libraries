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
            
            plot(obj.t, zeros(obj.N + 1, 1), 'b')
            axis equal
            
            axis([obj.minv - .1 * (obj.maxv - obj.minv(1)), obj.maxv(1) + .1 * (obj.maxv(1) - obj.minv(1)), - .1 * (obj.maxv(1) - obj.minv(1)), .1 * (obj.maxv(1) - obj.minv(1))])
            axis equal
            
        end
        function obj = createGrid(obj)
            
            obj.d = (obj.maxv - obj.minv) / (obj.N);
            obj.t = linspace(obj.minv, obj.maxv, obj.N + 1);
            
        end
    end
end
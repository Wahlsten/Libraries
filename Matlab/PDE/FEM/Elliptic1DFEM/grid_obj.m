classdef grid_obj
    % Grid object
    properties
        minv;
        maxv;
        Nx;
        dx;
        grid;
    end
    methods
        function plotGrid(obj)
            
            plot(obj.grid, zeros(obj.N(1) + 1, 1), 'b')
            axis equal
            
            axis([obj.minv(1) - .1 * (obj.maxv(1) - obj.minv(1)), obj.maxv(1) + .1 * (obj.maxv(1) - obj.minv(1)), - .1 * (obj.maxv(1) - obj.minv(1)), .1 * (obj.maxv(1) - obj.minv(1))])
            axis equal
            
        end
        function obj = createGrid(obj)
            
            % Spatial grid
            
            obj.grid = linspace(obj.minv(1), obj.maxv(1), obj.Nx(1)+1);
            
            obj.dx = obj.grid(2:end) - obj.grid(1:end-1);
                        
        end
    end
end
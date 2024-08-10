classdef grid_obj
    % Grid object
    properties
        minv;
        maxv;
        N;
        d;
        t;
    end
    methods
        function plotGrid(obj)
            
            plot(obj.grid, zeros(obj.N + 1, 1), 'b')
            axis equal
            
            axis([obj.minv(1) - .1 * (obj.maxv - obj.minv), obj.maxv + .1 * (obj.maxv - obj.minv), - .1 * (obj.maxv - obj.minv), .1 * (obj.maxv - obj.minv)])
            axis equal
            
        end
        function obj = createGrid(obj)
                        
            obj.d = (obj.maxv - obj.minv) / (obj.N);
            obj.t = obj.minv:obj.d:obj.maxv; %linspace(obj.minv, obj.maxv, obj.N + 1);
            
        end
    end
end
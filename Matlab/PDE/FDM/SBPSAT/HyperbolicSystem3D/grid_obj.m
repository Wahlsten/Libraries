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
            
            for k = 1:size(obj.grid, 3)
                
                h = plot3(obj.grid(:,:,k,2), obj.grid(:, :, k, 3), obj.grid(:, :, k, 1), 'g-');
                set(h,{'color'},num2cell(cool(length(h)),2));
                hold on
                h = plot3(obj.grid(:,:,k,3), obj.grid(:, :, k, 2), obj.grid(:, :, k, 1), 'r-');
                set(h,{'color'},num2cell(cool(length(h)),2));
                h = plot3(obj.grid(:,:,k,1), obj.grid(:, :, k, 2), obj.grid(:, :, k, 3), 'b-');
                set(h,{'color'},num2cell(cool(length(h)),2));
                
            end
            hold off
            
        end
        function obj = createGrid(obj)
            
            % Spatial grid
            
            x = linspace(obj.minv(1), obj.maxv(1), obj.N(1)+1);
            y = linspace(obj.minv(2), obj.maxv(2), obj.N(2)+1);
            z = linspace(obj.minv(3), obj.maxv(3), obj.N(3)+1);
            
            obj.d(1) = (obj.maxv(1)-obj.minv(1))/(obj.N(1));
            obj.d(2) = (obj.maxv(2)-obj.minv(2))/(obj.N(2));
            obj.d(3) = (obj.maxv(3)-obj.minv(3))/(obj.N(3));
            
            [obj.grid(:, :, :, 2) , obj.grid(:, :, :, 3),  obj.grid(:, :, :, 1)] = meshgrid(x, y, z);
            obj.dim = sign(obj.N(1)) + sign(obj.N(2)) + sign(obj.N(3));
            
            % Time grid
            
            obj.d(end + 1) = (obj.maxv(end) - obj.minv(end)) / (obj.N(end) - 1);
            obj.t          = linspace(obj.minv(end), obj.maxv(end), obj.N(end) + 1);
        end
    end
end
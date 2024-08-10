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
            
            plot(obj.grid(:,:,1), obj.grid(:, :, 2), 'b')
            hold on
            plot(obj.grid(:,:,1)', obj.grid(:, :, 2)', 'b')
            hold off
            axis equal
            
            dimG = 2;
            minG = zeros(1,dimG);
            maxG = zeros(1,dimG);
            lenG = zeros(1,dimG);
            
            for k = 1:dimG
                
                minG(k) = min(min(obj.grid(:, :, k)));
                maxG(k) = max(max(obj.grid(:, :, k)));
                
                if maxG(k) - minG(k) ~= 0
                    
                    lenG(k) = maxG(k) - minG(k);
                    
                else
                    
                    lenG(k) = 1;
                    
                end
            end
            
            axis([minG(1) - .1 * lenG(1), maxG(1) + .1 * lenG(1), minG(2) - .1 * lenG(2), maxG(2) + .1 * lenG(2)])
            axis equal
            
        end
        function obj = createGrid(obj, s)
            
            % Spatial grid
            
            x = linspace(obj.minv(1), obj.maxv(1), obj.N(1)+1);
            y = linspace(obj.minv(2), obj.maxv(2), obj.N(2)+1);
            
            obj.d(1) = (obj.maxv(1)-obj.minv(1))/(obj.N(1));
            obj.d(2) = (obj.maxv(2)-obj.minv(2))/(obj.N(2));
            
            [obj.grid(:, :, 1), obj.grid(:, :, 2)] = meshgrid(x, y);
            [obj.gridT(:, :, 1), obj.gridT(:, :, 2)] = meshgrid(x, y);
            obj.dim = sign(obj.N(1)) + sign(obj.N(2));
            
            % Time grid
            
            obj.d(end + 1) = (obj.maxv(end) - obj.minv(end)) / (obj.N(end) - 1);
            obj.t          = linspace(obj.minv(end), obj.maxv(end), obj.N(end) + 1);
            
            % Transformation
            if ~strcmp(obj.geometry, 'Square')

                obj = transformGrid(obj, s);
            end
        end
        function obj = transformGrid(obj, s)
            
            if strcmp(obj.geometry, 'Circle Sector')
                
                r0_0  = 1;
                r1_0  = 2;
                fi0_0 = pi/4;
                fi1_0 = 3*pi/4;
                
                obj.Trans.r0  = @(s) r0_0;
                obj.Trans.r1  = @(s) r1_0;
                obj.Trans.fi0 = @(s) fi0_0;
                obj.Trans.fi1 = @(s) fi1_0;
                
                r   = obj.Trans.r0(s)  + obj.grid(:, :, 1) .* (obj.Trans.r1(s)  - obj.Trans.r0(s));
                phi = obj.Trans.fi0(s) + obj.grid(:, :, 2) .* (obj.Trans.fi1(s) - obj.Trans.fi0(s));
                
                obj.grid(:, :, 1) = r.*cos(phi);
                obj.grid(:, :, 2) = r.*sin(phi);
                
            elseif strcmp(obj.geometry, 'Rectangle')
                
                obj.grid(:, :, 1) = obj.grid(:, :, 1);
                obj.grid(:, :, 2) = obj.grid(:, :, 2) * (1 + s);
                
            elseif strcmp(obj.geometry, 'RectangleSine')
                
                obj.grid(:, :, 1) = obj.grid(:, :, 1);
                obj.grid(:, :, 2) = obj.grid(:, :, 2) + (1 - obj.grid(:, :, 2)) .* 0.05 .* s .* sin(2 * pi * obj.grid(:, :, 1));
                
            end
        end
    end
end
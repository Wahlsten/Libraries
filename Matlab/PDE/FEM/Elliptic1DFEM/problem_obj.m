classdef problem_obj
    % problem object
    % ut + Aux + Buy = e(Cuxx + Duyy)
    %       Hu = g
    %        u = f
    
    properties
        % Parameters
        cont;    % Structure with continuous matrices: a, ax, b, by, c, d
        disc;    % Structure with discrete matrices: A, Ax, B, By, C, D
        epsilon; % epsilon
        dim;     % dimension of system
        
        % Problem type
        type;   % Type: problem type, data type, coeff type
        
        % Transformation
        Trans   % Transformation:
        
        % Data
        data;  % Boundary data: gE, gED, gW, gWD, ...
        
    end
    methods
        function obj = initialize(obj)
            
            obj = initializeProblem(obj);
            obj = initializeData(obj);
            
        end
        function obj = initializeProblem(obj)
            
            if strcmp(obj.type.problem, 'Scalar')
                
                obj.cont.a  = @(x) (1 + 0*x);
                obj.cont.b  = @(x) (0 + 0*x);
                            
            end
        end
        function obj = initializeData(obj)
            
            if strcmp(obj.type.data, 'sin')
                
                obj.data.u_a   = @(x) (          sin(2*pi*(x)));
                obj.data.uX_a  = @(x) (     2*pi*cos(2*pi*(x)));
                obj.data.uXX_a = @(x) (-(2*pi)^2*sin(2*pi*(x)));
                
            elseif strcmp(obj.type.data, 'constant')
                
                obj.data.u_a   = @(x) ((x.*(x - 1))/2);
                obj.data.uX_a  = @(x) (x - 1/2 + 0);
                obj.data.uXX_a = @(x) (1 + 0*x);
           
            end
        end
        function obj = createMatrices(obj, grid_obj)
            
            avals  = obj.cont.a (grid_obj.grid(:, :, 1))';
            bvals  = obj.cont.b (grid_obj.grid(:, :, 1))';

            obj.disc.A  = sparse(diag(avals(:)));
            obj.disc.B  = sparse(diag(bvals(:)));
            
        end
        function obj = createData(obj)

            
            obj.data.gE    = @(x) (obj.data.u_a(x));
            obj.data.gW    = @(x) (obj.data.u_a(x));
            obj.data.f     = @(x) (obj.data.u_a(x));
            obj.data.force = @(x) (-obj.cont.a(x) .* obj.data.uXX_a(x) + obj.cont.b(x) .* obj.data.u_a(x));

        end
        function obj = createBasis(obj)
           
            obj.data.hat1  = @(x, x1, x2) ((x - x1)/(x2 - x1));
            obj.data.hat2  = @(x, x1, x2) ((x2 - x)/(x2 - x1));
            
        end
    end
end
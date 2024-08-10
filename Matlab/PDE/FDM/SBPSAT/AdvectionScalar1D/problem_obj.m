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
        
        % Matrices
        Lambda; % Eigenvalue matrices
        X;      % Eigenvectors: C_plus, C_minus, D_plus...
        tilde;  % Matrices: X'CX, X'DX, ...
        
        
    end
    methods
        function obj = initialize(obj)
            
            obj = initializeProblem(obj);
            obj = initializeData(obj);
            
        end
        function obj = initializeProblem(obj)
            
            
            if strcmp(obj.type.problem, 'Advection1D')
                
                obj.dim        = 1;
                obj.cont.a     = @(x) (1 + 0*x);
                 
            end
        end
        function obj = initializeData(obj)
            
            if strcmp(obj.type.data, 'sin')
                
                obj.data.u_a   = @(x, t) (      sin(2*pi*(x - t)));
                obj.data.uT_a  = @(x, t) (-2*pi*cos(2*pi*(x - t)));
                obj.data.uX_a  = @(x, t) ( 2*pi*cos(2*pi*(x - t)));
            
            end
        end
        function obj = createMatrices(obj, grid_obj)
            
            avals  = obj.cont.a (grid_obj.grid(:, :, 1))';
            obj.disc.A  = sparse(diag(avals (:)));
            obj = transformation(obj);
            
        end
        function obj = createData(obj)
                
            obj.data.gE  = @(x, t) (obj.data.u_a (x, t));
            obj.data.gW  = @(x, t) (obj.data.u_a (x, t));
            obj.data.f   = @(x, t) (obj.data.u_a (x, t));

        end
        function obj = transformation(obj)
            
            obj.Trans.A_a    = sparse(obj.disc.A);

        end
    end
end
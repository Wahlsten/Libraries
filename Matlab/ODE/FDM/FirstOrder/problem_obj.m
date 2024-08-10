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
            
            
            if strcmp(obj.type.problem, 'ODE')
                 
            end
        end
        function obj = initializeData(obj)
            
            if strcmp(obj.type.data, 'sin')
                
                obj.data.u_a   = @(t) (      sin(2*pi*t));
                obj.data.uT_a  = @(t) (2*pi*cos(2*pi*t));
                obj.data.force = @(t) (2*pi*cos(2*pi*t));
                
            elseif strcmp(obj.type.data, 'constant')
                
                obj.data.u_a   = @(t) (      1 + 2*t.^2);
                obj.data.uT_a  = @(t) (-2*pi*cos(2*pi*t));
                obj.data.force = @(t) (0 + 4*t);
                
            end
        end
        function obj = createData(obj)
                
            obj.data.f   = @(t) (obj.data.u_a (t));

        end
    end
end
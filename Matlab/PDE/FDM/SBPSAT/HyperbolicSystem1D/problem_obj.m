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
            
           if strcmp(obj.type.problem, 'Euler')
                
                obj.dim     = 3;
                ubar        = 1;
                cbar        = 2;
                gamma       = 1.4;
                %obj.cont.a  = eye(3);
                obj.cont.a  = [1 cbar / sqrt(gamma) 0; ...
                    cbar / sqrt(gamma) ubar cbar * sqrt((gamma - 1)/gamma); ...
                    0 cbar * sqrt((gamma - 1)/gamma) ubar];
                
           elseif strcmp(obj.type.problem, 'Identity')
               
               obj.dim = 2;
               obj.cont.a = [2 0; 0 1];
               
           end
        end
        function obj = initializeData(obj)
            
            if strcmp(obj.type.data, 'sin')
                
                obj.data.u_a   = @(x, t) ([ ...
                    sin(2*pi*(x - t));...
                    sin(2*pi*(x - t));...
                    sin(2*pi*(x - t))]);
                obj.data.uT_a  = @(x, t) ([ ...
                    -2*pi*cos(2*pi*(x - t));...
                    -2*pi*cos(2*pi*(x - t));...
                    -2*pi*cos(2*pi*(x - t))]);
                obj.data.uX_a  = @(x, t, s) ([ ...
                    2*pi*cos(2*pi*(x - t));...
                    2*pi*cos(2*pi*(x - t));...
                    2*pi*cos(2*pi*(x - t))]);
            elseif strcmp(obj.type.data, 'sin2')
                
                obj.data.u_a   = @(x, t) ([ ...
                    sin(2*pi*(x - t));...
                    sin(2*pi*(x - t))]);
                obj.data.uT_a  = @(x, t) ([ ...
                    -2*pi*cos(2*pi*(x - t));...
                    -2*pi*cos(2*pi*(x - t))]);
                obj.data.uX_a  = @(x, t, s) ([ ...
                    2*pi*cos(2*pi*(x - t));...
                    2*pi*cos(2*pi*(x - t))]);
            end
        end
        function obj = createMatrices(obj, grid_obj)
                    
            % Diagonalization A and B
            
            [~, obj.Lambda.A]   = eig(obj.cont.a);
            obj.Lambda.A        = flip(flip(obj.Lambda.A)');

            obj.Lambda.posDimA  = sum(diag(obj.Lambda.A) > 0);
            obj.Lambda.negDimA  = sum(diag(obj.Lambda.A) < 0);
            obj.Lambda.zeroDimA = sum(diag(obj.Lambda.A) == 0);

            obj.Lambda.Aplus  = obj.Lambda.A(1:obj.Lambda.posDimA, 1:obj.Lambda.posDimA);
            obj.Lambda.Aminus = obj.Lambda.A(end-obj.Lambda.negDimA + 1:end, end-obj.Lambda.negDimA + 1: end);
            obj.Lambda.Azero  = zeros(obj.Lambda.zeroDimA);

            obj.Lambda.A = blkdiag(obj.Lambda.Aplus, obj.Lambda.Azero, obj.Lambda.Aminus);

            obj.disc.A  = obj.Lambda.A;
                
            obj = transformation(obj, grid_obj);
            
        end
        function obj = createData(obj)

            obj.data.gE  = @(x, t) (obj.data.u_a (x, t));
            obj.data.gED = @(x, t) (obj.data.uX_a(x, t));
            obj.data.gW  = @(x, t) (obj.data.u_a (x, t));
            obj.data.gWD = @(x, t) (obj.data.uX_a(x, t));
            obj.data.f   = @(x, t) (obj.data.u_a (x, t));

        end
        function obj = transformation(obj, grid_obj)
                        
            Ix            = eye(grid_obj.N(1) + 1);
            obj.Trans.A_a = sparse(kron(Ix, obj.disc.A));
            
        end
    end
end
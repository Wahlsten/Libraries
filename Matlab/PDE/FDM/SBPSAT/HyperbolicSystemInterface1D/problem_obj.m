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
                
                obj.cont.a     = [2 1; 1 2];
                obj.cont.b     = [3 1; 1 3];
                obj.dim        = length(obj.cont.a);
                
            end
        end
        function obj = initializeData(obj)
            
            if strcmp(obj.type.data, 'sin')
                
                obj.data.u_a   = @(x, t) ([...
                    sin(2*pi*(x - t)); ...
                    sin(2*pi*(x - t))]);
                obj.data.uT_a  = @(x, t) ([...
                    -2*pi*cos(2*pi*(x - t)); ...
                    -2*pi*cos(2*pi*(x - t))]);
                obj.data.uX_a  = @(x, t) ([ ...
                    2*pi*cos(2*pi*(x - t)); ...
                    2*pi*cos(2*pi*(x - t))]);
            
            end
        end
        function obj = createMatrices(obj)
            
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
                        
            [~, obj.Lambda.B]   = eig(obj.cont.b);
            obj.Lambda.B        = flip(flip(obj.Lambda.B)');

            obj.Lambda.posDimB  = sum(diag(obj.Lambda.B) > 0);
            obj.Lambda.negDimB  = sum(diag(obj.Lambda.B) < 0);
            obj.Lambda.zeroDimB = sum(diag(obj.Lambda.B) == 0);

            obj.Lambda.Bplus  = obj.Lambda.B(1:obj.Lambda.posDimB, 1:obj.Lambda.posDimB);
            obj.Lambda.Bminus = obj.Lambda.B(end-obj.Lambda.negDimB + 1:end, end-obj.Lambda.negDimB + 1: end);
            obj.Lambda.Bzero  = zeros(obj.Lambda.zeroDimB);

            obj.Lambda.B = blkdiag(obj.Lambda.Bplus, obj.Lambda.Bzero, obj.Lambda.Bminus);

            obj.disc.B  = obj.Lambda.B;
                
            obj = transformation(obj);
            
        end
        function obj = createData(obj)
                
            obj.data.gE  = @(x, t) (obj.data.u_a (x, t));
            obj.data.gW  = @(x, t) (obj.data.u_a (x, t));
            obj.data.f   = @(x, t) (obj.data.u_a (x, t));

        end
        function obj = transformation(obj)
            
            obj.Trans.A_a    = sparse(obj.disc.A);
            obj.Trans.B_a    = sparse(obj.disc.B);

        end
    end
end
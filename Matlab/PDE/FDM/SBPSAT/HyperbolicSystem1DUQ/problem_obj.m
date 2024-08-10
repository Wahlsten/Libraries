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
               
                obj.cont.a     = [2 0; 0 -3];
                obj.dim        = length(obj.cont.a);
                
            end
        end
        function obj = initializeData(obj)
            
            if strcmp(obj.type.data, 'sin')
                
                obj.data.u_a   = @(x, t, s) ([...
                    sin(2*pi*(x - t)) + 1*s; ...
                    sin(2*pi*(x - t)) + 1*s]);
                obj.data.uT_a  = @(x, t, s) ([...
                    -2*pi*cos(2*pi*(x - t)) + 0*s; ...
                    -2*pi*cos(2*pi*(x - t)) + 0*s]);
                obj.data.uX_a  = @(x, t, s) ([ ...
                    2*pi*cos(2*pi*(x - t)) + 0*s; ...
                    2*pi*cos(2*pi*(x - t)) + 0*s]);
                
            elseif strcmp(obj.type.data, 'constant')
                
                obj.data.u_a   = @(x, t, s) ([...
                    3 + 1*x + 0*s + 0*t; ...
                    3 + 1*x + 0*s + 0*t]);
                obj.data.uT_a  = @(x, t, s) ([...
                    0 + 0*x + 0*t + 0*s; ...
                    0 + 0*x + 0*t + 0*s]);
                obj.data.uX_a  = @(x, t, s) ([ ...
                    1 + 0*x + 0*t + 0*s; ...
                    1 + 0*x + 0*t + 0*s]);
                
            elseif strcmp(obj.type.data, 'sincos_smooth')
                
                obj.data.u_a   = @(x, t, s) ([ ...
                    sin(2*pi*(x - t)) + s .* sin(s/2 - t) .* cos(x - s); ...
                    sin(2*pi*(x - t)) + s .* sin(s/2 - t) .* cos(x - s)]);
                obj.data.uT_a  = @(x, t, s) ([ ...
                    -2*pi*cos(2*pi*(x - t)) - s .* cos(s/2 - t) .* cos(x - s); ...
                    -2*pi*cos(2*pi*(x - t)) - s .* cos(s/2 - t) .* cos(x - s)]);
                obj.data.uX_a  = @(x, t, s) ([ ...
                    2*pi*cos(2*pi*(x - t)) - s .* sin(s/2 - t) .* sin(x - s);
                    2*pi*cos(2*pi*(x - t)) - s .* sin(s/2 - t) .* sin(x - s)]);
                
            elseif strcmp(obj.type.data, 'smooth')
                
                eps = 1;
                
                obj.data.u_a   = @(x, t, s) ([ ...
                    sin(2*pi*(x - t)) + 1./(1 + exp(-s/(10*eps))); ...
                    sin(2*pi*(x - t)) + 1./(1 + exp(-s/(10*eps)))]);
                obj.data.uT_a  = @(x, t, s) ([ ...
                    -2*pi*cos(2*pi*(x - t)); ...
                    -2*pi*cos(2*pi*(x - t))]);
                obj.data.uX_a  = @(x, t, s) ([ ...
                    2*pi*cos(2*pi*(x - t)); ...
                    2*pi*cos(2*pi*(x - t))]);
            
            elseif strcmp(obj.type.data, 'non_smooth')
                
                eps = .01;
                
                obj.data.u_a   = @(x, t, s) ([ ...
                    sin(2*pi*(x - t)) + exp(-(x - 1/2).^2*100)./(1 + exp(-s/(10*eps))); ...
                    sin(2*pi*(x - t)) + exp(-(x - 1/2).^2*100)./(1 + exp(-s/(10*eps)))]);
                obj.data.uT_a  = @(x, t, s) ([ ...
                    -2*pi*cos(2*pi*(x - t)); ...
                    -2*pi*cos(2*pi*(x - t))]);
                obj.data.uX_a  = @(x, t, s) ([ ...
                    2*pi*cos(2*pi*(x - t)) - 200 * (x - 1/2) .* exp(-(x - 1/2).^2*100)./(1 + exp(-s/(10*eps))); ...
                    2*pi*cos(2*pi*(x - t)) - 200 * (x - 1/2) .* exp(-(x - 1/2).^2*100)./(1 + exp(-s/(10*eps)))]);
                
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

            obj = transformation(obj);
            
        end
        function obj = createData(obj)
                
            obj.data.gE  = @(x, t, s) (obj.data.u_a (x, t, s));
            obj.data.gW  = @(x, t, s) (obj.data.u_a (x, t, s));
            obj.data.f   = @(x, t, s) (obj.data.u_a (x, t, s));

        end
        function obj = transformation(obj)
            
            obj.Trans.A_a    = sparse(obj.disc.A);

        end
    end
end
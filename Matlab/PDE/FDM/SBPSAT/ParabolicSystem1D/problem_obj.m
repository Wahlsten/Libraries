classdef problem_obj
    % problem object
    % ut + Aux = eBuxx
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
            
            if strcmp(obj.type.problem, 'Navier-Stokes')
                
                obj.epsilon     = 0.01;
                obj.dim         = 3;
                obj.type.coeff  = 'constant';
                
                ubar   = 1;
                % pbar = 1;
                rhobar = 1;
                cbar   = 2;
                gamma  = 1.4;
                lam    = -2/3;
                mu     = 1;
                Pr     = 0.7;
                
                obj.cont.a = [1 cbar / sqrt(gamma) 0; ...
                    cbar / sqrt(gamma) ubar cbar * sqrt((gamma - 1)/gamma); ...
                    0 cbar * sqrt((gamma - 1)/gamma) ubar];
                obj.cont.b = [0 0 0; ...
                    0 (lam + 2 * mu) / rhobar 0; ...
                    0 0 gamma * mu / (Pr * rhobar)];
                
            end
        end
        function obj = initializeData(obj)
            
            if strcmp(obj.type.data, 'Sin')
                
                obj.data.u_a   = @(x, t) ([ ...
                    sin(2*pi*(x - t));...
                    sin(2*pi*(x - t));...
                    sin(2*pi*(x - t))]);
                obj.data.uT_a  = @(x, t) ([ ...
                    -2*pi*cos(2*pi*(x - t));...
                    -2*pi*cos(2*pi*(x - t));...
                    -2*pi*cos(2*pi*(x - t))]);
                obj.data.uX_a  = @(x, t) ([ ...
                    2*pi*cos(2*pi*(x - t));...
                    2*pi*cos(2*pi*(x - t));...
                    2*pi*cos(2*pi*(x - t))]);
                obj.data.uXX_a = @(x, t) ([ ...
                    -(2*pi)^2*sin(2*pi*(x - t));...
                    -(2*pi)^2*sin(2*pi*(x - t));...
                    -(2*pi)^2*sin(2*pi*(x - t))]);
            
            end
        end
        function obj = createMatrices(obj, grid_obj)

            if strcmp(obj.type.coeff, 'variable')

                if grid_obj.dim == 1

                    avals  = obj.cont.a (grid_obj.grid(:, :, 1))';
                    bvals  = obj.cont.b (grid_obj.grid(:, :, 1))';

                else

                    avals  = obj.cont.a (grid_obj.grid(:, :, 1), grid_obj.grid(:, :, 2))';
                    bvals  = obj.cont.b (grid_obj.grid(:, :, 1), grid_obj.grid(:, :, 2))';

                end

                obj.disc.A  = sparse(diag(avals (:)));
                obj.disc.B  = sparse(diag(bvals (:)));

            elseif strcmp(obj.type.coeff, 'constant')

                % Diagonalization A and B
                
                [X_A, obj.Lambda.A] = eig(obj.cont.a);
                obj.Lambda.A        = flip(flip(obj.Lambda.A)');
                obj.Lambda.posDimA  = sum(diag(obj.Lambda.A) > 0);
                obj.Lambda.negDimA  = sum(diag(obj.Lambda.A) < 0);
                obj.Lambda.zeroDimA = sum(diag(obj.Lambda.A) == 0);

                obj.Lambda.Aplus  = obj.Lambda.A(1:obj.Lambda.posDimA, 1:obj.Lambda.posDimA);
                obj.Lambda.Aminus = obj.Lambda.A(end-obj.Lambda.negDimA + 1:end, end-obj.Lambda.negDimA + 1: end);
                obj.Lambda.Azero  = zeros(obj.Lambda.zeroDimA);

                X_Aplus  = X_A(:, 1:obj.Lambda.posDimA);
                X_Aminus = X_A(:, end-obj.Lambda.negDimA + 1:end);
                X_Azero  = zeros(obj.dim, obj.Lambda.zeroDimA);
                obj.Lambda.A = blkdiag(obj.Lambda.Aplus, obj.Lambda.Azero, obj.Lambda.Aminus);

                X_A = [X_Aplus, X_Azero, X_Aminus];
                Btilde = X_A' * obj.cont.b * X_A;
                LB = -Btilde * pinv(obj.Lambda.A) * Btilde;
                [X_B, Lambda_B] = eig(LB);

                X_B = X_B(:, flip(1:end));
                Lambda_B = flip(flip(Lambda_B)');
                posDimB  = sum( diag(Lambda_B) > 1e-10);
                negDimB  = sum(-diag(Lambda_B) > 1e-10);

                obj.Lambda.zeroDimB = sum(abs(diag(Lambda_B)) < 1e-10);
                obj.Lambda.Bplus  = Lambda_B(1:posDimB, 1:posDimB);
                obj.Lambda.Bminus = Lambda_B(end-negDimB + 1:end, end-negDimB + 1: end);
                obj.X.Bplus  = X_B(:, 1:posDimB);
                obj.X.Bminus = X_B(:, end-negDimB + 1:end);
                
                obj.disc.A  = obj.Lambda.A;
                obj.disc.B  = Btilde;

            end
            
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
            
            if strcmp(obj.type.coeff, 'variable')

                obj.Trans.A_a    = sparse(obj.disc.A);
                obj.Trans.B_a    = sparse(obj.disc.B);

            elseif strcmp(obj.type.coeff, 'constant')

                Ix            = eye(grid_obj.N(1) + 1);
                obj.Trans.A_a = sparse(kron(Ix, obj.disc.A));
                obj.Trans.B_a = sparse(kron(Ix, obj.disc.B));

            end
        end
    end
end
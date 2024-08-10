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
            
            if strcmp(obj.type.problem, 'Maxwell')
                
                obj.dim         = 3;
                obj.type.coeff  = 'constant';
                obj.type.linear = 'linear';
                ep              = 1;
                mu              = 1;
                obj.cont.s      = [mu 0 0; 0 ep 0; 0 0 ep]; % not used
                obj.cont.a      = [0 0 1; 0 0 0; 1 0 0];
                obj.cont.b      = [0 -1 0; -1 0 0; 0 0 0];
                obj.cont.c      = zeros(obj.dim);
                obj.cont.d      = zeros(obj.dim);
                obj.epsilon     = 1;
                 
            end
        end
        function obj = initializeData(obj)

            if strcmp(obj.type.data, 'sin2DMaxwell')
                
                obj.data.u_a   = @(x, y, t, s) ([ ...
                    sin(2*pi*(x - t)) + sin(2*pi*(y - t));...
                    sin(2*pi*(x - t)) + sin(2*pi*(y - t));...
                    sin(2*pi*(x - t)) + sin(2*pi*(y - t))]);
                obj.data.uT_a  = @(x, y, t, s) ([ ...
                    -2*pi*cos(2*pi*(x - t)) + -2*pi*cos(2*pi*(y - t));...
                    -2*pi*cos(2*pi*(x - t)) + -2*pi*cos(2*pi*(y - t));...
                    -2*pi*cos(2*pi*(x - t)) + -2*pi*cos(2*pi*(y - t))]);
                obj.data.uX_a  = @(x, y, t, s) ([ ...
                    2*pi*cos(2*pi*(x - t)) + y*0;...
                    2*pi*cos(2*pi*(x - t)) + y*0;...
                    2*pi*cos(2*pi*(x - t)) + y*0]);
                obj.data.uY_a = @(x, y, t) ([ ...
                    2*pi*cos(2*pi*(y - t)) + x*0;...
                    2*pi*cos(2*pi*(y - t)) + x*0;...
                    2*pi*cos(2*pi*(y - t)) + x*0]);
                
            end
        end
        function obj = createMatrices(obj, grid_obj)

            % Diagonalization A and B
            [~, obj.Lambda.A]   = eig(obj.cont.a);
            [~, obj.Lambda.B]   = eig(obj.cont.b);
            obj.Lambda.A        = flip(flip(obj.Lambda.A)');
            obj.Lambda.B        = flip(flip(obj.Lambda.B)');

            obj.Lambda.posDimA  = sum(diag(obj.Lambda.A) > 0);
            obj.Lambda.posDimB  = sum(diag(obj.Lambda.B) > 0);
            obj.Lambda.negDimA  = sum(diag(obj.Lambda.A) < 0);
            obj.Lambda.negDimB  = sum(diag(obj.Lambda.B) < 0);
            obj.Lambda.zeroDimA = sum(diag(obj.Lambda.A) == 0);
            obj.Lambda.zeroDimB = sum(diag(obj.Lambda.B) == 0);

            obj.Lambda.Aplus  = obj.Lambda.A(1:obj.Lambda.posDimA, 1:obj.Lambda.posDimA);
            obj.Lambda.Bplus  = obj.Lambda.B(1:obj.Lambda.posDimB, 1:obj.Lambda.posDimB);
            obj.Lambda.Aminus = obj.Lambda.A(end-obj.Lambda.negDimA + 1:end, end-obj.Lambda.negDimA + 1: end);
            obj.Lambda.Bminus = obj.Lambda.B(end-obj.Lambda.negDimB + 1:end, end-obj.Lambda.negDimB + 1: end);
            obj.Lambda.Azero  = zeros(obj.Lambda.zeroDimA);
            obj.Lambda.Bzero  = zeros(obj.Lambda.zeroDimB);

            obj.Lambda.A = blkdiag(obj.Lambda.Aplus, obj.Lambda.Azero, obj.Lambda.Aminus);
            obj.Lambda.B = blkdiag(obj.Lambda.Bplus, obj.Lambda.Bzero, obj.Lambda.Bminus);

            obj.disc.A  = obj.Lambda.A;
            obj.disc.B  = obj.Lambda.B;
            obj.disc.Ax = zeros(size(obj.disc.A));
            obj.disc.By = zeros(size(obj.disc.B));

            obj = transformation(obj, grid_obj);
            
        end
        function obj = createData(obj)
            
            obj.data.gE   = @(x, y, t) (obj.data.u_a (x, y, t));
            obj.data.gEDx = @(x, y, t) (obj.data.uX_a(x, y, t));
            obj.data.gEDy = @(x, y, t) (obj.data.uY_a(x, y, t));
            obj.data.gW   = @(x, y, t) (obj.data.u_a (x, y, t));
            obj.data.gWDx = @(x, y, t) (obj.data.uX_a(x, y, t));
            obj.data.gWDy = @(x, y, t) (obj.data.uY_a(x, y, t));
            obj.data.gN   = @(x, y, t) (obj.data.u_a (x, y, t));
            obj.data.gNDx = @(x, y, t) (obj.data.uX_a(x, y, t));
            obj.data.gNDy = @(x, y, t) (obj.data.uY_a(x, y, t));
            obj.data.gS   = @(x, y, t) (obj.data.u_a (x, y, t));
            obj.data.gSDx = @(x, y, t) (obj.data.uX_a(x, y, t));
            obj.data.gSDy = @(x, y, t) (obj.data.uY_a(x, y, t));
            obj.data.f    = @(x, y, t) (obj.data.u_a (x, y, t));

        end
        function obj = transformation(obj, grid_obj, s)
            
            if strcmp(grid_obj.geometry, 'Square')
                
                Ix            = eye(grid_obj.N(1) + 1);
                Iy            = eye(grid_obj.N(2) + 1);
                obj.Trans.A_a = sparse(kron(kron(Ix, Iy), obj.disc.A));
                obj.Trans.B_a = sparse(kron(kron(Ix, Iy), obj.disc.B));
                
            elseif strcmp(grid_obj.geometry, 'RectangleSine')
                
                A1p  = diag(obj.disc.A);
                A2p  = diag(obj.disc.B);
                NxNy = (grid_obj.N(1) + 1) * (grid_obj.N(2) + 1);
                ix   = ones(1, grid_obj.N(1) + 1)';
                iy   = ones(1, grid_obj.N(2) + 1)';
                
                X = grid_obj.grid(:,:,1);
                Y = grid_obj.grid(:,:,2);
                
                X0 = kron(kron(0, ix), iy);
                X1 = kron(kron(1, ix), iy);
                Y0 = kron(kron(0, ix), iy);
                Y1 = kron(kron(1, ix), iy);
                
                X_XI  = X1 - X0;
                Y_XI  = 0*Y(:);
                
                X_ETA = 0*X(:);
                Y_ETA = Y1 - Y0;
                
                Xi_X = 1 ./ (X1-X0);
                Xi_Y = 0 * Y(:);
                Yi_X = 0 * X(:);
                Yi_Y = 1 ./ (Y1-Y0);
                
                % Analytical
                J_a       =   1 ./ (Xi_X .* Yi_Y - Yi_X .* Xi_Y);
                J_xi_x_a  =   Y_ETA;
                J_xi_y_a  = - X_ETA;
                J_yi_x_a  = - Y_XI;
                J_yi_y_a  =   X_XI;
                
                %% New domain
                
                J_aM = spdiags(J_a, 0, NxNy, NxNy);
                A_a1 = spdiags(J_xi_x_a .* A1p + J_xi_y_a .* A2p, 0, NxNy, NxNy);
                A_a2 = spdiags(J_yi_x_a .* A1p + J_yi_y_a .* A2p, 0, NxNy, NxNy);

                obj.Trans.J_aM = J_aM;
                obj.Trans.A_a  = A_a1;
                obj.Trans.B_a  = A_a2;
                
            else
                
                A1p  = diag(obj.disc.A);
                A2p  = diag(obj.disc.B);
                NxNy = (grid_obj.N(1) + 1) * (grid_obj.N(2) + 1);
                ix   = ones(grid_obj.N(1) + 1, 1);
                iy   = ones(grid_obj.N(2) + 1, 1);
                
                Rn  = kron((linspace(grid_obj.Trans.r0(s), grid_obj.Trans.r1(s), grid_obj.N(1) + 1)'), iy );
                FIn = kron(ix, (linspace(grid_obj.Trans.fi0(s), grid_obj.Trans.fi1(s), grid_obj.N(2) + 1)'));
                
                R   = Rn  - 0.05 * s * sin(2 * pi * (FIn - grid_obj.Trans.fi0(s))./(grid_obj.Trans.fi1(s) - grid_obj.Trans.fi0(s)));
                FI  = FIn + 0.05 * s * sin(2 * pi * (Rn  - grid_obj.Trans.r0 (s))./(grid_obj.Trans.r1 (s) - grid_obj.Trans.r0 (s)));
                
                R0  = kron(kron(grid_obj.Trans.r0 (s), ix), iy);
                R1  = kron(kron(grid_obj.Trans.r1 (s), ix), iy);
                FI0 = kron(kron(grid_obj.Trans.fi0(s), ix), iy);
                FI1 = kron(kron(grid_obj.Trans.fi1(s), ix), iy);
                
                R_XI =  R1 - R0;
                X_XI =  R_XI .* cos(FI);
                Y_XI =  R_XI .* sin(FI);
                
                FI_ETA = FI1 - FI0;
                X_ETA  = -R .* FI_ETA .* sin(FI);
                Y_ETA  =  R .* FI_ETA .* cos(FI);
                
                Xi_X = cos(FI)  ./ (R1-R0);
                Xi_Y = sin(FI)  ./ (R1-R0);
                Yi_X = ((-1./R) .* sin(FI)) ./ (FI1-FI0);
                Yi_Y = (( 1./R) .* cos(FI)) ./ (FI1-FI0);
                
                % Analytical
                J_a       =   1 ./ (Xi_X .* Yi_Y - Yi_X .* Xi_Y);
                J_xi_x_a  =   Y_ETA;
                J_xi_y_a  = - X_ETA;
                J_yi_x_a  = - Y_XI;
                J_yi_y_a  =   X_XI;
                
                %% New domain
                
                J_aM = spdiags(J_a, 0, NxNy, NxNy);
                A_a1 = spdiags(J_xi_x_a .* A1p + J_xi_y_a .* A2p, 0, NxNy, NxNy);
                A_a2 = spdiags(J_yi_x_a .* A1p + J_yi_y_a .* A2p, 0, NxNy, NxNy);

                obj.Trans.J_aM = J_aM;
                obj.Trans.A_a  = A_a1;
                obj.Trans.B_a  = A_a2;
                
            end
        end
    end
end
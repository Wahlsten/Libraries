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
            
            if strcmp(obj.type.problem, 'Navier-Stokes2D')
                
                obj.epsilon    = 0.01;
                obj.dim        = 4;
                obj.type.coeff = 'constant';
                obj.type.linear = 'linear';
                ubar   = 1;
                vbar   = 1;
                % pbar = 1;
                rhobar = 1;
                cbar   = 2;
                gamma  = 1.4;
                lam    = -2/3;
                mu     = 1;
                Pr     = 0.7;
                
                obj.cont.a = [ubar cbar / sqrt(gamma) 0 0; ...
                    cbar / sqrt(gamma) ubar 0 cbar * sqrt((gamma - 1)/gamma); ...
                    0 0 ubar 0; ...
                    0 cbar * sqrt((gamma - 1)/gamma) 0 ubar];
                obj.cont.b = [vbar 0 cbar / sqrt(gamma) 0; ...
                    0 vbar 0 0; ...
                    cbar / sqrt(gamma) 0 vbar cbar * sqrt((gamma - 1)/gamma); ...
                    0 0 cbar * sqrt((gamma - 1)/gamma) vbar];
                obj.cont.c = [0 0 0 0; ...
                    0 (lam + 2 * mu) / rhobar 0 0; ...
                    0 0 mu 0; ...
                    0 0 0 gamma * mu / (Pr * rhobar)];
                obj.cont.d = [0 0 0 0; ...
                    0 mu 0 0; ...
                    0 0 (lam + 2 * mu) / rhobar 0; ...
                    0 0 0 gamma * mu / (Pr * rhobar)];
                obj.cont.f = [0 0 0 0; ...
                    0 0 lam/rhobar 0; ...
                    0 mu 0 0; ...
                    0 0 0 0];
                obj.cont.g = [0 0 0 0; ...
                    0 0 lam/rhobar 0; ...
                    0 mu/rhobar 0 0; ...
                    0 0 0 0];
            end
        end
        function obj = initializeData(obj)
            
            if strcmp(obj.type.data, 'NSXY')
                
                obj.data.u_a   = @(x, y, t) ([   ...
                    sin(x) .* sin(y - t), ...
                    sin(x) .* sin(y - t), ...
                    sin(x) .* sin(y - t), ...
                    sin(x) .* sin(y - t)]);
                obj.data.uT_a  = @(x, y, t) ([ ...
                    - sin(x) .* cos(y - t), ...
                    - sin(x) .* cos(y - t), ...
                    - sin(x) .* cos(y - t), ...
                    - sin(x) .* cos(y - t)]);
                obj.data.uY_a  = @(x, y, t) ([   ...
                    sin(x) .* cos(y - t), ...
                    sin(x) .* cos(y - t), ...
                    sin(x) .* cos(y - t), ...
                    sin(x) .* cos(y - t)]);
                obj.data.uYY_a = @(x, y, t) ([   ...
                    - sin(x) .* sin(y - t), ...
                    - sin(x) .* sin(y - t), ...
                    - sin(x) .* sin(y - t), ...
                    - sin(x) .* sin(y - t)]);
                obj.data.uX_a = @(x, y, t) ([ ...
                    cos(x) .* sin(y - t), ...
                    cos(x) .* sin(y - t), ...
                    cos(x) .* sin(y - t), ...
                    cos(x) .* sin(y - t)]);
                obj.data.uXX_a = @(x, y, t) ([ ...
                     - sin(x) .* sin(y - t), ...
                     - sin(x) .* sin(y - t), ...
                     - sin(x) .* sin(y - t), ...
                     - sin(x) .* sin(y - t)]);
                obj.data.uXY_a = @(x, y, t) ([ ...
                    cos(x).*cos(y - t), ...
                    cos(x).*cos(y - t), ...
                    cos(x).*cos(y - t), ...
                    cos(x).*cos(y - t)]);
                obj.data.uYX_a = @(x, y, t) ([ ...
                    cos(x).*cos(y - t), ...
                    cos(x).*cos(y - t), ...
                    cos(x).*cos(y - t), ...
                    cos(x).*cos(y - t)]);

            end
        end
        function obj = createMatrices(obj, grid_obj)
            
            % Diagonalization A and B
            [X_A, obj.Lambda.A] = eig(obj.cont.a);
            [X_B, obj.Lambda.B] = eig(obj.cont.b);
            %X_A                 = X_A(:, flip(1:end));
            %X_B                 = X_B(:, flip(1:end));
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

            X_Aplus  = X_A(:, 1:obj.Lambda.posDimA);
            X_Bplus  = X_B(:, 1:obj.Lambda.posDimB);
            X_Aminus = X_A(:, end-obj.Lambda.negDimA + 1:end);
            X_Bminus = X_B(:, end-obj.Lambda.negDimB + 1:end);
            X_Azero  = zeros(obj.dim, obj.Lambda.zeroDimA);
            X_Bzero  = zeros(obj.dim, obj.Lambda.zeroDimB);

            obj.Lambda.A = blkdiag(obj.Lambda.Aplus, obj.Lambda.Azero, obj.Lambda.Aminus);
            obj.Lambda.B = blkdiag(obj.Lambda.Bplus, obj.Lambda.Bzero, obj.Lambda.Bminus);

            X_A = [X_Aplus, X_Azero, X_Aminus];
            X_B = [X_Bplus, X_Bzero, X_Bminus];

            Ctilde = X_A' * obj.cont.c * X_A;
            Dtilde = X_B' * obj.cont.d * X_B;
            Ftilde = X_A' * obj.cont.f * X_A;
            Gtilde = X_B' * obj.cont.g * X_B;

            LC = -Ctilde * pinv(obj.Lambda.A) * Ctilde;
            LD = -Dtilde * pinv(obj.Lambda.B) * Dtilde;
            LF = -Ftilde * pinv(obj.Lambda.A) * Ftilde;
            LG = -Gtilde * pinv(obj.Lambda.B) * Gtilde;

            [X_C, Lambda_C] = eig(LC);
            [X_D, Lambda_D] = eig(LD);
            [X_F, Lambda_F] = eig(LF);
            [X_G, Lambda_G] = eig(LG);

            X_C = X_C(:, flip(1:end));
            X_D = X_D(:, flip(1:end));
            X_F = X_F(:, flip(1:end));
            X_G = X_G(:, flip(1:end));

            Lambda_C = flip(flip(Lambda_C)');
            Lambda_D = flip(flip(Lambda_D)');
            Lambda_F = flip(flip(Lambda_F)');
            Lambda_G = flip(flip(Lambda_G)');

            posDimC  = sum( diag(Lambda_C) > 1e-10);
            posDimD  = sum( diag(Lambda_D) > 1e-10);
            posDimF  = sum( diag(Lambda_F) > 1e-10);
            posDimG  = sum( diag(Lambda_G) > 1e-10);
            negDimC  = sum(-diag(Lambda_C) > 1e-10);
            negDimD  = sum(-diag(Lambda_D) > 1e-10);
            negDimF  = sum(-diag(Lambda_F) > 1e-10);
            negDimG  = sum(-diag(Lambda_G) > 1e-10);

            obj.Lambda.zeroDimC = sum(abs(diag(Lambda_C)) < 1e-10);
            obj.Lambda.zeroDimD = sum(abs(diag(Lambda_D)) < 1e-10);
            obj.Lambda.zeroDimF = sum(abs(diag(Lambda_F)) < 1e-10);
            obj.Lambda.zeroDimG = sum(abs(diag(Lambda_G)) < 1e-10);

            obj.Lambda.Cplus  = Lambda_C(1:posDimC, 1:posDimC);
            obj.Lambda.Dplus  = Lambda_D(1:posDimD, 1:posDimD);
            obj.Lambda.Fplus  = Lambda_F(1:posDimF, 1:posDimF);
            obj.Lambda.Gplus  = Lambda_G(1:posDimG, 1:posDimG);

            obj.Lambda.Cminus = Lambda_C(end-negDimC + 1:end, end-negDimC + 1: end);
            obj.Lambda.Dminus = Lambda_D(end-negDimD + 1:end, end-negDimD + 1: end);
            obj.Lambda.Fminus = Lambda_F(end-negDimF + 1:end, end-negDimF + 1: end);
            obj.Lambda.Gminus = Lambda_G(end-negDimG + 1:end, end-negDimG + 1: end);

            %Lambda_Dzero  = Lambda_D(posDimD + 1:posDimD + obj.Lambda.zeroDimD, posDimD + 1:posDimD + obj.Lambda.zeroDimD);

            obj.X.Cplus  = X_C(:, 1:posDimC);
            obj.X.Dplus  = X_D(:, 1:posDimD);
            obj.X.Fplus  = X_F(:, 1:posDimF);
            obj.X.Gplus  = X_G(:, 1:posDimG);

            obj.X.Cminus = X_C(:, end-negDimC + 1:end);
            obj.X.Dminus = X_D(:, end-negDimD + 1:end);
            obj.X.Fminus = X_F(:, end-negDimF + 1:end);
            obj.X.Gminus = X_G(:, end-negDimG + 1:end);

            %X_Dzero  = X_D(:, posDimD + 1: posDimD + obj.Lambda.zeroDimD);
            %X_D = [obj.X_Dplus, obj.X_Dminus];

            obj.disc.A  = obj.Lambda.A;
            obj.disc.B  = obj.Lambda.B;
            obj.disc.C  = Ctilde;
            obj.disc.D  = Dtilde;
            obj.disc.F  = Ftilde;
            obj.disc.G  = Gtilde;
            obj.disc.Ax = zeros(size(obj.disc.A));
            obj.disc.By = zeros(size(obj.disc.B));
            
            s = 0;
            obj = transformation(obj, grid_obj, s);
            
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
                obj.Trans.C_a = sparse(kron(kron(Ix, Iy), obj.disc.C));
                obj.Trans.D_a = sparse(kron(kron(Ix, Iy), obj.disc.D));
                obj.Trans.F_a = sparse(kron(kron(Ix, Iy), obj.disc.F));
                obj.Trans.G_a = sparse(kron(kron(Ix, Iy), obj.disc.G));
                
            elseif strcmp(grid_obj.geometry, 'RectangleSine')
                
                A1p  = diag(obj.disc.A);
                A2p  = diag(obj.disc.B);
                B1   = diag(obj.disc.C);
                B2   = diag(obj.disc.D);
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
                B_a1 = spdiags(J_xi_x_a .* J_xi_x_a ./ J_a .* B1 ...
                    +          J_xi_y_a .* J_xi_y_a ./ J_a .* B2, 0, NxNy, NxNy);
                B_a2 = spdiags(J_yi_x_a .* J_yi_x_a ./ J_a .* B1 ...
                    +          J_yi_y_a .* J_yi_y_a ./ J_a .* B2, 0, NxNy, NxNy);
                
                obj.Trans.J_aM = J_aM;
                obj.Trans.A_a  = A_a1;
                obj.Trans.B_a  = A_a2;
                obj.Trans.C_a  = B_a1;
                obj.Trans.D_a  = B_a2;
                obj.Trans.F_a  = zeros(size(B_a1));
                obj.Trans.G_a  = zeros(size(B_a2));
                
            else
                
                A1p  = diag(obj.disc.A);
                A2p  = diag(obj.disc.B);
                B1   = diag(obj.disc.C);
                B2   = diag(obj.disc.D);
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
                B_a1 = spdiags(J_xi_x_a .* J_xi_x_a ./ J_a .* B1 ...
                    +          J_xi_y_a .* J_xi_y_a ./ J_a .* B2, 0, NxNy, NxNy);
                B_a2 = spdiags(J_yi_x_a .* J_yi_x_a ./ J_a .* B1 ...
                    +          J_yi_y_a .* J_yi_y_a ./ J_a .* B2, 0, NxNy, NxNy);
                
                obj.Trans.J_aM = J_aM;
                obj.Trans.A_a  = A_a1;
                obj.Trans.B_a  = A_a2;
                obj.Trans.C_a  = B_a1;
                obj.Trans.D_a  = B_a2;
                obj.Trans.F_a  = zeros(size(B_a1));
                obj.Trans.G_a  = zeros(size(B_a2));
                
            end
        end
    end
end
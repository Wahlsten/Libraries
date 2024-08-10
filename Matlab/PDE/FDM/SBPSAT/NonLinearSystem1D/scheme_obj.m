classdef scheme_obj
    % Numerical formulation
    % ut + Au = 0
    properties
        % Numerical scheme
        A;                      % A in  cA = Au + pen*g
        cA;                     % cA in cA = Au + pen*g
        Sigma;                  % Penalty matrices
        Matrices;               % Matrices: I, Isys
        Lambda;                 % Eigenvalues: W_p, E_m
        Boundary_Operators_n;   % Discrete boundary operators: L_W, L_E
        Boundary_Operators_a;   % Continuous boundary operators:
        Scheme;                 % Scheme matrices:
        Trans;                  % Transformation: A, B, C, D, JdetM
        
        % Numerical solution
        Solution;               % Structure with solutions: u_a, u_n, u_ns, e_L2
        
    end
    methods
        function obj = createAMatrix(obj, grid_obj, sbp_obj, problem_obj, ~)
            
            % Computes matrix A in : Ut + AU = Pen
            
            if strcmp(problem_obj.type.coeff, 'constant') && nargin == 4
                
                E_0            = zeros(grid_obj.N(1) + 1);
                E_N            = zeros(grid_obj.N(1) + 1);
                E_0(1, 1)      = 1;
                E_N(end, end)  = 1;
                obj.Matrices.I = speye(grid_obj.N(1) + 1);

                obj.Matrices.Isys = speye(problem_obj.dim);

                obj.Scheme.E_East = sparse(kron(E_N, obj.Matrices.Isys));
                obj.Scheme.E_West = sparse(kron(E_0, obj.Matrices.Isys));
               
                obj.Matrices.Pinv_sys = sparse(kron(sbp_obj.Pinv, obj.Matrices.Isys));
                obj.Matrices.Dx       = sparse(kron(sbp_obj.D1, obj.Matrices.Isys));
                
                if strcmp(problem_obj.type.dissipation, 'dissipation')
                    
                    obj.Matrices.Pinv_Diss = sparse(kron(grid_obj.d(1) * sbp_obj.Pinv * sbp_obj.DI, obj.Matrices.Isys));
                    
                else
                    
                    obj.Matrices.Pinv_Diss = 0;
                    
                end
                
            elseif strcmp(problem_obj.type.coeff, 'constant') && nargin == 5
                
                obj.Scheme.PinvENSigmaE = sparse(obj.Matrices.Pinv_sys * obj.Sigma.E * obj.Scheme.E_East);
                obj.Scheme.PinvE0SigmaW = sparse(obj.Matrices.Pinv_sys * obj.Sigma.W * obj.Scheme.E_West);

                obj.Scheme.PinvE0SigmaWH0 = sparse(obj.Scheme.PinvE0SigmaW * obj.Boundary_Operators_n.W);
                obj.Scheme.PinvENSigmaEH0 = sparse(obj.Scheme.PinvENSigmaE * obj.Boundary_Operators_n.E);

                obj.A = - (1/3 * obj.Scheme.Am * obj.Matrices.Dx + 1/3 * obj.Matrices.Dx * obj.Scheme.Am) ...
                    + obj.Scheme.PinvE0SigmaW * obj.Boundary_Operators_n.W ...
                    + obj.Scheme.PinvENSigmaE * obj.Boundary_Operators_n.E;
 
%                 diss = -min(real(eig(full(-obj.A)))) * (min(real(eig(full(-obj.A)))) < 0);
%                 
%                 obj.A = obj.A - diss*obj.Matrices.Pinv_Diss;
                
                obj.A = sparse(obj.A);

            end
        end
        function obj = createScheme(obj, grid_obj, sbp_obj, problem_obj, uq_obj)
            
            X = grid_obj.grid;        
            xlen     = grid_obj.N(1) + 1;
            xzeroVec = grid_obj.minv(1) * ones(1, xlen);
            xoneVec  = grid_obj.maxv(1) * ones(1, xlen);

            if strcmp(uq_obj.method, 'NI')
                
                s = uq_obj.GridData.grid;

                obj.cA = @(t,u)(cAfun(obj, grid_obj, sbp_obj, uq_obj, problem_obj, X, xzeroVec, xoneVec, t, u, s));

            elseif strcmp(uq_obj.method, 'PC')

                obj.cA = @(t,u)(cAfun(obj, grid_obj, sbp_obj, uq_obj, problem_obj, X, xzeroVec, xoneVec, t, u));

            elseif strcmp(uq_obj.method, 'DM')

                obj.cA = @(t,u)(cAfun(obj, grid_obj, sbp_obj, uq_obj, problem_obj, X, xzeroVec, xoneVec, t, u, 0));

            end
        end
        function obj = createPenalties(obj)
            
            obj.Sigma.E =  obj.Boundary_Operators_n.E' * obj.Lambda.E_m;
            obj.Sigma.W = -obj.Boundary_Operators_n.W' * obj.Lambda.W_p;

        end
        function obj = createBoundaryOperator(obj, problem_obj, u)
            
            % Diagonalization A
            dim  = size(problem_obj.disc.A, 3);
            l    = length(u);
            Am   = sparse(l, l);
            
            for k = 1:dim

                Am   = Am + kron(diag(sparse(u(k:dim:end))), sparse(problem_obj.disc.A(:,:,k)));
               
            end

%             [~, obj.Lambda.A]   = eig(full(Am));
%             [lambdaA, index]    = sort(diag(obj.Lambda.A));
%             %obj.Lambda.A        = flip(flip(sparse(obj.Lambda.A))');
%             posDimA             = sum(lambdaA > 1e-12);
%             negDimA             = sum(lambdaA < -1e-12);
%             zeroDimA            = sum(abs(lambdaA) < 1e-12);
%             Aplus               = sparse(lambdaA(end - posDimA + 1:end));
%             Aminus              = sparse(lambdaA(1:negDimA));
%             Azero               = sparse(zeros(zeroDimA));
%             obj.Lambda.A        = sparse(blkdiag(diag(Aplus), Azero, diag(Aminus)));
% 
%             obj.Lambda.W_p = 2/3 * sparse(blkdiag(diag(Aplus) , sparse(l - length(Aplus), l - length(Aplus))));
%             obj.Lambda.E_m = 2/3 * sparse(blkdiag(sparse(l - length(Aminus), l - length(Aminus)), diag(Aminus)));
%             
%             obj.Boundary_Operators_n.E = sparse(blkdiag(sparse(posDimA, posDimA), sparse(zeroDimA, zeroDimA), speye(negDimA)));
%             obj.Boundary_Operators_n.W = sparse(blkdiag(speye(posDimA, posDimA), sparse(zeroDimA, zeroDimA), sparse(negDimA, negDimA)));
%             
%             obj.Scheme.Am = sparse(obj.Lambda.A);
%

            dim = problem_obj.dim;
            Am_W = full(Am(1:dim, 1:dim));
            [X_W, obj.Lambda.A_W] = eig(Am_W);
            [lambdaA_W, index_W]  = sort(diag(obj.Lambda.A_W));
            X_W            = sparse(X_W(:, index_W));
            posDimA_W      = sum(lambdaA_W > 0);
            negDimA_W      = sum(lambdaA_W < 0);
            zeroDimA_W     = sum(abs(lambdaA_W) == 0);
            Aplus_W        = lambdaA_W(end - posDimA_W + 1:end);
            Aminus_W       = lambdaA_W(1:negDimA_W);
            Azero_W        = sparse(zeroDimA_W, zeroDimA_W);
            obj.Lambda.A_W = blkdiag(diag(Aplus_W), Azero_W, diag(Aminus_W));
            X_Wplus        = X_W(:, end-posDimA_W + 1:end);
            obj.Lambda.W_p = sparse(kron(obj.Matrices.I, obj.Lambda.A_W));
            
            
            Am_E = full(Am(end-dim + 1:end, end-dim + 1:end));
            [X_E, obj.Lambda.A_E] = eig(Am_E);
            [lambdaA_E, index_E]  = sort(diag(obj.Lambda.A_E));
            X_E                   = sparse(X_E(:, index_E));
            obj.Lambda.A_E        = sparse(lambdaA_E);
            posDimA_E             = sum(lambdaA_E > 0);
            negDimA_E             = sum(lambdaA_E < 0);
            zeroDimA_E            = sum(abs(lambdaA_E) == 0);
            Aplus_E               = lambdaA_E(end - posDimA_E + 1:end);
            Aminus_E              = lambdaA_E(1:negDimA_E);
            Azero_E               = sparse(zeroDimA_E, zeroDimA_E);
            obj.Lambda.A_E        = blkdiag(diag(Aplus_E), Azero_E, diag(Aminus_E));
            X_Eminus              = X_E(:, 1:negDimA_E);
            obj.Lambda.E_m        = sparse(kron(obj.Matrices.I, obj.Lambda.A_E));
            
            obj.Boundary_Operators_n.W = 2/3*sparse(kron(obj.Matrices.I, [X_Wplus, zeros(dim, negDimA_W + zeroDimA_W)]'));
            obj.Boundary_Operators_n.E = 2/3*sparse(kron(obj.Matrices.I, [zeros(dim, posDimA_E + zeroDimA_E), X_Eminus]'));
            
            obj.Scheme.Am = sparse(Am);
            
        end
        function obj = solveUsingRK4(obj, grid_obj, problem_obj, uq_obj)
            
            
            X  = grid_obj.grid;
            opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

            if strcmp(uq_obj.method, 'NI')
                
                obj.Solution.u_ns = zeros(length(X), length(grid_obj.t), length(uq_obj.GridData.grid));
            
                s = uq_obj.GridData.grid';

                [~, obj.Solution.u_n] = ode45(obj.cA, grid_obj.t, ...
                    reshape(problem_obj.data.f(X, 0, s), [], 1), opts);

                obj.Solution.u_ns = zeros(length(X), length(grid_obj.t), problem_obj.dim);
                obj.Solution.u_n  = obj.Solution.u_n';
                
                for k = 1:problem_obj.dim
                   
                    obj.Solution.u_ns(:,:,k) = obj.Solution.u_n(k:problem_obj.dim:end, :);
                    
                end                

            elseif strcmp(uq_obj.method, 'PC')
                
                [~, obj.Solution.u_n] = ode45(obj.cA, grid_obj.t, ...
                        reshape(obj.innerproduct(uq_obj, problem_obj.data.f, X, 0), [], 1), opts);
                
                obj.Solution.u_ns = zeros(length(X), length(grid_obj.t), problem_obj.dim);
                obj.Solution.u_n  = obj.Solution.u_n';
                
                for k = 1:problem_obj.dim
                   
                    obj.Solution.u_ns(:,:,k) = obj.Solution.u_n(k:problem_obj.dim:end, :);
                    
                end
            
            elseif strcmp(uq_obj.method, 'DM')
                
                [~, obj.Solution.u_n] = ode45(obj.cA, grid_obj.t, ...
                        reshape(problem_obj.data.f(X, 0, 0), [], 1), opts);
                    
                obj.Solution.u_n = obj.Solution.u_n';
                
            end
        end
        function obj = transformation(obj, grid_obj, problem_obj)
            
            Ix = eye(grid_obj.N(1) + 1);

            obj.Trans.A    = sparse(kron(Ix, problem_obj.disc.A));
            obj.Trans.B    = sparse(kron(Ix, problem_obj.disc.B));
            obj.Trans.C    = sparse(kron(Ix, problem_obj.disc.C));
            obj.Trans.D    = sparse(kron(Ix, problem_obj.disc.D));
            obj.Trans.F    = sparse(kron(Ix, problem_obj.disc.F));
            obj.Trans.G    = sparse(kron(Ix, problem_obj.disc.G));
            obj.Trans.ADxi = sparse(kron(Ix, problem_obj.disc.Ax));
            obj.Trans.BDyi = sparse(kron(Ix, problem_obj.disc.By));
            
        end
        function plotSolution(obj, grid_obj, problem_obj)
            
            for k = 1:problem_obj.dim
                
                figure
                plot(grid_obj.grid, obj.Solution.u_n(k:problem_obj.dim:end - problem_obj.dim + k, 1))
                
            end
            
        end
        function obj = computeError(obj, grid_obj, sbp_obj, problem_obj)
            
            obj   = computeAnalyticalSolution(obj, grid_obj, problem_obj, grid_obj.maxv(end));
%            error = abs(obj.Solution.u_n(:, end) - obj.Solution.u_a);
%            obj.Solution.e_L2 = sqrt(error' * (kron(sbp_obj.Pinv, obj.Matrices.Isys)  \ error));
            
        end
        function obj = computeAnalyticalSolution(obj, grid_obj, problem_obj, t)
            
            X = grid_obj.grid;
            obj.Solution.u_a = reshape(problem_obj.data.u_a(X, t, 0), [], 1);
                
        end
        function c   = innerproduct(~, uq_obj, f, x, t)
            
            weight = NodesWeights2(uq_obj, uq_obj.Distribution.name, uq_obj.MaxPol + 1);
            r = uq_obj.IPfact * f(x, t, weight(:, 2));
            c = r(:);
            
        end
        function cA  = cAfun(obj, grid_obj, sbp_obj, uq_obj, problem_obj, X, xzeroVec, xoneVec, t, u, s)
            
            % Boundary conditions
            obj = createBoundaryOperator(obj, problem_obj, u);
            
            % Penalty coefficients
            obj = createPenalties(obj);
            
            % Create A matrix
            obj = createAMatrix(obj, grid_obj, sbp_obj, problem_obj, u);
            
            if strcmp(uq_obj.method, 'PC')
                
                force = @(t) (        reshape(obj.innerproduct(uq_obj, problem_obj.data.uT_a,  X, t), [], 1) ...
                    + 1/3*obj.Scheme.Am * reshape(obj.innerproduct(uq_obj, problem_obj.data.uX_a,  X, t), [], 1) ...
                    + 1/3*obj.Matrices.Dx * obj.Scheme.Am * reshape(obj.innerproduct(uq_obj, problem_obj.data.u_a,  X, t), [], 1));
            
                cAtemp = @(t, u) ( obj.A * u ...
                   - obj.Scheme.PinvE0SigmaWH0 * reshape(obj.innerproduct(uq_obj, problem_obj.data.gW,  xzeroVec, t), [], 1) ...
                   - obj.Scheme.PinvENSigmaEH0 * reshape(obj.innerproduct(uq_obj, problem_obj.data.gE,  xoneVec,  t), [], 1) ...
                   + force(t));
               
            elseif strcmp(uq_obj.method, 'NI')
            
                force = @(t) (        reshape(problem_obj.data.uT_a(X, t, s'), [], 1) ...
                    + obj.Scheme.Am * reshape(problem_obj.data.uX_a(X, t, s'), [], 1));
            
                cAtemp = @(t, u) ( obj.A * u ...
                    - obj.Scheme.PinvE0SigmaWH0 * reshape(problem_obj.data.gW(xzeroVec, t, s'), [], 1) ...
                    - obj.Scheme.PinvENSigmaEH0 * reshape(problem_obj.data.gE(xoneVec,  t, s'), [], 1) ...
                    + force(t));
                
            elseif strcmp(uq_obj.method, 'DM')

                force = @(t) (        reshape(problem_obj.data.uT_a(X, t, 0), [], 1) ...
                    + obj.Scheme.Am * reshape(problem_obj.data.uX_a(X, t, 0), [], 1));

                cAtemp = @(t, u) ( obj.A * u ...
                    - obj.Scheme.PinvE0SigmaWH0 * reshape(problem_obj.data.gW(xzeroVec, t, 0), [], 1) ...
                    - obj.Scheme.PinvENSigmaEH0 * reshape(problem_obj.data.gE(xoneVec,  t, 0), [], 1) ...
                    + force(t));
                
            end
            
            cA = cAtemp(t, u);
            
        end
    end
end
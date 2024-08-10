classdef problem_obj
    % problem object
    % ut + Aux + Buy = e(Cuxx + Duyy)
    %       Hu = g
    %        u = f
    
    properties
        % Parameters
        cont;  % Structure with continuous matrices: a, ax, b, by, c, d
        disc;  % Structure with discrete matrices: A, Ax, B, By, C, D
        dim;   % dimension of system
        
        % Problem type
        type;  % Type: problem type, data type, coeff type
        
        % Transformation
        Trans  % Transformation:
        
        % Data
        data;  % Boundary data: gE, gED, gW, gWD, ...
        
    end
    methods
        function obj = initialize(obj)
            
            obj = initializeProblem(obj);
            obj = initializeData(obj);
            
        end
        function obj = initializeProblem(obj)
            
            if strcmp(obj.type.problem, 'Burger')
                
                obj.dim           = 1;
                obj.type.coeff    = 'constant';
                obj.type.linear   = 'non-linear';
                obj.type.dissipation = 'None';
                obj.cont.a        = 1;
                
            elseif strcmp(obj.type.problem, 'BurgersSystem')
                
                obj.dim           = 3;
                obj.type.coeff    = 'constant';
                obj.type.linear   = 'non-linear';
                obj.type.uq       = 'deterministic';
                obj.cont.a(:,:,1) = [2 0 3; 0 2 0; 3 0 1];
                obj.cont.a(:,:,2) = [0 1 0; 1 2 0; 0 0 0];
                obj.cont.a(:,:,3) = [1 0 3; 0 2 0; 3 0 2];
                
            elseif strcmp(obj.type.problem, 'BurgersSystemConstant')
                
                obj.dim           = 3;
                obj.type.coeff    = 'constant';
                obj.type.linear   = 'non-linear';
                obj.type.uq       = 'deterministic';
                obj.cont.a(:,:,1) = [1 0 0; 0 0 0; 0 0 0];
                obj.cont.a(:,:,2) = [0 0 0; 0 1 0; 0 0 0];
                obj.cont.a(:,:,3) = [0 0 0; 0 0 0; 0 0 -1];
               
            elseif strcmp(obj.type.problem, 'BurgersSystemPC6')
                
                M = 5;
                obj.cont.a = zeros(M + 1, M + 1, M + 1);
                obj.dim           = 6;
                obj.type.coeff    = 'constant';
                obj.type.linear   = 'non-linear';
                obj.type.uq       = 'deterministic';
                obj.type.dissipation = 'dissipation';
                
                for i = 0:M

                    for j = 0:M

                        for k = 0:M

                            s = (i + j + k) / 2;

                            if rem(i + j + k, 2) == 1 || abs(i - j) > k || k > i + j || max([i,j,k]) > s

                                obj.cont.a(i + 1, j + 1, k + 1) = 0;

                            else

                                obj.cont.a(i + 1, j + 1, k + 1) = sqrt((2 * i + 1) ...
                                    * (2 * j + 1) * (2 * k + 1)) / (i + j + k + 1) ...
                                    * obj.A_L(s - i) * obj.A_L(s - j) * obj.A_L(s - k) / obj.A_L(s);

%                                 obj.cont.a(i + 1, j + 1, k + 1) = ...
%                                     factorial(i) * factorial(j) * factorial(k) / ...
%                                     (factorial(s - i) * factorial(s - j) * factorial(s - k));

                            end
                        end
                    end
                end
                
                elseif strcmp(obj.type.problem, 'BurgersSystemPC7')
                
                M = 6;
                obj.cont.a           = zeros(M + 1, M + 1, M + 1);
                obj.dim              = 7;
                obj.type.coeff       = 'constant';
                obj.type.linear      = 'non-linear';
                obj.type.uq          = 'deterministic';
                obj.type.dissipation = 'dissipation';
                for i = 0:M

                    for j = 0:M

                        for k = 0:M

                            s = (i + j + k) / 2;

                            if rem(i + j + k, 2) == 1 || abs(i - j) > k || k > i + j || max([i,j,k]) > s

                                obj.cont.a(i + 1, j + 1, k + 1) = 0;

                            else

                                obj.cont.a(i + 1, j + 1, k + 1) = sqrt((2 * i + 1) ...
                                    * (2 * j + 1) * (2 * k + 1)) / (i + j + k + 1) ...
                                    * obj.A_L(s - i) * obj.A_L(s - j) * obj.A_L(s - k) / obj.A_L(s);

%                                 obj.cont.a(i + 1, j + 1, k + 1) = ...
%                                     factorial(i) * factorial(j) * factorial(k) / ...
%                                     (factorial(s - i) * factorial(s - j) * factorial(s - k));

                            end
                        end
                    end
                end
                %obj.cont.a = obj.cont.a > 0;
                
            end
        end
        function obj = initializeData(obj)
            
            if strcmp(obj.type.data, 'Sin')
                
                obj.data.u_a   = @(x, t, s) (          sin(2*pi*(x - t)) + 0*x + 0*t + 0*s);
                obj.data.uT_a  = @(x, t, s) (-2 * pi * cos(2*pi*(x - t)) + 0*x + 0*t + 0*s);
                obj.data.uX_a  = @(x, t, s) ( 2 * pi * cos(2*pi*(x - t)) + 0*x + 0*t + 0*s);
                
            elseif strcmp(obj.type.data, 'constant')
                
                obj.data.u_a   = @(x, t, s) (1 + 0*x + 0*t + 0*s);
                obj.data.uT_a  = @(x, t, s) (0 + 0*x + 0*t + 0*s);
                obj.data.uX_a  = @(x, t, s) (0 + 0*x + 0*t + 0*s);
                
            elseif strcmp(obj.type.data, 'Constant')
            
                obj.data.u_a   = @(x, t, s) (3 + 0*x + 0*t + 0*s);
                obj.data.uT_a  = @(x, t, s) (2 + 0*x + 0*t + 0*s);
                obj.data.uX_a  = @(x, t, s) (4 + 0*x + 0*t + 0*s);
                
            elseif strcmp(obj.type.data, 'BurgersSystem')
                
                obj.data.u_a   = @(x, t, s) ([ ...
                    sin(2*pi*(x - t)) + 0*t + 0*s; ...
                    sin(2*pi*(x - t)) + 0*t + 0*s; ...
                    sin(2*pi*(x - t)) + 0*t + 0*s]);
                obj.data.uT_a  = @(x, t, s) ([ ...
                    -2 * pi * cos(2*pi*(x - t)) + 0*x + 0*t + 0*s; ...
                    -2 * pi * cos(2*pi*(x - t)) + 0*x + 0*t + 0*s; ...
                    -2 * pi * cos(2*pi*(x - t)) + 0*x + 0*t + 0*s]);
                obj.data.uX_a  = @(x, t, s) ([ ...
                    2 * pi * cos(2*pi*(x - t)) + 0*x + 0*t + 0*s; ...
                    2 * pi * cos(2*pi*(x - t)) + 0*x + 0*t + 0*s; ...
                    2 * pi * cos(2*pi*(x - t)) + 0*x + 0*t + 0*s]);
            elseif strcmp(obj.type.data, 'BurgersSystemPC6')
                
                obj.data.u_a   = @(x, t, s) ([ ...
                    sin(2*pi*(x - t)) + 0*t + 0*s; ...
                    sin(2*pi*(x - t)) + 0*t + 0*s; ...
                    sin(2*pi*(x - t)) + 0*t + 0*s; ...
                    sin(2*pi*(x - t)) + 0*t + 0*s; ...
                    sin(2*pi*(x - t)) + 0*t + 0*s; ...
                    sin(2*pi*(x - t)) + 0*t + 0*s]);
                obj.data.uT_a  = @(x, t, s) ([ ...
                    -2 * pi * cos(2*pi*(x - t)) + 0*x + 0*t + 0*s; ...
                    -2 * pi * cos(2*pi*(x - t)) + 0*x + 0*t + 0*s; ...
                    -2 * pi * cos(2*pi*(x - t)) + 0*x + 0*t + 0*s; ...
                    -2 * pi * cos(2*pi*(x - t)) + 0*x + 0*t + 0*s; ...
                    -2 * pi * cos(2*pi*(x - t)) + 0*x + 0*t + 0*s; ...
                    -2 * pi * cos(2*pi*(x - t)) + 0*x + 0*t + 0*s]);
                obj.data.uX_a  = @(x, t, s) ([ ...
                    2 * pi * cos(2*pi*(x - t)) + 0*x + 0*t + 0*s; ...
                    2 * pi * cos(2*pi*(x - t)) + 0*x + 0*t + 0*s; ...
                    2 * pi * cos(2*pi*(x - t)) + 0*x + 0*t + 0*s; ...
                    2 * pi * cos(2*pi*(x - t)) + 0*x + 0*t + 0*s; ...
                    2 * pi * cos(2*pi*(x - t)) + 0*x + 0*t + 0*s; ...
                    2 * pi * cos(2*pi*(x - t)) + 0*x + 0*t + 0*s]);
                 
            elseif strcmp(obj.type.data, 'BurgersSystemPC7')
                
                obj.data.u_a   = @(x, t, s) ([ ...
                    sin(2*pi*(x - t)) + 0*t + 0*s; ...
                    sin(2*pi*(x - t)) + 0*t + 0*s; ...
                    sin(2*pi*(x - t)) + 0*t + 0*s; ...
                    sin(2*pi*(x - t)) + 0*t + 0*s; ...
                    sin(2*pi*(x - t)) + 0*t + 0*s; ...
                    sin(2*pi*(x - t)) + 0*t + 0*s; ...
                    sin(2*pi*(x - t)) + 0*t + 0*s]);
                obj.data.uT_a  = @(x, t, s) ([ ...
                    -2 * pi * cos(2*pi*(x - t)) + 0*x + 0*t + 0*s; ...
                    -2 * pi * cos(2*pi*(x - t)) + 0*x + 0*t + 0*s; ...
                    -2 * pi * cos(2*pi*(x - t)) + 0*x + 0*t + 0*s; ...
                    -2 * pi * cos(2*pi*(x - t)) + 0*x + 0*t + 0*s; ...
                    -2 * pi * cos(2*pi*(x - t)) + 0*x + 0*t + 0*s; ...
                    -2 * pi * cos(2*pi*(x - t)) + 0*x + 0*t + 0*s; ...
                    -2 * pi * cos(2*pi*(x - t)) + 0*x + 0*t + 0*s]);
                obj.data.uX_a  = @(x, t, s) ([ ...
                    2 * pi * cos(2*pi*(x - t)) + 0*x + 0*t + 0*s; ...
                    2 * pi * cos(2*pi*(x - t)) + 0*x + 0*t + 0*s; ...
                    2 * pi * cos(2*pi*(x - t)) + 0*x + 0*t + 0*s; ...
                    2 * pi * cos(2*pi*(x - t)) + 0*x + 0*t + 0*s; ...
                    2 * pi * cos(2*pi*(x - t)) + 0*x + 0*t + 0*s; ...
                    2 * pi * cos(2*pi*(x - t)) + 0*x + 0*t + 0*s; ...
                    2 * pi * cos(2*pi*(x - t)) + 0*x + 0*t + 0*s]);
            elseif strcmp(obj.type.data, 'BurgersSystemPC7Constant')
                
                obj.data.u_a   = @(x, t, s) ([ ...
                    sin(2*pi*(x - t)) + 0*t + 0*s; ...
                    sin(2*pi*(x - t)) + 0*t + 0*s; ...
                    sin(2*pi*(x - t)) + 0*t + 0*s; ...
                    sin(2*pi*(x - t)) + 0*t + 0*s; ...
                    sin(2*pi*(x - t)) + 0*t + 0*s; ...
                    sin(2*pi*(x - t)) + 0*t + 0*s; ...
                    sin(2*pi*(x - t)) + 0*t + 0*s]);
                obj.data.uT_a  = @(x, t, s) ([ ...
                    -2 * pi * cos(2*pi*(x - t)) + 0*x + 0*t + 0*s; ...
                    -2 * pi * cos(2*pi*(x - t)) + 0*x + 0*t + 0*s; ...
                    -2 * pi * cos(2*pi*(x - t)) + 0*x + 0*t + 0*s; ...
                    -2 * pi * cos(2*pi*(x - t)) + 0*x + 0*t + 0*s; ...
                    -2 * pi * cos(2*pi*(x - t)) + 0*x + 0*t + 0*s; ...
                    -2 * pi * cos(2*pi*(x - t)) + 0*x + 0*t + 0*s; ...
                    -2 * pi * cos(2*pi*(x - t)) + 0*x + 0*t + 0*s]);
                obj.data.uX_a  = @(x, t, s) ([ ...
                    2 * pi * cos(2*pi*(x - t)) + 0*x + 0*t + 0*s; ...
                    2 * pi * cos(2*pi*(x - t)) + 0*x + 0*t + 0*s; ...
                    2 * pi * cos(2*pi*(x - t)) + 0*x + 0*t + 0*s; ...
                    2 * pi * cos(2*pi*(x - t)) + 0*x + 0*t + 0*s; ...
                    2 * pi * cos(2*pi*(x - t)) + 0*x + 0*t + 0*s; ...
                    2 * pi * cos(2*pi*(x - t)) + 0*x + 0*t + 0*s; ...
                    2 * pi * cos(2*pi*(x - t)) + 0*x + 0*t + 0*s]);
                
            elseif strcmp(obj.type.data, 'BurgersSystemConstant')
                
                obj.data.u_a   = @(x, t, s) ([ ...
                    1 + 0*x + 0*t + 0*s; ...
                    1 + 0*x + 0*t + 0*s; ...
                    1 + 0*x + 0*t + 0*s]);
                obj.data.uT_a  = @(x, t, s) ([ ...
                    0 + 0*x + 0*t + 0*s; ...
                    0 + 0*x + 0*t + 0*s; ...
                    0 + 0*x + 0*t + 0*s]);
                obj.data.uX_a  = @(x, t, s) ([ ...
                    0 + 0*x + 0*t + 0*s; ...
                    0 + 0*x + 0*t + 0*s; ...
                    0 + 0*x + 0*t + 0*s]);
                
            elseif strcmp(obj.type.data, 'sincos_smooth')
                
                obj.data.u_a   = @(x, t, s) ( ...
                    sin(2*pi*(x - t)) + s .* sin(s/2 - t) .* cos(x - s));
                obj.data.uT_a  = @(x, t, s) ( ...
                    -2*pi*cos(2*pi*(x - t)) - s .* cos(s/2 - t) .* cos(x - s));
                obj.data.uX_a  = @(x, t, s) ( ...
                    2*pi*cos(2*pi*(x - t)) - s .* sin(s/2 - t) .* sin(x - s));

            elseif strcmp(obj.type.data, 'sincos_nonsmooth')
                
                obj.data.u_a   = @(x, t, s) ( ...
                    sin(2*pi*(x - t)) +   sin(5*(s/2 - t)) .* cos(2*(x - 5*s)));
                obj.data.uT_a  = @(x, t, s) ( ...
                    -2*pi*cos(2*pi*(x - t)) -   5*cos(5*(s/2 - t)) .* cos(2*(x - 5*s)));
                obj.data.uX_a  = @(x, t, s) ( ...
                    2*pi*cos(2*pi*(x - t)) -   2*sin(5*(s/2 - t)) .* sin(2*(x - 5*s)));
                
            end
        end
        function obj = createMatrices(obj, grid_obj, uq_obj)
                            
            if strcmp(obj.type.coeff, 'constant')

                if strcmp(obj.type.uq, 'deterministic') || strcmp(obj.type.uq, 'NI')

                    obj.disc.A = obj.cont.a;
                    obj.dim    = length(uq_obj.GridData.grid);
                    
                elseif strcmp(obj.type.uq, 'PC')

                    obj.disc.A = uq_obj.TripleProduct(obj.cont.a);
                    obj.disc.S = uq_obj.DoubleProduct();
                    obj.dim    = size(obj.disc.A, 3);
                    
                end

            elseif strcmp(obj.type.coeff, 'variable')

                avals      = obj.cont.a (grid_obj.grid(:, :, 1))';
                obj.disc.A = sparse(diag(avals (:)));

            end
            
            obj = transformation(obj);
            
        end
        function obj = createData(obj)

            obj.data.gE  = @(x, t, s) (obj.data.u_a (x, t, s));
            obj.data.gED = @(x, t, s) (obj.data.uX_a(x, t, s));
            obj.data.gW  = @(x, t, s) (obj.data.u_a (x, t, s));
            obj.data.gWD = @(x, t, s) (obj.data.uX_a(x, t, s));
            obj.data.f   = @(x, t, s) (obj.data.u_a (x, t, s));

        end
        function obj = transformation(obj)
                            
            if strcmp(obj.type.coeff, 'variable')

                obj.Trans.A_a = sparse(obj.disc.A);

            elseif strcmp(obj.type.coeff, 'constant')

                obj.Trans.A_a = obj.disc.A;

            end
        end
        function A        = A_L(~, n)
           
            if n == 0
                
                A = 1;
                
            elseif n < 0
                
                A = 0;
                
            elseif n >= 1
                
                A = 1;
                
                for j = 1:n
                    
                    A = A * (2 * j - 1);
                    
                end
                
                A = A / factorial(n);
                
            end
        end
    end
end
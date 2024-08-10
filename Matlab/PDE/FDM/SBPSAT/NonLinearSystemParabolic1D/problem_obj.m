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
            
            if strcmp(obj.type.problem, 'Burger001')
                
                obj.dim           = 1;
                obj.type.coeff    = 'constant';
                obj.type.linear   = 'non-linear';
                obj.type.dissipation = 'None';
                obj.cont.a        = 1;
                obj.cont.b        = .01;
            elseif strcmp(obj.type.problem, 'Burger01')
                
                obj.dim           = 1;
                obj.type.coeff    = 'constant';
                obj.type.linear   = 'non-linear';
                obj.type.dissipation = 'None';
                obj.cont.a        = 1;
                obj.cont.b        = .1;
                
            elseif strcmp(obj.type.problem, 'Burger1')
                
                obj.dim              = 1;
                obj.type.coeff       = 'constant';
                obj.type.linear      = 'non-linear';
                obj.type.dissipation = 'None';
                obj.cont.a           = 1;
                obj.cont.b           = 1;
                
            elseif strcmp(obj.type.problem, 'BurgersSystem')
                
                obj.dim           = 3;
                obj.type.coeff    = 'constant';
                obj.type.linear   = 'non-linear';
                obj.type.uq       = 'deterministic';
                obj.cont.a(:,:,1) = [2 0 3; 0 2 0; 3 0 1];
                obj.cont.a(:,:,2) = [0 1 0; 1 2 0; 0 0 0];
                obj.cont.a(:,:,3) = [1 0 3; 0 2 0; 3 0 2];
                obj.cont.b        = 0.1 * eye(3);
                
            elseif strcmp(obj.type.problem, 'BurgersSystemConstant')
                
                obj.dim           = 3;
                obj.type.coeff    = 'constant';
                obj.type.linear   = 'non-linear';
                obj.type.uq       = 'deterministic';
                obj.cont.a(:,:,1) = [1 0 0; 0 0 0; 0 0 0];
                obj.cont.a(:,:,2) = [0 0 0; 0 1 0; 0 0 0];
                obj.cont.a(:,:,3) = [0 0 0; 0 0 0; 0 0 -1];
                obj.cont.b        = 0.1 * eye(3);
                
            elseif strcmp(obj.type.problem, 'ViscousBurger')
                
                obj.dim           = 1;
                obj.type.coeff    = 'constant';
                obj.type.linear   = 'non-linear';
                obj.type.dissipation = 'None';
                obj.cont.a        = 0;
                obj.cont.b        = 0.01;
                
            end
        end
        function obj = initializeData(obj)
            
            if strcmp(obj.type.data, 'Sin')
                
                obj.data.u_a   = @(x, t, s) (10 +              sin(2*pi*(x - t)) + 0*x + 0*t + 0*s);
                obj.data.uT_a  = @(x, t, s) (- 2 * pi    * cos(2*pi*(x - t)) + 0*x + 0*t + 0*s);
                obj.data.uX_a  = @(x, t, s) (  2 * pi    * cos(2*pi*(x - t)) + 0*x + 0*t + 0*s);
                obj.data.uXX_a = @(x, t, s) (-(2 * pi)^2 * sin(2*pi*(x - t)) + 0*x + 0*t + 0*s);
                
            elseif strcmp(obj.type.data, 'Sin2')
                
                obj.data.u_a   = @(x, t, s) (10+             sin(2*pi*(x)) + 0*x + 0*t + 0*s);
                obj.data.uT_a  = @(x, t, s) ( 0*x + 0*t + 0*s);
                obj.data.uX_a  = @(x, t, s) (  2 * pi    * cos(2*pi*(x)) + 0*x + 0*t + 0*s);
                obj.data.uXX_a = @(x, t, s) (-(2 * pi)^2 * sin(2*pi*(x)) + 0*x + 0*t + 0*s);
                    
            elseif strcmp(obj.type.data, 'zero')
                
                obj.data.u_a   = @(x, t, s) (0.01 + 1*x + 0*t + 0*s);
                obj.data.uT_a  = @(x, t, s) (0 + 0*x + 0*t + 0*s);
                obj.data.uX_a  = @(x, t, s) (1 + 0*x + 0*t + 0*s);
                obj.data.uXX_a = @(x, t, s) (0 + 0*x + 0*t + 0*s);
                    
            elseif strcmp(obj.type.data, 'constant')
                
                obj.data.u_a   = @(x, t, s) (-10 + 1*x.^2 + 0*t + 0*s);
                obj.data.uT_a  = @(x, t, s) (0 + 0*x + 0*t + 0*s);
                obj.data.uX_a  = @(x, t, s) (0 + 2*x + 0*t + 0*s);
                obj.data.uXX_a = @(x, t, s) (2 + 0*x + 0*t + 0*s);
                
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
                obj.data.uXX_a = @(x, t, s) ( ...
                    -(2*pi)^2*sin(2*pi*(x - t)) - s .* sin(s/2 - t) .* cos(x - s));

            elseif strcmp(obj.type.data, 'sincos_nonsmooth')
                
                obj.data.u_a   = @(x, t, s) ( ...
                    sin(2*pi*(x - t)) +   sin(5*(s/2 - t)) .* cos(2*(x - 5*s)));
                obj.data.uT_a  = @(x, t, s) ( ...
                    -2*pi*cos(2*pi*(x - t)) -   5*cos(5*(s/2 - t)) .* cos(2*(x - 5*s)));
                obj.data.uX_a  = @(x, t, s) ( ...
                    2*pi*cos(2*pi*(x - t)) -   2*sin(5*(s/2 - t)) .* sin(2*(x - 5*s)));
                obj.data.uXX_a  = @(x, t, s) ( ...
                    -(2*pi)^2*sin(2*pi*(x - t)) - 4*sin(5*(s/2 - t)) .* cos(2*(x - 5*s)));
                
            elseif strcmp(obj.type.data, 'Exact')
                
                obj.data.u_a   = @(x, t, s) ( ...
                    exp(-1*abs((x).^2)) + sin(2*pi*s) + 0*t + 0*s + 0*x);
                obj.data.uT_a  = @(x, t, s) ( ...
                    0*x + 0*t + 0*s);
                obj.data.uX_a  = @(x, t, s) ( ...
                    0*x + 0*t + 0*s);
                obj.data.uXX_a  = @(x, t, s) ( ...
                    0*x + 0*t + 0*s);
                
            elseif strcmp(obj.type.data, 'Heaviside')
                
                obj.data.u_a   = @(x, t, s) ( ...
                    sin(2*pi*(x - t)) + 1./(1 + exp(-1*s/(20*obj.cont.b))) + 2);
                obj.data.uT_a  = @(x, t, s) ( ...
                    -2*pi*cos(2*pi*(x - t)) + 0*s);
                obj.data.uX_a  = @(x, t, s) ( ...
                    2*pi*cos(2*pi*(x - t)) + 0*s);
                obj.data.uXX_a  = @(x, t, s) ( ...
                    -(2*pi)^2*sin(2*pi*(x - t)) + 0*s);
                
            elseif strcmp(obj.type.data, 'Sin_fast')
                
                a = @(s) (1*pi*exp(-.5*s.^2));
                b = @(s) (exp(-s/2));
                
                obj.data.u_a   = @(x, t, s) ( ...
                    5 + b(s).*sin(a(s).*x - t));
                obj.data.uT_a  = @(x, t, s) ( ...
                    -b(s).*cos(a(s).*x - t));
                obj.data.uX_a  = @(x, t, s) ( ...
                    a(s).*b(s).*cos(a(s).*x - t));
                obj.data.uXX_a  = @(x, t, s) ( ...
                    -a(s).^2.*b(s).*sin(a(s).*x - t));
                
            elseif strcmp(obj.type.data, 'Sin_slow')
                
                a = @(s) (5*pi*exp(-0.5*s.^2));
                b = @(s) (exp(-s/2));
                
                obj.data.u_a   = @(x, t, s) ( ...
                    5 + b(s).*sin(a(s).*x - t));
                obj.data.uT_a  = @(x, t, s) ( ...
                    -b(s).*cos(a(s).*x - t));
                obj.data.uX_a  = @(x, t, s) ( ...
                    a(s).*b(s).*cos(a(s).*x - t));
                obj.data.uXX_a  = @(x, t, s) ( ...
                    -a(s).^2.*b(s).*sin(a(s).*x - t));
            end
        end
        function obj = createMatrices(obj, grid_obj, uq_obj)
                            
            if strcmp(obj.type.coeff, 'constant')

                if strcmp(obj.type.uq, 'deterministic') || strcmp(obj.type.uq, 'NI')

                    obj.disc.A = obj.cont.a;
                    obj.dim    = length(uq_obj.GridData.grid);
                    
                elseif strcmp(obj.type.uq, 'PC')

                    obj.disc.A = uq_obj.TripleProduct(obj.cont.a);
                    B = uq_obj.TripleHermite(obj.cont.a);
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

            if strcmp(obj.type.data, 'Exact')
                
                obj.data.gE  = @(x, t, s) (0*x + 0*t + 0*s);
                obj.data.gED = @(x, t, s) (0*x + 0*t + 0*s);
                obj.data.gW  = @(x, t, s) (0*x + 0*t + 0*s);
                obj.data.gWD = @(x, t, s) (0*x + 0*t + 0*s);
                obj.data.f   = @(x, t, s) (obj.data.u_a (x, t, s));
            
            else
                
                obj.data.gE  = @(x, t, s) (obj.data.u_a (x, t, s));
                obj.data.gED = @(x, t, s) (obj.data.uX_a(x, t, s));
                obj.data.gW  = @(x, t, s) (obj.data.u_a (x, t, s));
                obj.data.gWD = @(x, t, s) (obj.data.uX_a(x, t, s));
                obj.data.f   = @(x, t, s) (obj.data.u_a (x, t, s));
                
            end
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
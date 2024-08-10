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
            
           if strcmp(obj.type.problem, 'Advection3D')
                
                obj.cont.a     = @(x, y, z) ( sin(x) .* cos(y));
                obj.cont.ax    = @(x, y, z) ( cos(x) .* cos(y));
                obj.cont.b     = @(x, y, z) (-cos(x) .* sin(y));
                obj.cont.by    = @(x, y, z) (-cos(x) .* cos(y));
                
            elseif strcmp(obj.type.problem, 'Scalar')
                
                obj.cont.a     = 2;
                obj.cont.b     = 3;
                obj.cont.c     = 4;

            end
        end
        function obj = initializeData(obj)
            
            if strcmp(obj.type.data, 'Sin')
                
                obj.data.u_a   = @(x, y, z, t) (            sin(2*pi*(x - t)) + sin(2*pi*(y - t)) + sin(2*pi*(z - t))     );
                obj.data.uX_a  = @(x, y, z, t) (  2*pi    * cos(2*pi*(x - t)) + 0*y + 0*z                   );
                obj.data.uY_a  = @(x, y, z, t) (  2*pi    * cos(2*pi*(y - t)) + 0*x + 0*z                    );
                obj.data.uZ_a  = @(x, y, z, t) (  2*pi    * cos(2*pi*(z - t)) + 0*x + 0*y                    );
                obj.data.uT_a  = @(x, y, z, t) (- 2*pi    * cos(2*pi*(x - t)) - 2*pi* cos(2*pi*(y - t)) - 2*pi* cos(2*pi*(z - t)) );
                
            elseif strcmp(obj.type.data, 'constant')
                
                obj.data.u_a   = @(x, y, z, t) (0 + sin(x- t) + sin(y - t) + 0*z);
                obj.data.uT_a  = @(x, y, z, t) (- cos(x - t) - cos(y - t) + 0*x + 0*y + 0*z);
                obj.data.uX_a  = @(x, y, z, t) (cos(x - t) + 0*x + 0*y + 0*z);
                obj.data.uY_a  = @(x, y, z, t) (cos(y - t) + 0*x + 0*y + 0*z);
                obj.data.uZ_a  = @(x, y, z, t) (0 + 0*x + 0*y + 0*z);
                
            end
        end
        function obj = createMatrices(obj)
            
            obj.disc.A  = obj.cont.a;
            obj.disc.B  = obj.cont.b;
            obj.disc.C  = obj.cont.c;
                                
        end
        function obj = createData(obj)
            
            obj.data.gE   = @(x, y, z, t) (obj.data.u_a (x, y, z, t));
            obj.data.gW   = @(x, y, z, t) (obj.data.u_a (x, y, z, t));
            obj.data.gN   = @(x, y, z, t) (obj.data.u_a (x, y, z, t));
            obj.data.gS   = @(x, y, z, t) (obj.data.u_a (x, y, z, t));
            obj.data.gF   = @(x, y, z, t) (obj.data.u_a (x, y, z, t));
            obj.data.gB   = @(x, y, z, t) (obj.data.u_a (x, y, z, t));
            obj.data.f    = @(x, y, z, t) (obj.data.u_a (x, y, z, t));

        end
    end
end
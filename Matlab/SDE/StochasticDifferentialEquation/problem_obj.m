classdef problem_obj
    % problem object
    % ut + Aux + Buy = e(Cuxx + Duyy)
    %       Hu = g
    %        u = f
    
    properties
        % Parameters
        cont;    % Structure with continuous matrices: a, ax, b, by, c, d

        % Problem type
        type;   % Type: problem type, data type, coeff type
        
        % Data
        data;  % Boundary data: gE, gED, gW, gWD, ...        
        
    end
    methods
        function obj = initialize(obj)
            
            obj = initializeProblem(obj);
            obj = initializeData(obj);
            
        end
        function obj = initializeProblem(obj)
            
            % dX = a dt + b dWt
            
            if strcmp(obj.type.problem, 'SDE')
                
                obj.cont.a   = @(x, t)(2*x.*(1 - x) + 0*t);
                obj.cont.b   = @(x, t)(0.25*x + 0*t);
                obj.cont.b_x = @(x, t)(0.25.^2*x + 0*t);
                 
            elseif strcmp(obj.type.problem, 'Black-Scholes')
                
                obj.cont.a   = @(x, t)(2*x + 0*t);
                obj.cont.b   = @(x, t)(1*x + 0*t);
                obj.cont.b_x = @(x, t)(1^2*x + 0*t);
                
            elseif strcmp(obj.type.problem, 'Black-Scholes2')
                
                obj.cont.a   = @(x, t)(2*x + 0*t);
                obj.cont.b   = @(x, t)(0.1*x + 0*t);
                obj.cont.b_x = @(x, t)(0.1^2*x + 0*t);
                
            elseif strcmp(obj.type.problem, 'Langevin')
                
                obj.cont.a  = 0;
                obj.cont.aW = -10;
                obj.cont.b  = 1;
                obj.cont.bW = 0;
               
            elseif strcmp(obj.type.problem, 'Zero')
                
                obj.cont.a   = @(x, t)(1 + 0*x + 0*t);
                obj.cont.b   = @(x, t)(10 + 0*x + 0*t);
                obj.cont.b_x = @(x, t)(0*x + 0*t);
                
            end
        end
        function obj = initializeData(obj)
            
            if strcmp(obj.type.data, 'constant')
                
                obj.data.u_a = @(t) (1 + 0*t);
            
            elseif strcmp(obj.type.data, 'SDE')
                
                obj.data.u_a = @(t) (0.5 + 0*t);
            
            elseif strcmp(obj.type.data, 'zero')
                
                obj.data.u_a = @(t) (0 + 0*t);
            
            end
        end
        function obj = createData(obj)
                
            obj.data.f   = @(t) (obj.data.u_a(t));

        end
    end
end
% clear all
% Nt = 2;
% Nx = 2;
% h = 1/Nx;
% k = 1/Nt;
% a = -1;
% 
% block1 = eye(Nt+1);
% block2 = eye(Nt+1) + diag(ones(Nt,1)*a*k, 1)/(h-a*k);
% block3 = diag(-ones(Nt+1,1)*h, 0)/(h-a*k);
% block3(end,end) = 0;
% 
% B = eye(Nx+1);
% B(1,1) = 0;
% B = kron(B, block2);
% C = diag(ones(Nx,1), -1);
% C = kron(C, block3);
% D = zeros(Nx+1);
% D(1,1) = 1;
% D = kron(D, block1);
% 
% A = B+C+D;
% 
% f = @(x) (sin(2*pi*x));
% F = zeros(Nt+1,1);
% F(1,1) = 1;
% x = (0:h:1)';
% F = kron(F, f(x));
% 
% g = @(t) (sin(2*pi*(1 - a*t)));
% G = zeros(Nx+1,1);
% G(end,end) = 1;
% t = (0:k:1)';
% G = kron(g(t), G);
% G(Nx+1) = 0;
% 
% H = F + G;
% 
% A = sparse(A);
% H = sparse(H);
% 
% u = A\H;
% %u = gmres(A,H);
% u = reshape(u, Nx+1, Nt+1);
% [X, T] = meshgrid(x,t);
% 
% figure
% surf(X, T, u)
% 
% ua = @(x, t) (sin(2*pi*(x - a*t)));
% figure
% surf(X, T, ua(X,T)')
% figure
% surf(X, T, ua(X,T)' - u)

%%

clear all
Nt = 120;
Nx = 120;
h = 1/Nx;
k = 1/Nt;
a = -1;

block1 = eye(Nt+1);
block2 = eye(Nt+1) + diag(ones(Nt,1)*a*k, 1)/(2*h) - diag(ones(Nt,1)*a*k, -1)/(2*h);
block3 = -diag(ones(Nt+1,1), 0);

if a < 0
    block3(end, end) = 0;
    block2(end, end - 1) = 0;
    block2(1, 2) = a*k/h;
    block2(1, 1) = 1 - a*k/h;
elseif a > 0
    block3(1,1) = 0;
end

B = eye(Nx+1);
B(1,1) = 0;
B = kron(B, block2);
C = diag(ones(Nx,1), -1);
C = kron(C, block3);
D = zeros(Nx+1);
D(1,1) = 1;
D = kron(D, block1);

A = B+C+D;

f = @(x) (sin(2*pi*x));
F = zeros(Nt+1,1);
F(1,1) = 1;
x = (0:h:1)';
F = kron(F, f(x));

g = @(t) (sin(2*pi*(1 - a*t)));
G = zeros(Nx+1,1);
G(end,end) = 1;
t = (0:k:1)';
G = kron(g(t), G);
G(Nx+1) = 0;

H = F + G;

A = sparse(A);
H = sparse(H);

u = A\H;
%u = gmres(A,H);
u = reshape(u, Nx+1, Nt+1);
[X, T] = meshgrid(x,t);

figure
surf(X, T, u)

ua = @(x, t) (sin(2*pi*(x - a*t)));
figure
surf(X, T, ua(X,T)')
figure
surf(X, T, ua(X,T)' - u)

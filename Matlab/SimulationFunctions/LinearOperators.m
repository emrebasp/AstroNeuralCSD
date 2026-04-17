function [x, kappa, Dxx] = LinearOperators(nx, Lx, check)

% [X, kappa, DXX] = LINEAROPERATORS(NX, LX) provides the discretized points X, 
% connectivity kernel WCONNECTIVITY either in spatial domain or in Fourier domain, 
% and the second-order finite differences Dxx. Periodic boundary conditions 
% are used for Dxx.
 
% LINEAROPERATORS(NX, LX, CHECK): NX is the number of discretization nodes on the
% cortical layer. LX is the length of the cortical layer. CHECK is
% either 'spatial' or 'fourier'. This determines if wConnectivity will be
% provided in Fourier domain (CHECK='fourier') or in spatial domain (CHECK='spatial').

% Grid parameters
hx = Lx/(nx-1);             % distance between two discretization nodes 
x = -Lx/2 + [0:nx-1]'*hx;   % discretization nodes

% Connectivity kernel function
kappaFun = @(x) 0.5*exp(-abs(x));

if check == 'spatial'
    % Connectivity kernel in spatial domain
    kappa = kappaFun(x);
elseif check == 'fourier'
    % Connectivity kernel in Fourier domain
    kappa = fft(kappaFun(x));
end

% Ancillary vector of ones for differentiation matrix
e = ones(nx,1); 

% Differentiation matrix in x: second-order central finite differences with periodic BCs
Dxx = spdiags([e -2*e e], [-1 0 1], nx, nx); 
Dxx(1,nx) = 1; 
Dxx(nx,1) = 1; 
Dxx = Dxx/hx^2;

end


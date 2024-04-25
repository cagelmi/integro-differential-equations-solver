function nominales = idsolver(lim,n,IC,c,d,k,alpha,beta,Tol,Flag)
% IDSOLVER solves a wide variety of integro-differential equations (IDE) of
% arbitrary order, including the Volterra and Fredholm IDE, variable limits
% on the integral, and non-linear IDE.
%
% Last update: 01/19/2013.
% Copyright (c) 2012 Claudio Gelmi & Héctor Jorquera.

home
disp(' ')
disp('***********************************************************************')
disp('**   IDSOLVER: A general purpose solver for nth-order ID equations   **')
disp('**   Copyright (2012) <a href = "http://www.iiq.cl" >Dr. Claudio Gelmi</a> and <a href = "http://www.ing.puc.cl/jorquera" >Dr. Héctor Jorquera</a>      **')
disp('***********************************************************************')
disp(' ')

% Tolerance options
if Flag == 0
    options = [];
    TolQuad = [];
elseif Flag == 1
    options = odeset('AbsTol',1e-8,'RelTol',1e-8);
    TolQuad = 1e-8;
else
    disp(' ')
    disp('Flag = 0 or Flag = 1. No other options are supported.')
    return
end

% Interval partition
M = 100;
interval = linspace(lim(1),lim(2),M);
% Initial guess generator
warning('off')
[x nominales] = ode113(@model0,interval,IC,options,c,n);
nominales = [x nominales(:,1)];

% Iterative solution
error = 1e3;
iteration = 1;
fprintf('  Error convergence\n ');
fprintf(' ================= \n');
fprintf(' Iteration    Error \n');
while error > Tol
    [x y] = ode113(@model,interval,IC,options,nominales,c,n,k,alpha,beta,d,TolQuad);
    error = sum((y(:,1)-nominales(:,2)).^2);
    fprintf('  %4i      %8.2e\n',[iteration error]);
    a = 0.5;
    nominales = [x (1-a)*nominales(:,2)+a*y(:,1)];
    iteration = iteration+1;
end
warning('on')
disp(' ')

% Final message
if Flag == 0
    disp('Note: Problem solved using MATLAB''s default tolerances.')
elseif Flag == 1
    disp('Note: Problem solved using RelTol = AbsTol = 1e-8.')
end

% Saves and plot the final answer. The first column of ansidsolver shows the
% independent variable x, and the second one shows y(x).
save ansidsolver.txt nominales -ascii
plot(nominales(:,1),nominales(:,2))
xlabel('X'), ylabel('Y')

function dy = model0(x,y,c,n)
% Initial guess generator
dy = zeros(n,1);
if n > 1
    dy(1:n-1) = y(2:n);
end
dy(n) = c(x,y(1));

function dy = model(x,y,nominales,c,n,k,alpha,beta,d,TolQuad)
% Interpolation step
ys = @(s) interp1(nominales(:,1),nominales(:,2),s);
% Integro-differential equation
dy = zeros(n,1);
if n > 1
    dy(1:n-1) = y(2:n);
end
dy(n) = c(x,y(1)) + d(x)*quadl(@(s) k(x,s).*ys(s),alpha(x),beta(x),TolQuad);
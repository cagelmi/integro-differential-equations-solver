function examples
% EXAMPLES illustrates the use of IDSOLVER. 
% 
% The method and nomenclature are described in the following manuscript:
%
% Gelmi C. and H. Jorquera. IDSOLVER: A general purpose solver for nth-order
% integro-differential equations.
% 
% Last update: 01/19/2013.

%% Example 1.

xinterval = [0 1];          % Interval of the independent variable x
n = 1;                      % Order of the IDE
InitCond = 0;               % Initial condition
c = @(x,y) y-0.5*x+1/(1+x)-log(1+x); % See Eq. (1) in manuscript
d = @(x) 1/log(2)^2;                 % See Eq. (1) in manuscript
k = @(x,s) x./(1+s);        % Kernel
alpha = @(x) 0;             % Lower limit of the integral
beta = @(x) 1;              % Upper limit of the integral
Tol = 1e-8;                 % Tolerance
Flag = 0;                   % Setting MATLAB’s default tolerances

% IDSOLVER call:
idsolver(xinterval,n,InitCond,c,d,k,alpha,beta,Tol,Flag)

%% Example 2.

xinterval = [0 1];          % Interval of the independent variable x
n = 1;                      % Order of the IDE
InitCond = 1;               % Initial condition
c = @(x,y) y-cos(2*pi*x)-2*pi*sin(2*pi*x)-0.5*sin(4*pi*x);  % See Eq. (1) in manuscript
d = @(x) 1;                     % See Eq. (1) in manuscript
k = @(x,s) sin(4*pi*x+2*pi*s);  % Kernel
alpha = @(x) 0;             % Lower limit of the integral
beta = @(x) 1;              % Upper limit of the integral
Tol = 1e-8;                 % Tolerance
Flag = 0;                   % Setting MATLAB’s default tolerances

% IDSOLVER call:
idsolver(xinterval,n,InitCond,c,d,k,alpha,beta,Tol,Flag)

%% Example 3.

xinterval = [0 1];          % Interval of the independent variable x
n = 1;                      % Order of the IDE
InitCond = 1;               % Initial condition
c = @(x,y) 1-29*x/60;       % See Eq. (1) in manuscript
d = @(x) 1;                 % See Eq. (1) in manuscript
k = @(x,s) x.*s;            % Kernel
alpha = @(x) 0;             % Lower limit of the integral
beta = @(x) 1;              % Upper limit of the integral
Tol = 1e-8;                 % Tolerance
Flag = 0;                   % Setting MATLAB’s default tolerances

% IDSOLVER call:
idsolverEx3(xinterval,n,InitCond,c,d,k,alpha,beta,Tol,Flag)

%% Example 4.

xinterval = [0 1];          % Interval of the independent variable x
n = 1;                      % Order of the IDE
InitCond = 1;
c = @(x,y) x*(1+x^0.5)*exp(-x^0.5)-(x^2+x+1)*exp(-x);  % See Eq. (1) in manuscript
d = @(x) 1;                 % See Eq. (1) in manuscript
k = @(x,s) x.*s;            % Kernel
alpha = @(x) x;             % Lower limit of the integral
beta = @(x) x^0.5;          % Upper limit of the integral
Tol = 1e-8;                 % Tolerance
Flag = 0;                   % Setting MATLAB’s default tolerances

% IDSOLVER call:
idsolver(xinterval,n,InitCond,c,d,k,alpha,beta,Tol,Flag)

%% Example 5.

xinterval = [0 1];          % Interval of the independent variable x
n = 4;                      % Order of the IDE
InitCond = [1 1 1 1];       % Initial condition
c = @(x,y) exp(x)-x;        % See Eq. (1) in manuscript
d = @(x) 1;                 % See Eq. (1) in manuscript
k = @(x,s) x.*s;            % Kernel
alpha = @(x) 0;             % Lower limit of the integral
beta = @(x) 1;              % Upper limit of the integral
Tol = 1e-8;                 % Tolerance
Flag = 0;                   % Setting MATLAB’s default tolerances

% IDSOLVER call:
idsolver(xinterval,n,InitCond,c,d,k,alpha,beta,Tol,Flag)


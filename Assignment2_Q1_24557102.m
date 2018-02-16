clc
clear
clear all
%% Initialisation of parameters
m = 64; % Number of grid points
init_point = 0; % Initial point
end_point = 1; % End point
range = end_point-init_point; % Interval Length

%% Calculation of h and x,y, depending on the interval
h = range/(m+1); % Grid Interval
x = init_point:h:end_point; % Generating Nx grid point.
y = init_point:h:end_point; % Generating Ny grid point
[xnew,ynew] = meshgrid(x,y);

%% Definition of Boundary Condition
gb = 0*x + 0*y; % u(x,0)
gt = 0*x + 0*y; % u(x,end_point)
gl = 0*x + 0*y; % u(0,y)
gr = 0*x + 0*y; % u(init_point,y)

%% Definition of Equation (A*T = b) Components
f = -2*xnew; % Generating f in Poisson equation
A = Form_A(m); % Generating A (LHS of the problem)
b = Form_B(m,h,f,gb,gt,gl,gr); % Generating b (RHS of the problem)

%% Solve for T
% Performing Fourier Transform upon b
b_T = FTransform(m,b)';
[b_T_mod] = mat_mod(m,b_T);

% Solve for T, T = A^-1 * b
T = A\b_T_mod;
[T_mod] = mat_mod(m,T);

% Performing Fourier Transform upon T
U = FTransform(m,T_mod)';

% Rearranging U to be a square matrix for plotting purpose
U_mat1 = zeros(m+2);
for j = 2:m+1
    U_mat1(2:m+1,j) = U((j-2)*m + 1:(j-2)*m + m);
end

%% Graphical Plotting
figure(1)
mesh(xnew,ynew,U_mat1);
title('U_X_X + U_Y_Y = -2X');
xlabel('x'); ylabel('y'); zlabel('U^h')








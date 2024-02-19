% Project 1
clc; clear;

syms x y mu;

P1 = sqrt((x+mu)^2 + y^2);
P2 = sqrt((x-1+mu)^2 + y^2);

U = 0.5*(x^2 + y^2) + (1-mu)/P1 + mu/P2;

% Partial derivatives
U_x = diff(U, x);
U_xx = diff(U, x, 2);
U_xy = diff(U_x, y);
U_yy = diff(U, y, 2);

A = [0 0 1 0;
    0 0 0 1;
    U_xx U_xy 0 2;
    U_xy U_yy -2 0];

% Lists the values of mu for sun-earth, earth-moon, and saturn-titan respectively
Mu = [3.0039*(10^(-7)); 1.2151*(10^(-2)); 2.366*(10^(-4))];

% X and Y value pairs of the lagrange points for various CR3BP systems
L4_Points = [0.5*Mu(1), sqrt(3)/2; % Sun-Earth system (system 1)
             0.5*Mu(2), sqrt(3)/2; % Earth-Moon system (system 2)
             0.5*Mu(3), sqrt(3)/2;]; % Saturn-Titan System (system 3)

L5_Points = [0.5*Mu(1), -sqrt(3)/2; % Sun-Earth system 
             0.5*Mu(2), -sqrt(3)/2; % Earth-Moon system
             0.5*Mu(3), -sqrt(3)/2;]; % Saturn-Titan System

%% Calculating the eigenvalues
% Step 1: Plug x & y values into matrix A
% Step 2: Compute eigenvalues

% L4 Points
fprintf("L4 Points: \n");
for i = 1: 3
    A_subbed = subs(A, [x y mu], [L4_Points(i,:) Mu(i)] );
    eigenvalues = eig(A_subbed);
    fprintf("System %d: %f\t%f\t%f\t%f\n\n", i, eigenvalues);
end

% L5 Points
fprintf("L5 Points: \n");
for i = 1: 3
    A_subbed = subs(A, [x y mu], [L5_Points(i,:) Mu(i)] );
    eigenvalues = eig(A_subbed);
    fprintf("System %d: %f\t%f\t%f\t%f\n\n", i, eigenvalues);
end
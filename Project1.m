% Project 1
% Group Members: Nicholas Luis, Shawn Watkins, Matthew Ominsky
clc; clear;

syms x y mu;

P1 = sqrt((x+mu)^2 + y^2);
P2 = sqrt((x-1+mu)^2 + y^2);

U = 0.5*(x^2 + y^2) + ((1-mu)/P1) + (mu/P2);

% Partial derivatives
U_x = diff(U, x);
U_xx = diff(U, x, 2);
U_xy = diff(U_x, y);
U_yy = diff(U, y, 2);

A = [0 0 1 0;
    0 0 0 1;
    U_xx U_xy 0 2;
    U_xy U_yy -2 0;];

% Lists the values of mu for sun-earth, earth-moon, and saturn-titan respectively
Mu = [3.0039*(10^(-7)); 1.2151*(10^(-2)); 2.366*(10^(-4))];

% X and Y value pairs of the lagrange points for various CR3BP systems
L1_Points = [0.995363, 0;
             0.836915, 0;
             0.9575,   0;];
L2_Points = [1.004637, 0;
             1.15568,  0;
             1.0425,   0;];
L3_Points = [-1.00001, 0;
             -1.00506, 0;
             -1.0001,  0;];
L4_Points = [0.5-Mu(1), sqrt(3)/2; % Sun-Earth system (system 1)
             0.5-Mu(2), sqrt(3)/2; % Earth-Moon system (system 2)
             0.5-Mu(3), sqrt(3)/2;]; % Saturn-Titan System (system 3)

L5_Points = [0.5-Mu(1), -sqrt(3)/2; % Sun-Earth system 
             0.5-Mu(2), -sqrt(3)/2; % Earth-Moon system
             0.5-Mu(3), -sqrt(3)/2;]; % Saturn-Titan System

%% Calculating the eigenvalues
% Step 1: Plug x & y values into matrix A
% Step 2: Compute eigenvalues

% L1 Points
fprintf("L1 Points: \n");
for i = 1: 3
    A_subbed = subs(A, [x y mu], [L1_Points(i,1),L1_Points(i,2), Mu(i)] );
    eigenvalues = eig(A_subbed);
    eigenvaluePairs = [real(eigenvalues), imag(eigenvalues)]; % Separates real and imaginary
    fprintf('System %d: %f%+fi\t%f%+fi\t%f%+fi\t%f%+fi\n\n', i, eigenvaluePairs');
end

% L2 Points
fprintf("L2 Points: \n");
for i = 1: 3
    A_subbed = subs(A, [x y mu], [L2_Points(i,1),L2_Points(i,2), Mu(i)] );
    eigenvalues = eig(A_subbed); 
    eigenvaluePairs = [real(eigenvalues), imag(eigenvalues)]; % Separates real and imaginary
    fprintf('System %d: %f%+fi\t%f%+fi\t%f%+fi\t%f%+fi\n\n', i, eigenvaluePairs');
end

% L3 Points
fprintf("L3 Points: \n");
for i = 1: 3
    A_subbed = subs(A, [x y mu], [L3_Points(i,1),L3_Points(i,2), Mu(i)] );
    eigenvalues = eig(A_subbed);
    eigenvaluePairs = [real(eigenvalues), imag(eigenvalues)]; % Separates real and imaginary
    fprintf('System %d: %f%+fi\t%f%+fi\t%f%+fi\t%f%+fi\n\n', i, eigenvaluePairs');
end

% L4 Points
fprintf("L4 Points: \n");
for i = 1: 3
    A_subbed = subs(A, [x y mu], [L4_Points(i,1),L4_Points(i,2), Mu(i)] );
    eigenvalues = eig(A_subbed);
    eigenvaluePairs = [real(eigenvalues), imag(eigenvalues)]; % Separates real and imaginary
    fprintf('System %d: %f%+fi\t%f%+fi\t%f%+fi\t%f%+fi\n\n', i, eigenvaluePairs');
end

% L5 Points
fprintf("L5 Points: \n");
for i = 1: 3
    A_subbed = subs(A, [x y mu], [L5_Points(i,1),L5_Points(i,2), Mu(i)] );
    eigenvalues = eig(A_subbed);
    eigenvaluePairs = [real(eigenvalues), imag(eigenvalues)];
    fprintf('System %d: %f%+fi\t%f%+fi\t%f%+fi\t%f%+fi\n\n', i, eigenvaluePairs');
end
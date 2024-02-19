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

% Plugging into A matrix
A = [0 0 1 0;
    0 0 0 1;
    U_xx U_xy 0 2;
    U_xy U_yy -2 0];

% 
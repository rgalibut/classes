clear; clc;

%% Map each pole and zero into z-plane using z = exp(s*h)
sym s;              % set up symbolic 's'
h = 0.01;           % timestep 1/100

% Set up y(s) and x(s) vectors
% Test case using 2a/2b
a = RR_poly([-1 1 -3 3 -6 6],1); 
b = RR_poly([-2 2 -5 5],1); 
%b = RR_poly(0,1); 
f = RR_poly([1 1 3 3 6 6 20 20 20 20 20 20],1);
[xs,ys] = RR_diophantine(a,b,f);

% Create D(s)
Ds = RR_tf(ys,xs);

% Call function
Dz = RJG_C2D_matched(Ds,h)

% Check with MATLAB c2d function
fprintf('\nCheck with MATLAB c2d function: \n');
c2d(tf(b.poly,a.poly),0.01,'matched')

%% FUNCTIONS
function [Dz] = RJG_C2D_matched(Ds,h)
% function [Dz] = RJG_C2D_matched(Ds,h)
% Convert D(s) to D(z) using the matched z-transform method with timestep
% h.
    % Part (i): Map each pole and zero into the z-plane using z = exp(s*h)
    % Zeroes (ys)
    for i = 1:length(RR_roots(Ds.num))
        z_zeros = exp(RR_roots(Ds.num)*h);
    end % for loop
    % Poles (xs)
    for i = 1:length(RR_roots(Ds.den))
        z_poles = exp(RR_roots(Ds.den)*h);
    end % for loop

    % Part (ii): Map infinite zeros to z = -1 in D(z)
    for i = 1:length(RR_roots(Ds.num))
        if RR_roots(Ds.num) == 0
            z_zeros(i) = -1;
        else
            z_zeros(i) = z_zeros(i);
        end % if statement
    end % for loop

    % Part (iii): Adjust gain of D(z) by setting s=0 and z=1 and solve for
    % gain K
    % Zeros (ys)
    % for i = 1:length(z_zeros)
    % 
    % end
    % % Poles (xs)
    % for j = 1:length(z_poles)
    % 
    % end

    % Create Dz
    Dz = RR_tf(z_zeros,z_poles);

end % function RJG_C2D_matched
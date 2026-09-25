function [X] = NSPD(A)
% NSPD computes the nearest symmetric positive definite matrix to 
% the input. This does not assume that the input is symmetric, even if it
% commonly will be.
% Based on the blog post:
% https://nhigham.com/2021/01/26/what-is-the-nearest-positive-semidefinite-matrix/
% And the algorithm described in the paper by Higham (1988):
% https://doi.org/10.1016/0024-3795(88)90223-6
% The algorithm has been slightly modified from the original by setting
% the minimum possible value for eigenvalues to a small number (1e-6)
% ensuring positive definiteness not just "semi" positive definiteness.

B=(A+A')./2;
[Q,d] =eig(B,'vector'); 
X=Q*(max(d,1e-6).*Q'); 
X=(X+X')/2;



end


function [prox_sol,num]=prox_l1(x,lam,t,N)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the proximal mapping of g(x)=\norm{x}_{1}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
%
% x ....................... the given vector
% t ....................... the parameter of the prox mapping
%
% OUTPUT:
%
% prox_sol.................  the solution of the prox computation

lambda=lam/t;

prox_sol=sign(x).*max(abs(x)-lambda,0);
num=1;

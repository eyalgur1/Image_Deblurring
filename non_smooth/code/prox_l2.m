function [prox_sol,num]=prox_l2(x,lam,t,N)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the proximal mapping of g(x)=\norm{x}_{2}
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

c=norm(x);
prox_sol=x*max((c-lambda)/c,0);
num=1;

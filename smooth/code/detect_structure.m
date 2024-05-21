function [p,S,d]=detect_structure(A)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function detects the structure of a given matrix A. 
%
% INPUT
%
% A .............................. A matrix
%
% OUPUT
% 
% p .............................. number of structure components
% S .............................. three dimensional array containing the p
%                                    structure matrices. 
%
% EXAMPLE
% Suppose that the input matrix is
% A=[0.234 0.234
%       0.14    0.234
%       0         0.14]
% After running the commans
% [p,S]=detect_structure(A)
% we obtain p=2 and 
% S(:,:,1)= [0.234 0.234
%                     0    0.234
%                     0        0]
% S(:,:,2)=[    0     0
%                  0.14   0
%                   0      0.14]


d=unique(A(find(A)));
p=length(d);
[m,n]=size(A);
P=zeros(m,n,p);
for i=1:p
    S(:,:,i)=d(i)*(A==d(i));
end

% =========================================================================
% 					Bravyi-Kitaev Superfast Simulation
% Author: 	Kanav Setia*
% Date:     2016/11/08
% Version: 	0.1
%
% *Whitfield Group, Department of Physics and Astronomy,
%  Dartmouth College
% =========================================================================
% HELP SECTION
% bj_general is the function to calculate B_k for the hamiltonian in bravyi-kitaev 
% algorithm(the one with just the nearest neighbor hopping terms in the 
% Hamiltonian).
% NOTE: This is not the general algorithm since the hamiltonian does not
% contain all the hopping terms specified by the G matrix
% 
% Input: totalqubits, k (in B_k)
% Output: B_k
% 
% A similar function bj_genstring generates the algebraic string for Bj
% See also Ajk_general, Ajk_genstring, bj_genstring
% =========================================================================

function y=bj_general(edges,j)
x=[0  1 ; 1   0];
y=[0 -1j; 1j  0];
z=[1  0 ; 0  -1];
I=eye(2);
A=eye;
S=[];
for d1=1:length(edges)
    if edges(d1,1)==j|edges(d1,2)==j
        A=kron(z,A);
        edstr= strcat('Z',strcat(num2str(edges(d1,1)),num2str(edges(d1,2))));
        edstr=strcat(strcat('(',edstr),')');
        S=strcat(edstr,S);
    else        
        A=kron(I,A);
        edstr=strcat('I',strcat(num2str(edges(d1,1)),num2str(edges(d1,2))));
        edstr=strcat(strcat('(',edstr),')');
        S=strcat(edstr,S);
    end
end

y=A;
S
end

        
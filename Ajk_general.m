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
% Ajk_general is the function to calculate A_jk for the hamiltonian in bravyi-
% kitaev algorithm(the one with just the nearest neighbor hopping terms 
% in the Hamiltonian).
% Input: edges (matrix specifying edges. this matrix contains as its rows,
%        j,k), en (row number of edges to get jk for Ajk)
% Output: A_jk
% 
% A similar function Ajk_genstring generates the algebraic string for Ajk.
% 
% Convention for qubits --> qubits in 1st row--> qubits in 2nd row--> ...
% 
% See also, bj_general, bj_genstring, Ajk_genstring, bk_general
% =========================================================================
function y=Ajk_general(edges,en)
A=eye;
x=[0  1 ; 1   0];
y=[0 -1j; 1j  0];
z=[1  0 ; 0  -1];
I=eye(2);
j=edges(en,1);
k=edges(en,2);
edstr=[];
S=[];
for d1=1:length(edges)
    if edges(d1,1)==j & edges(d1,2)==k
        A=kron(x,A);
        edstr= strcat('X',strcat(num2str(edges(d1,1)),num2str(edges(d1,2))));
        edstr=strcat(strcat('(',edstr),')');
        S=strcat(edstr,S);
    elseif edges(d1,1)==j
        if edges(d1,2)<k
            A=kron(z,A);
            edstr= strcat('Z',strcat(num2str(edges(d1,1)),num2str(edges(d1,2))));
            edstr=strcat(strcat('(',edstr),')');
            S=strcat(edstr,S);
        else
            A=kron(I,A);
            edstr= strcat('I',strcat(num2str(edges(d1,1)),num2str(edges(d1,2))));
            edstr=strcat(strcat('(',edstr),')');
            S=strcat(edstr,S);
        end
    elseif  edges(d1,2)==j
        if edges(d1,1)<k
            A=kron(z,A);
            edstr= strcat('Z',strcat(num2str(edges(d1,1)),num2str(edges(d1,2))));
            edstr=strcat(strcat('(',edstr),')');
            S=strcat(edstr,S);
        else
            A=kron(I,A);
            edstr= strcat('I',strcat(num2str(edges(d1,1)),num2str(edges(d1,2))));
            edstr=strcat(strcat('(',edstr),')');
            S=strcat(edstr,S);
        end
    elseif edges(d1,1)==k
        if edges(d1,2)<j
            A=kron(z,A);
            edstr= strcat('Z',strcat(num2str(edges(d1,1)),num2str(edges(d1,2))));
            edstr=strcat(strcat('(',edstr),')');
            S=strcat(edstr,S);
        else
            A=kron(I,A);
            edstr= strcat('I',strcat(num2str(edges(d1,1)),num2str(edges(d1,2))));
            edstr=strcat(strcat('(',edstr),')');
            S=strcat(edstr,S);
        end
    elseif edges(d1,2)==k
        if edges(d1,1)<j
            A=kron(z,A);
            edstr= strcat('Z',strcat(num2str(edges(d1,1)),num2str(edges(d1,2))));
            edstr=strcat(strcat('(',edstr),')');
            S=strcat(edstr,S);
        else
            A=kron(I,A);
            edstr= strcat('I',strcat(num2str(edges(d1,1)),num2str(edges(d1,2))));
            edstr=strcat(strcat('(',edstr),')');
            S=strcat(edstr,S);
        end
    else
        A=kron(I,A);
        edstr= strcat('I',strcat(num2str(edges(d1,1)),num2str(edges(d1,2))));
        edstr=strcat(strcat('(',edstr),')');
        S=strcat(edstr,S);
    end    
end
 y=A;
 S
end

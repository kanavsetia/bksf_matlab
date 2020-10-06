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
% bj is the function to calculate B_k for the hamiltonian in bravyi-kitaev 
% algorithm(the one with just the nearest neighbor hopping terms in the 
% Hamiltonian).
% NOTE: This is not the general algorithm since the hamiltonian does not
% contain all the hopping terms specified by the G matrix
% 
% Input: totalqubits, k (in B_k)
% Output: B_k
% 
% See also Ajk2
% =========================================================================

function y=bj(totqub,j)
x=[0  1 ; 1   0];
y=[0 -1j; 1j  0];
z=[1  0 ; 0  -1];
A=eye;
if j==1
    A=kron(z,A);
    for s1=2:totqub-1
        A=kron(eye(2),A);
    end
    A=kron(z,A);
elseif j==2    
    A=kron(z,A);
    A=kron(z,A);
    for s2=3:totqub
        A=kron(eye(2),A);
    end
else
    for s3=1:j-2
        A=kron(eye(2),A);
    end
    A=kron(z,A);
    A=kron(z,A);
    for s4=j+1:totqub
        A=kron(eye(2),A);
    end
end
y=A;
end

        
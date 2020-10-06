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
% Ajk2 is the function to calculate A_jk for the hamiltonian in bravyi-
% kitaev algorithm(the one with just the nearest neighbor hopping terms 
% in the Hamiltonian).
% NOTE: This is not the general algorithm since the hamiltonian does not
% contain all the hopping terms specified by the G matrix
% 
% Input: totalqubits, j (in A_jk) (only j is required because k=j+1)
% Output: A_jk
% 
% See also bj
% =========================================================================
function y=Ajk2(totqub,j)
A=eye;
x=[0  1 ; 1   0];
y=[0 -1j; 1j  0];
z=[1  0 ; 0  -1];
if j==1
    A=kron(x,A);
    for i=2:totqub
        A=kron(eye(2),A);
    end
elseif j==2 && j==totqub
    A=kron(z,A);
    A=kron(x,A);    
elseif j==2
    A=kron(z,A);
    A=kron(x,A);
    if totqub==3
        A=kron(z,A);
    else
        for t=3:totqub
            A=kron(eye(2),A);
        end
    end    
elseif j==totqub
    A=kron(z,A);
    for s=2:totqub-1
        A=kron(eye(2),A);
    end
    A=kron(x,A);
elseif j==totqub-1    
    for k=1:j-2
        A=kron(eye(2),A);
    end
    A=kron(z,A);
    A=kron(x,A);
    A=kron(z,A);    
else
    for s1=1:j-2
        A=kron(eye(2),A);
    end
    A=kron(z,A);
    A=kron(x,A);
    for s2=j+1:totqub
        A=kron(eye(2),A);
    end    
end
y=A;
end


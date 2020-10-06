% =========================================================================
% Calculates the Annhilation Operator for the Jordan Wigner Simulation
% inputs: 
% 	tot_input = Total Number of qubits
%   quj = # of the qubit/fermion for which the annhilation operator is
%   being calculated
% 	
% output: Annhilataion Operator
%==========================================================================
function y=annihilator(tot_qub,quj)
x=[0  1 ; 1   0];
y=[0 -1j; 1j  0];
z=[1  0 ; 0  -1];

cre_dum=eye;
q = (1/2).*(x+1j.*y);
if quj-1>0
    for i=1:quj-1
        cre_dum=kron(z,cre_dum);
%         cre_dum=kron(z,cre_dum);
    end
    cre_dum=kron(q,cre_dum);
    if tot_qub-quj>0
        for i=quj+1:tot_qub
            cre_dum=kron(eye(2),cre_dum);
        end
    end
else
    cre_dum=kron(q,cre_dum);
    if tot_qub>1
        for i=2 :tot_qub
            cre_dum=kron(eye(2),cre_dum);
        end
    end
end
y=cre_dum;
end
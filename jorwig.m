% =========================================================================
%               Jordan-Wigner Simulation for fermions
% 
% Author: 	Kanav Setia*
% Date:     2016/11/08
% Version: 	0.1
% 
% *Whitfield Group, Department of Physics and Astronomy,
%  Dartmouth College
% =========================================================================
% HELP SECTION
% jorwig is the function to run the jordan-wigner algorithm. This algorithm
% considers the general case and the hamiltonian is constructed from the
% hopping matrix, G.
% 
% Input: totalqubits, Hamiltonian (G matrix)
% Output: Data_JW (All the data for the algorithm in a data structure)
% Data Format for the Data structure
% 
% Operator(1).T --> Kinetic Operator
% Operator(1).V --> Potential Operator
% Operator(1).Ntotal--> Number Operator
% Operator(1).H --> Hamiltonian
% Operator(1).expectation.Time--> time vector
% Operator(1).expectation.tau--> Kinetic Energy Expectation Values
% corresponding to time vector
% Operator(1).expectation.nu--> Potential Energy Expectation Values
% corresponding to time vector
% Operator(1).expectation.probability--> Column vectors of this matrix
% contains probability values of respective particles corresponding to
% time vector
% Operator(1).expectation.n--> Number Expectation Values of different
% particles corresponding to time vector
% See also, annihilator
% =========================================================================

function f2=jorwig(totalqubits,Hamiltonian)

% =========================================================================
% Defining the Pauli Matrices
% =========================================================================

x=[0  1 ; 1   0];
y=[0 -1j; 1j  0];
z=[1  0 ; 0  -1];

% =========================================================================
% Defining the Structure for Operators. This Matlab data structure will
% hold all the Operators in one place.
% Operator(i).Cre holds the creation operator for the i'th fermion. It is a
% tensor which has the sigma z's included to get the right parity.
% Similarly Operator(i).Ann gives the annihilation operator.
% Operator(1).H is a place holder for the hamiltonian.
% 
% =========================================================================

Operator=[];

% =========================================================================
% Initialization Matrix and the number of qubits.
% =========================================================================


G=Hamiltonian;
tot_qub=length(G);

% =========================================================================
% Initialization of different Operators.
% =========================================================================
Operator(1).T=zeros(2^tot_qub);         %Kinetic Operator
Operator(1).V=zeros(2^tot_qub);         %Potential Operator
Operator(1).Ntotal=zeros(2^tot_qub);    %Number Operator
Operator(1).H=zeros(2^tot_qub);
Operator(1).expectation.Time=linspace(0,10,200);
Operator(1).expectation.tau=[];         %Kinetic Expectation Values (at different times)
Operator(1).expectation.nu=[];          %Potential Expectation Values (at different times)
Operator(1).expectation.probability=[];
Operator(1).expectation.n=[];           %Number Expectation Values (at different times)
% =========================================================================
% Calculation of Operators
% annihilator.m is the function file that calculates the Annhilation
% Operator.
% Note: Annhilation and creation operators are stored in i'th place in the
% structure, where i is the number of the fermion on which we are
% operating. But H, V, and T are stored in the 1st branch of the structure.
% =========================================================================
for d1=1:tot_qub
Operator(d1).Ann = annihilator(tot_qub,d1);
Operator(d1).Cre = Operator(d1).Ann';
end
for d2=1:length(Operator)
    for d3=1:length(Operator)
        if d2==d3
            Operator(1).V=Operator(1).V+G(d2,d3)*Operator(d2).Cre*...
                Operator(d3).Ann;
            Operator(d2).N=Operator(d2).Cre*Operator(d3).Ann;
            Operator(1).Ntotal=Operator(1).Ntotal+Operator(d2).Cre*...
                Operator(d3).Ann;
            
        else
            Operator(1).T=Operator(1).T+G(d2,d3)*Operator(d2).Cre*...
                Operator(d3).Ann;
        end
    end
end
Operator(1).H=Operator(1).T+Operator(1).V;
% =========================================================================
% Interaction between two particles.
% =========================================================================
Operator(1).W=-7.*Operator(1).Cre*Operator(1).Ann*Operator(3).Cre*Operator(3).Ann;
Operator(1).H=Operator(1).H+Operator(1).W;
% eig(Operator(1).H)


% =========================================================================
% Hydrogen Hamiltonian
% =========================================================================
% 
h11=-1.252477;
h22=h11;
h33=-.475934;
h44=h33;
h1221=.674493;
h2112=h1221;
h3443=.697397;
h4334=h3443;
h1331=.663472;
h1441=h1331; h2332=h1331; h2442=h1331;h3113=h1331;h4114=h1331;h3223=h1331;
h4224=h1331;
h1313=0.181287;
h2424=h1313; h3241=h1313; h3421=h1313; h1423=h1313;h1243=h1313;

H_1=h11*Operator(1).Cre*Operator(1).Ann +...
    h22*Operator(2).Cre*Operator(2).Ann +...
    h33*Operator(3).Cre*Operator(3).Ann +...
    h44*Operator(4).Cre*Operator(4).Ann;
H_2=   h1221*Operator(1).Cre*Operator(2).Cre*Operator(2).Ann*Operator(1).Ann+...
       h3443*Operator(3).Cre*Operator(4).Cre*Operator(4).Ann*Operator(3).Ann+...
       h1441*Operator(1).Cre*Operator(4).Cre*Operator(4).Ann*Operator(1).Ann+...
       h2332*Operator(2).Cre*Operator(3).Cre*Operator(3).Ann*Operator(2).Ann+...
(h1331-h1313)*Operator(1).Cre*Operator(3).Cre*Operator(3).Ann*Operator(1).Ann+...
(h2442-h2424)*Operator(2).Cre*Operator(4).Cre*Operator(4).Ann*Operator(2).Ann;

H_2=H_2+h1243*Operator(1).Cre*Operator(2).Cre*Operator(4).Ann*Operator(3).Ann+...
        h1243*Operator(3).Cre*Operator(4).Cre*Operator(2).Ann*Operator(1).Ann+...
        h1423*Operator(1).Cre*Operator(4).Cre*Operator(2).Ann*Operator(3).Ann+...
        h1423*Operator(3).Cre*Operator(2).Cre*Operator(4).Ann*Operator(1).Ann;
 Operator(1).H=H_1+H_2;
        



% =========================================================================
% Separating Eigenvalues into sectors
% 
% This section is not required but we used it while developing to  
% understand the problem better.
% =========================================================================
[H_evec,H_eval]=eig(Operator(1).H);
% [H_evec,H_R]=qr(H_evec);
% H_evec_rearranged=[];
% H_eval_rearranged=[];
% H_evec_particles=[];
% % for d7=1:length(H_eval)    
% %         H_evec_particles(d7)=H_evec(:,d7)'*Operator(1).Ntotal*H_evec(:,d7);           
% % end
% H_evec_particles=diag(H_evec'*Operator(1).Ntotal*H_evec);
% [Hval,Hindex]=sort(H_evec_particles);
% for d8=1:length(Hindex)
%     H_evec_rearranged(:,d8)=H_evec(:,Hindex(d8));
%     H_eval_rearranged(d8)=H_eval(Hindex(d8),Hindex(d8));
% end
% H_evec_particles=Hval;
% H_eval_rearranged;
% % H_evec_rearranged
% =========================================================================
% State Initialization
% =========================================================================
psi0=zeros(2^tot_qub,1);
psi0(1)=1;
% psi0=(Operator(2).Cre*Operator(3).Cre*psi0);%+exp(-2i)*Operator(4).Cre*psi0)/sqrt(2);
psi0=(Operator(2).Cre*Operator(3).Cre*psi0);
%psid=0;
% 12
% 13
% 14
% 24
% 23
% % % 34
% coef=[1;1;1;1;1;1];
% A=(coef(3)*Operator(1).Cre*Operator(4).Cre*psi0)+...
%     (coef(5)*Operator(2).Cre*Operator(3).Cre*psi0)+...
%     (coef(1)*Operator(1).Cre*Operator(2).Cre*psi0)+...
%     (coef(6)*Operator(3).Cre*Operator(4).Cre*psi0)+...
%     (coef(2)*Operator(1).Cre*Operator(3).Cre*psi0)+...
%     (coef(4)*Operator(2).Cre*Operator(4).Cre*psi0);
%psi0=psi0-psid/norm(psid);
% (Operator(1).Cre*Operator(2).Cre*Operator(3).Cre*Operator(4).Cre*psi0);
% keyboard
% psi0=A;
psi0=psi0/norm(psi0);
% tau1=psi0'*(Operator(1).T)*psi0
% nu1=psi0'*(Operator(1).V)*psi0
% psi0=(Operator(5).Cre*psi0);%+exp(-2i)*Operator(4).Cre*psi0)/sqrt(2);
% psi0=psi0/sqrt(norm(psi0));
% keyboard
% psi0=(Operator(1).Cre*Operator(2).Cre*psi0);%+exp(-2i)*Operator(4).Cre*psi0)/sqrt(2);
% psi0=psi0/sqrt(norm(psi0));
%norm(psi0)
% =========================================================================
% Time Evolution of the State
% =========================================================================
for d4=1:length(Operator(1).expectation.Time(:))
 psit=expm(-1i*Operator(1).H*Operator(1).expectation.Time(d4))*psi0;
 psit=expm(-1i*Operator(1).H*(Operator(1).expectation.Time(2)-...
     Operator(1).expectation.Time(1)))*psit;
 Operator(1).expectation.tau(d4)=psit'*(Operator(1).T)*psit;
 Operator(1).expectation.nu(d4)=psit'*Operator(1).V*psit;
 for d5=1:length(Operator)
     Operator(1).expectation.probability(d4,d5)=psit'*Operator(d5).N*psit;     
 end
Operator(1).expectation.n(d4)=psit'*(Operator(1).Ntotal)*psit;
end
f2=Operator;
end
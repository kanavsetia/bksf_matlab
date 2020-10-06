% =========================================================================
% 					Bravyi-Kitaev Superfast Simulation
%
% Author: 	Kanav Setia*
% Date:     2016/11/08
% Version: 	0.1
% *Whitfield Group, Department of Physics and Astronomy,
%  Dartmouth College
% =========================================================================
% HELP SECTION
% bk_general is the function to run the general bravyi-kitaev algorithm. 
% 
% Input: totalqubits, Hamiltonian (G matrix)
% Output: Data_BK (All the data for the algorithm in a data structure)
% 
% Data Format for the Data structure:
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
% Operator(1).expectation.C--> Stabilizer
% Operator(1).edges--> edge matrix, the order of the edges in this matrix
% also defines the convention for the qubits
% Convention for qubits --> qubits in 1st row--> qubits in 2nd row--> ...
% 
% See also, bj_general, Ajk_general, bkgif.gif
% =========================================================================
function f1=bk_general(totalqubits,Hamiltonian)

% =============================================================================
% Defining the Pauli Matrices
% =============================================================================

x=[0  1 ; 1   0];
y=[0 -1j; 1j  0];
z=[1  0 ; 0  -1];

% =============================================================================
% Defining the Structure for Operators. This Matlab data structure will
% hold all the Operators in one place.
% Operator(i).Ajk operator Ajk for the edge jk. It is a
% tensor which has the sigma z's included to get the right parity.
% Similarly Operator(i).b gives the annhilation operator.
% Operator(1).H is a place holder for the hamiltonian.
% =============================================================================

Operator=[];

% =============================================================================
% Initialization Matrix and the number of qubits.
% =============================================================================


G=Hamiltonian;
edgematrix=triu(G,1);
[e_row,ecol]=find(triu(G,1));
em=[e_row,ecol];
em=sortrows(em,1);
tot_qub=length(e_row);
% =============================================================================
% Initialization of different Operators.
% =============================================================================


Operator(1).T=zeros(2^tot_qub);         %Kinetic Operator
Operator(1).V=zeros(2^tot_qub);         %Potential Operator
Operator(1).Ntotal=zeros(2^tot_qub);    %Number Operator
Operator(1).N=zeros(2^tot_qub);
Operator(1).H=zeros(2^tot_qub);
Operator(1).expectation.Time=linspace(0,10,200);
Operator(1).expectation.tau=[];         %Kinetic Expectation Values (at different times)
Operator(1).expectation.nu=[];          %Potential Expectation Values (at different times)
Operator(1).expectation.probability=[];
Operator(1).expectation.n=[];           %Number Expectation Values (at different times)
Operator(1).expectation.C=[];           %Stabilizer
Operator(1).edges=em;   
% =============================================================================
% Calculation of Operators
% annhilator.m is the function file that calculates the Annhilation
% Operator.
% Note: Annhilation and creation operators are stored in i'th place in the
% structure, where i is the number of the fermion on which we are
% operating. But H, V, and T are stored in the 1st branch of the structure.
% =============================================================================
for d9=1:length(G)
   Operator(d9).B = bj_general(em,d9);
   Operator(1).V=Operator(1).V+G(d9,d9)*(eye(length(Operator(1).B))...
        -Operator(d9).B)/2;
    Operator(d9).N=(eye(length(Operator(d9).B))-Operator(d9).B)/2;
    Operator(1).Ntotal=Operator(1).Ntotal+Operator(d9).N;
end 

for d1=1:tot_qub
    Operator(d1).Ajk = Ajk_general(em,d1);   
end

for d2 =1:length(em)
    
        Operator(1).T=Operator(1).T+-1i/2*G(em(d2,1),em(d2,2))...
            *(Operator(d2).Ajk*Operator(em(d2,2)).B+...
            Operator(em(d2,1)).B  *Operator(d2).Ajk);           
end
Operator(1).H=Operator(1).T+Operator(1).V;



% =============================================================================
% Interaction between two particles.
% =============================================================================


% Operator(1).W=-7.*Operator(1).N*Operator(3).N;
% Operator(1).H=Operator(1).H+Operator(1).W;
% eig(Operator(1).H)

% =============================================================================
% Separating Eigenvalues into sectors
% This section is not required but we used it while developing to  
% understand the problem better.
% =============================================================================
% 
[H_evec,H_eval]=eig(Operator(1).H);
% 
% for d7=1:length(H_eval)    
%         bk_n(d7)=H_evec(:,d7)'*Operator(1).Ntotal*H_evec(:,d7);           
% end


% 
% [N_evecs' diag(H_eval)]
% 
% H_evec_rearranged=[];
% H_eval_rearranged=[];
% H_evec_particles=[];
% [H_evec,d]=qr(H_evec);
% H_evec_particles=diag(H_evec'*Operator(1).Ntotal*H_evec);
% [Hval,Hindex]=sort(H_evec_particles);
% for d8=1:length(Hindex)
%     H_evec_rearranged(:,d8)=H_evec(:,Hindex(d8));
%     H_eval_rearranged(d8)=H_eval(Hindex(d8),Hindex(d8));
% end
% H_evec_particles=Hval
% H_eval_rearranged
% H_evec_rearranged;


% =============================================================================
% State Initialization
% =============================================================================
psi0=zeros(2^tot_qub,1);
% psi0(1)=1;
% =============================================================================
% Defining Stabilizer
% =============================================================================
Operator(1).C=eye(2^tot_qub);
% stabfound=false;
% while(~stabfound)
%     loop=[];
%     loopfound=false;
%     [d10,d11]=find(G~=0);
%     edged=sortrows([d10,d11],1);
%     for d12=1:length(G)
%         d13=find(edged(:,1)==d12);
%         if length(d13)<2
%             edged(d13,:)=[];
%         end                        
%     end
%     while(~loopfound)
%         
%         
%     end
%     
% end
fprintf('\n\n\n\nHelp me pick the loop\n');
em
loop=input('Enter the row numbers: \n  ');




for d3=1:length(loop)
    Operator(1).C = 1i.*Operator(1).C*Operator(loop(d3)).Ajk;
end

% =============================================================================
% Figuring out the vacuum state for Bravyi Kitaev
% =============================================================================


C_eig=eig(Operator(1).C);
[y1,z1]=eig(Operator(1).C);


for d6=1:length(C_eig)
    if C_eig(d6)==1
        np=y1(:,d6)'*Operator(1).Ntotal*y1(:,d6);
        if np==0
            psi0=y1(:,d6); %State Initialized to the Vacuum State
      
        end
    else 
        
    end
end

% psi0=y1(:,7);

d14=input('Row for A23\n');

psi0=-1i/2*(  (Operator(d14).Ajk*Operator(3).B)...
    -(Operator(2).B  *Operator(d14).Ajk) )  * psi0;
psi0=psi0/sqrt(norm(psi0));

% =============================================================================
% Time Evolution of the State
% =============================================================================

psit=psi0;
for d4=1:length(Operator(1).expectation.Time(:))
    psit=expm(-1i*Operator(1).H*(Operator(1).expectation.Time(2)-...
        Operator(1).expectation.Time(1)))*psit;
    Operator(1).expectation.tau(d4)=psit'*(Operator(1).T)*psit;
    Operator(1).expectation.nu(d4)=psit'*Operator(1).V*psit;
    %  b1(k)=psit'*Operator(1).B*psit;
    %  b2(k)=psit'*Operator(2).B*psit;
    %  b3(k)=psit'*Operator(3).B*psit;
    for d5=1:length(G)
        Operator(1).expectation.probability(d4,d5)=psit'*Operator(d5).N*psit;
    end
    
    Operator(1).expectation.n(d4) =psit'*Operator(1).Ntotal*psit;
    
end

f1=Operator;
% plot(t,n)
end

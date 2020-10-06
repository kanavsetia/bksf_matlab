%
%
function y=hydrogen(ordering)

FLAG_PLOT=true;
ORDER_FLAG=2; % 1 is for Seeley's order, 2 is custom

x=[0  1 ; 1   0];
y=[0 -1j; 1j  0];
z=[1  0 ; 0  -1];
i=[1  0 ; 0   1];
tot_qub=4;
G=ones(4,4);
%%% We'll be using 6 qubits
%       34       24          23          14          13          12
B1=kron(i       ,kron(i     ,kron(i     ,kron(z     ,kron(z     ,z  )))));
B2=kron(i       ,kron(z     ,kron(z     ,kron(i     ,kron(i     ,z  )))));
B3=kron(z       ,kron(i     ,kron(z     ,kron(i     ,kron(z     ,i)))));
B4=kron(z       ,kron(z     ,kron(i     ,kron(z     ,kron(i     ,i)))));

%        34           24          23          14          13          12
A12=kron(i      ,kron(i     ,kron(i     ,kron(i     ,kron(i     ,x)))));
A13=kron(i      ,kron(i     ,kron(i     ,kron(i     ,kron(x     ,z)))));
A14=kron(i      ,kron(i     ,kron(i     ,kron(x     ,kron(z     ,z)))));
A23=kron(i      ,kron(i     ,kron(x     ,kron(i     ,kron(z     ,z)))));
A24=kron(i      ,kron(x     ,kron(z     ,kron(z     ,kron(i,    z)))));
A34=kron(x      ,kron(z     ,kron(z     ,kron(z     ,kron(z     ,i)))));
%A34=kron(x     ,kron(i     ,kron(i     ,kron(z     ,kron(z     ,i)))));



%% taking out extra qubits
%
% %       34       24          23          14          13          12
% B1=kron(i       ,kron(i     ,kron(i     ,kron(i     ,kron(i     ,i  )))));
% B2=kron(i       ,kron(z     ,kron(z     ,kron(i     ,kron(i     ,i  )))));
% B3=kron(z       ,kron(i     ,kron(z     ,kron(i     ,kron(i     ,i)))));
% B4=kron(z       ,kron(z     ,kron(i     ,kron(i     ,kron(i     ,i)))));
%
% %        34           24          23          14          13          12
% A12=kron(i      ,kron(i     ,kron(i     ,kron(i     ,kron(i     ,i)))));
% A13=kron(i      ,kron(i     ,kron(i     ,kron(i     ,kron(i     ,i)))));
% A14=kron(i      ,kron(i     ,kron(i     ,kron(i     ,kron(i     ,i)))));
% A23=kron(i      ,kron(i     ,kron(x     ,kron(i     ,kron(i     ,i)))));
% A24=kron(i      ,kron(x     ,kron(z     ,kron(i     ,kron(i,    i)))));
% A34=kron(x      ,kron(z     ,kron(z     ,kron(i     ,kron(i     ,i)))));
% %A34=kron(x     ,kron(i     ,kron(i     ,kron(z     ,kron(z     ,i)))));





V=G(1,1)*(eye(length(B1))-B1)/2 ...
    +G(2,2)*(eye(length(B2))-B2)/2 ...
    +G(3,3)*(eye(length(B3))-B3)/2 ...
    +G(4,4)*(eye(length(B4))-B4)/2;

T=  -1i/2*G(1,2)*(A12*B2+B1*A12)...
    -1i/2*G(1,3)*(A13*B3+B1*A13)...
    -1i/2*G(1,4)*(A14*B4+B1*A14)...
    -1i/2*G(2,3)*(A23*B3+B2*A23)...
    -1i/2*G(2,4)*(A24*B4+B2*A24)...
    -1i/2*G(3,4)*(A34*B4+B3*A34);
H=T+V;
N=   (eye(length(B1))-B1)/2+(eye(length(B2))-B2)/2 ...
    +(eye(length(B3))-B3)/2+(eye(length(B4))-B4)/2;

W=-7.*((eye(length(B1))-B1)/2)*((eye(length(B3))-B3)/2);
H=H+W;





% =========================================================================
%% Hydrogen Hamiltonian
% =========================================================================

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

H_1=h11*(eye(length(B1))-B1)/2+...
    h22*(eye(length(B2))-B2)/2+...
    h33*(eye(length(B3))-B3)/2+...
    h44*(eye(length(B4))-B4)/2;
H_2=   h1221*(eye(length(B1))-B1)/2*(eye(length(B2))-B2)/2+...
    h3443*(eye(length(B3))-B3)/2*(eye(length(B4))-B4)/2+...
    h1441*(eye(length(B1))-B1)/2*(eye(length(B4))-B4)/2+...
    h2332*(eye(length(B2))-B2)/2*(eye(length(B3))-B3)/2+...
    (h1331-h1313)*(eye(length(B1))-B1)/2*(eye(length(B3))-B3)/2+...
    (h2442-h2424)*(eye(length(B2))-B2)/2*(eye(length(B4))-B4)/2;

H_2=H_2+h1243*(1/8)*A12*(-1)*A34*(-eye(64)-B1*B2+B1*B4+B1*B3+B2*B4+B2*B3-B4*B3+B1*B2*B4*B3)+...
    h1423*(1/8)*A14*A23*(-eye(64)-B1*B4+B1*B2+B1*B3+B4*B2+B4*B3-B2*B3+B1*B4*B2*B3);
H=H_1+H_2;
%from openfermion
% H=H+0.71375399.*eye(64);


% =========================================================================
%% Hydrogen Hamiltonian (BKSF) in terms of pauli operators
% =========================================================================
%                            34       24          23          14          13     12
% T_1_c = 0.25*(2*h11+2*h22+2*h33+2*h44+h1221+h3443+h1441+h2332+h1331-h1313+...
%     h2442-h2424)
% T_2_c = 0.25*(-2*h11-h1221-h1441-h1331+h1313)
% T_3_c = 0.25*(-2*h22-h1221-h2332-h2442+h2424)
% T_4_c = 0.25*(-2*h33-h3443-h2332-h1331+h1313)
% T_5_c = 0.25*(-2*h44-h3443-h1441-h2442+h2424)
% T_6_c = 0.25*(h1221+h3443)
% T_7_c = 0.25*(h1441+h2332)
% T_8_c = 0.25*(h1331-h1313+h2442-h2424)
% T_9_c = h1243/4
% 
% 
% T_1_BK = -0.812610*    kron(i  ,kron(i     ,kron(i     ,kron(i     ,kron(i     ,i  )))));
% T_2_BK = 0.171201*     kron(i  ,kron(i     ,kron(i     ,kron(z     ,kron(z     ,z  )))));
% T_3_BK = 0.171201*     kron(i  ,kron(z     ,kron(z     ,kron(i     ,kron(i     ,z  )))));
% T_4_BK = -0.2227965*   kron(z  ,kron(i     ,kron(z     ,kron(i     ,kron(z     ,i  )))));
% T_5_BK = -0.2227965*   kron(z  ,kron(z     ,kron(i     ,kron(z     ,kron(i     ,i  )))));
% T_8_BK = 0.3429725*    kron(i  ,kron(z     ,kron(z     ,kron(z     ,kron(z     ,i  )))));
% T_7_BK = 0.331736*     kron(z  ,kron(z     ,kron(i     ,kron(i     ,kron(z     ,z  )))));
% T_6_BK = 0.2410925*    kron(z  ,kron(i     ,kron(z     ,kron(z     ,kron(i     ,z  )))));
% 
% T_9_BK =  0.04532175*  kron(x  ,kron(i     ,kron(i     ,kron(i     ,kron(i     ,x  )))));
% T_10_BK = 0.04532175*  kron(y  ,kron(i     ,kron(z     ,kron(z     ,kron(i     ,y  )))));
% T_11_BK = 0.04532175*  kron(y  ,kron(z     ,kron(i     ,kron(i     ,kron(z     ,y  )))));
% T_12_BK = -0.04532175*  kron(z  ,kron(z     ,kron(x     ,kron(x     ,kron(z     ,z  )))));
% T_13_BK = -0.04532175*  kron(i  ,kron(z     ,kron(y     ,kron(y     ,kron(z     ,i  )))));
% T_14_BK = -0.04532175*  kron(z  ,kron(i     ,kron(y     ,kron(y     ,kron(i     ,z  )))));

% =========================================================================
%% Different order (BKSF) (current best)
% =========================================================================
%
% %

T_1_BK = -0.812610*    kron(i  ,kron(i     ,kron(i     ,kron(i     ,kron(i     ,i  )))));
% T_1_BK = -0.098863.*eye(64);
%changing the coefficient also changes the ground state, so make sure to
%pick the right ground state for the Hamiltonian.
T_7_BK = 0.171201*     kron(i  ,kron(i     ,kron(i     ,kron(z     ,kron(z     ,z  )))));
T_9_BK = 0.171201*     kron(i  ,kron(z     ,kron(z     ,kron(i     ,kron(i     ,z  )))));
T_3_BK = -0.2227965*   kron(z  ,kron(i     ,kron(z     ,kron(i     ,kron(z     ,i  )))));
T_5_BK = -0.2227965*   kron(z  ,kron(z     ,kron(i     ,kron(z     ,kron(i     ,i  )))));
T_14_BK = 0.3429725*    kron(i  ,kron(z     ,kron(z    ,kron(z     ,kron(z     ,i  )))));
T_13_BK = 0.331736*     kron(z  ,kron(z     ,kron(i    ,kron(i     ,kron(z     ,z  )))));
T_11_BK = 0.2410925*    kron(z  ,kron(i     ,kron(z    ,kron(z     ,kron(i     ,z  )))));



T_2_BK =  0.04532175*  kron(x  ,kron(i     ,kron(i     ,kron(i     ,kron(i     ,x  )))));
T_6_BK =  0.04532175*  kron(y  ,kron(i     ,kron(z      ,kron(z     ,kron(i     ,y  )))));
T_10_BK = 0.04532175*  kron(y  ,kron(z     ,kron(i     ,kron(i     ,kron(z     ,y  )))));
T_4_BK = -0.04532175*  kron(z  ,kron(z     ,kron(x     ,kron(x     ,kron(z     ,z  )))));
T_8_BK = -0.04532175*  kron(i  ,kron(z     ,kron(y     ,kron(y     ,kron(z     ,i  )))));
T_12_BK = -0.04532175*  kron(z  ,kron(i     ,kron(y    ,kron(y     ,kron(i     ,z  )))));




% =========================================================================
%% Different order (BKSF)
% =========================================================================
%

% 
% T_13_BK = -0.812610*    kron(i  ,kron(i     ,kron(i     ,kron(i     ,kron(i     ,i  )))));
% T_7_BK = 0.171201*     kron(i  ,kron(i     ,kron(i     ,kron(z     ,kron(z     ,z  )))));
% T_9_BK = 0.171201*     kron(i  ,kron(z     ,kron(z     ,kron(i     ,kron(i     ,z  )))));
% T_3_BK = -0.2227965*   kron(z  ,kron(i     ,kron(z     ,kron(i     ,kron(z     ,i  )))));
% T_5_BK = -0.2227965*   kron(z  ,kron(z     ,kron(i     ,kron(z     ,kron(i     ,i  )))));
% T_1_BK = 0.3429725*    kron(i  ,kron(z     ,kron(z     ,kron(z     ,kron(z     ,i  )))));
% T_14_BK = 0.331736*     kron(z  ,kron(z     ,kron(i     ,kron(i     ,kron(z     ,z  )))));
% T_11_BK = 0.2410925*    kron(z  ,kron(i     ,kron(z     ,kron(z     ,kron(i     ,z  )))));
% 
% T_4_BK =  0.04532175*  kron(x  ,kron(i     ,kron(i     ,kron(i     ,kron(i     ,x  )))));
% T_8_BK = 0.04532175*  kron(y  ,kron(i     ,kron(z     ,kron(z     ,kron(i     ,y  )))));
% T_12_BK = 0.04532175*  kron(y  ,kron(z     ,kron(i     ,kron(i     ,kron(z     ,y  )))));
% T_2_BK = -0.04532175*  kron(z  ,kron(z     ,kron(x     ,kron(x     ,kron(z     ,z  )))));
% T_6_BK = -0.04532175*  kron(i  ,kron(z     ,kron(y     ,kron(y     ,kron(z     ,i  )))));
% T_10_BK = -0.04532175*  kron(z  ,kron(i     ,kron(y     ,kron(y     ,kron(i     ,z  )))));
% 


% =========================================================================
%% Order for 2nd order Trotterization (BKSF)
% =========================================================================
%

% 
% T_7_BK = -0.812610*    kron(i  ,kron(i     ,kron(i     ,kron(i     ,kron(i     ,i  )))));
% T_8_BK = 0.171201*     kron(i  ,kron(i     ,kron(i     ,kron(z     ,kron(z     ,z  )))));
% T_9_BK = 0.171201*     kron(i  ,kron(z     ,kron(z     ,kron(i     ,kron(i     ,z  )))));
% T_10_BK = -0.2227965*   kron(z  ,kron(i     ,kron(z     ,kron(i     ,kron(z     ,i  )))));
% T_11_BK = -0.2227965*   kron(z  ,kron(z     ,kron(i     ,kron(z     ,kron(i     ,i  )))));
% T_12_BK = 0.3429725*    kron(i  ,kron(z     ,kron(z     ,kron(z     ,kron(z     ,i  )))));
% T_13_BK = 0.331736*     kron(z  ,kron(z     ,kron(i     ,kron(i     ,kron(z     ,z  )))));
% T_14_BK = 0.2410925*    kron(z  ,kron(i     ,kron(z     ,kron(z     ,kron(i     ,z  )))));
% 
% T_1_BK =  0.04532175*  kron(x  ,kron(i     ,kron(i     ,kron(i     ,kron(i     ,x  )))));
% T_2_BK = 0.04532175*  kron(y  ,kron(i     ,kron(z     ,kron(z     ,kron(i     ,y  )))));
% T_3_BK = 0.04532175*  kron(y  ,kron(z     ,kron(i     ,kron(i     ,kron(z     ,y  )))));
% T_4_BK = -0.04532175*  kron(z  ,kron(z     ,kron(x     ,kron(x     ,kron(z     ,z  )))));
% T_5_BK = -0.04532175*  kron(i  ,kron(z     ,kron(y     ,kron(y     ,kron(z     ,i  )))));
% T_6_BK = -0.04532175*  kron(z  ,kron(i     ,kron(y     ,kron(y     ,kron(i     ,z  )))));
% 
% 
% 

% =========================================================================
%% Hydrogen Hamiltonian (Jordan-Wigner) in terms of pauli operators
% =========================================================================
% %
% T_1_JW = -0.81261.*eye(16);
% T_2_JW = 0.171201.*     kron(i  ,kron(i,    kron(i  ,z)));
% T_3_JW = 0.171201.*     kron(i  ,kron(i,    kron(z  ,i)));
% T_4_JW = -0.2227965.*   kron(i  ,kron(z,    kron(i  ,i)));
% T_5_JW = -0.2227965.*   kron(z  ,kron(i,    kron(i  ,i)));
% T_6_JW = 0.16862325.*   kron(i  ,kron(i,    kron(z  ,z)));
% T_7_JW = 0.12054625.*   kron(i  ,kron(z,    kron(i  ,z)));
% T_8_JW = 0.165868.*     kron(i  ,kron(z,    kron(z  ,i)));
% T_9_JW = 0.165868.*     kron(z  ,kron(i,    kron(i  ,z)));
% T_10_JW = 0.12054625.*  kron(z  ,kron(i,    kron(z  ,i)));
% T_11_JW = 0.17434925.*  kron(z  ,kron(z,    kron(i  ,i)));
% 
% T_12_JW = -0.04532175.* kron(x  ,kron(x,    kron(y  ,y)));
% T_13_JW = 0.04532175.*  kron(x  ,kron(y,    kron(y  ,x)));
% T_14_JW = 0.04532175.*  kron(y  ,kron(x,    kron(x  ,y)));
% T_15_JW = -0.04532175.* kron(y  ,kron(y,    kron(x  ,x)));


if ORDER_FLAG==1
% % =========================================================================
% % Hydrogen Hamiltonian (Jordan-Wigner) (Best Order)
% % =========================================================================

T_1_JW = -0.81261.*eye(16);
T_13_JW = 0.171201.*     kron(i  ,kron(i,    kron(i  ,z)));
T_14_JW = 0.171201.*     kron(i  ,kron(i,    kron(z  ,i)));
T_3_JW = -0.2227965.*   kron(i  ,kron(z,    kron(i  ,i)));
T_5_JW = -0.2227965.*   kron(z  ,kron(i,    kron(i  ,i)));
T_12_JW = 0.16862325.*   kron(i  ,kron(i,    kron(z  ,z)));
T_7_JW = 0.12054625.*   kron(i  ,kron(z,    kron(i  ,z)));
T_10_JW = 0.165868.*     kron(i  ,kron(z,    kron(z  ,i)));
T_11_JW = 0.165868.*     kron(z  ,kron(i,    kron(i  ,z)));
T_9_JW = 0.12054625.*  kron(z  ,kron(i,    kron(z  ,i)));
T_15_JW = 0.17434925.*  kron(z  ,kron(z,    kron(i  ,i)));

T_6_JW = -0.04532175.* kron(x  ,kron(x,    kron(y  ,y)));
T_2_JW = 0.04532175.*  kron(x  ,kron(y,    kron(y  ,x)));
T_4_JW = 0.04532175.*  kron(y  ,kron(x,    kron(x  ,y)));
T_8_JW = -0.04532175.* kron(y  ,kron(y,    kron(x  ,x)));

elseif ORDER_FLAG==2
% =========================================================================
% Hydrogen Hamiltonian (Jordan-Wigner) (Better than Peter's paper)
% =========================================================================
%
T_13_JW = -0.81261.*eye(16);
% T_13_JW = -0.098863.*eye(16);
T_1_JW = 0.171201.*     kron(i  ,kron(i,    kron(i  ,z)));
T_14_JW = 0.171201.*     kron(i  ,kron(i,    kron(z  ,i)));
T_3_JW = -0.2227965.*   kron(i  ,kron(z,    kron(i  ,i)));
T_5_JW = -0.2227965.*   kron(z  ,kron(i,    kron(i  ,i)));
T_12_JW = 0.16862325.*   kron(i  ,kron(i,    kron(z  ,z)));
T_7_JW = 0.12054625.*   kron(i  ,kron(z,    kron(i  ,z)));
T_10_JW = 0.165868.*     kron(i  ,kron(z,    kron(z  ,i)));
T_11_JW = 0.165868.*     kron(z  ,kron(i,    kron(i  ,z)));
T_9_JW = 0.12054625.*  kron(z  ,kron(i,    kron(z  ,i)));
T_15_JW = 0.17434925.*  kron(z  ,kron(z,    kron(i  ,i)));

T_2_JW = -0.04532175.* kron(x  ,kron(x,    kron(y  ,y)));
T_4_JW = 0.04532175.*  kron(x  ,kron(y,    kron(y  ,x)));
T_6_JW = 0.04532175.*  kron(y  ,kron(x,    kron(x  ,y)));
T_8_JW = -0.04532175.* kron(y  ,kron(y,    kron(x  ,x)));

else
end



% =========================================================================
%% Order for 2nd order trotterization (Jordan-Wigner)
% =========================================================================
% %
% T_5_JW = -0.81261.*eye(16);
% T_6_JW = 0.171201.*     kron(i  ,kron(i,    kron(i  ,z)));
% T_7_JW = 0.171201.*     kron(i  ,kron(i,    kron(z  ,i)));
% T_8_JW = -0.2227965.*   kron(i  ,kron(z,    kron(i  ,i)));
% T_9_JW = -0.2227965.*   kron(z  ,kron(i,    kron(i  ,i)));
% T_10_JW = 0.16862325.*   kron(i  ,kron(i,    kron(z  ,z)));
% T_11_JW = 0.12054625.*   kron(i  ,kron(z,    kron(i  ,z)));
% T_12_JW = 0.165868.*     kron(i  ,kron(z,    kron(z  ,i)));
% T_13_JW = 0.165868.*     kron(z  ,kron(i,    kron(i  ,z)));
% T_14_JW = 0.12054625.*  kron(z  ,kron(i,    kron(z  ,i)));
% T_15_JW = 0.17434925.*  kron(z  ,kron(z,    kron(i  ,i)));
% 
% T_1_JW = -0.04532175.* kron(x  ,kron(x,    kron(y  ,y)));
% T_2_JW = 0.04532175.*  kron(x  ,kron(y,    kron(y  ,x)));
% T_3_JW = 0.04532175.*  kron(y  ,kron(x,    kron(x  ,y)));
% T_4_JW = -0.04532175.* kron(y  ,kron(y,    kron(x  ,x)));
% 




% =========================================================================
%% Jordan Wigner check
% =========================================================================

string_H_JW='';
for iter=1:15
    string_H_JW=strcat(strcat(strcat(string_H_JW,'(T_'),num2str(iter)),'_JW)');
    if iter~=15
        string_H_JW=strcat(string_H_JW,'+');
    end
end
H_JW=eval(string_H_JW);
[H_JW_evec,H_JW_eval]=eig(H_JW);
psi_gs_JW=H_JW_evec(:,1);

% =========================================================================
%% Hydrogen Hamiltonian (Bravyi-Kitaev)
% =========================================================================
% %
% T_1_KB = -0.81261.*      eye(16);
% T_2_KB = 0.171201.*     kron(i  ,kron(i,    kron(i  ,z)));
% T_3_KB = 0.171201.*     kron(i  ,kron(i,    kron(z  ,z)));
% T_4_KB = -0.2227965.*    kron(i  ,kron(z,    kron(i  ,i)));
% T_5_KB = -0.2227965.*    kron(z  ,kron(z,    kron(z  ,i)));
% T_6_KB = 0.16862325.*   kron(i  ,kron(i,    kron(z  ,i)));
% T_7_KB = 0.12054625.*    kron(i  ,kron(z,    kron(i  ,z)));
% T_8_KB = 0.165868.*     kron(i  ,kron(z,    kron(z  ,z)));
% T_9_KB = 0.165868.*     kron(z  ,kron(z,    kron(z  ,z)));
% T_10_KB =  0.12054625.*   kron(z  ,kron(z,    kron(i  ,z)));
% T_11_KB = 0.17434925.*   kron(z  ,kron(i,    kron(z  ,i)));
% 
% 
% T_12_KB = 0.04532175.*    kron(i  ,kron(x,    kron(z  ,x)));
% T_13_KB = 0.04532175.*    kron(i  ,kron(y,    kron(z  ,y)));
% T_14_KB = 0.04532175.*    kron(z  ,kron(x,    kron(z  ,x)));
% T_15_KB = 0.04532175.*    kron(z  ,kron(y,    kron(z  ,y)));



if ORDER_FLAG==1
% =========================================================================
%% Hydrogen Hamiltonian (Bravyi-Kitaev) (Best)
% =========================================================================

T_1_KB = -0.81261*       eye(16);
T_13_KB = 0.171201.*     kron(i  ,kron(i,    kron(i  ,z)));
T_14_KB = 0.171201.*     kron(i  ,kron(i,    kron(z  ,z)));
T_3_KB = -0.2227965.*    kron(i  ,kron(z,    kron(i  ,i)));
T_5_KB = -0.2227965.*    kron(z  ,kron(z,    kron(z  ,i)));
T_12_KB = 0.16862325.*   kron(i  ,kron(i,    kron(z  ,i)));
T_7_KB = 0.12054625.*    kron(i  ,kron(z,    kron(i  ,z)));
T_10_KB = 0.165868.*     kron(i  ,kron(z,    kron(z  ,z)));
T_11_KB = 0.165868.*     kron(z  ,kron(z,    kron(z  ,z)));
T_9_KB =  0.12054625.*   kron(z  ,kron(z,    kron(i  ,z)));
T_15_KB = 0.17434925.*   kron(z  ,kron(i,    kron(z  ,i)));


T_6_KB = 0.04532175.*    kron(i  ,kron(x,    kron(z  ,x)));
T_2_KB = 0.04532175.*    kron(i  ,kron(y,    kron(z  ,y)));
T_4_KB = 0.04532175.*    kron(z  ,kron(x,    kron(z  ,x)));
T_8_KB = 0.04532175.*    kron(z  ,kron(y,    kron(z  ,y)));
% 


elseif ORDER_FLAG==2
%
% =========================================================================
% Hydrogen Hamiltonian (Bravyi-Kitaev) (better than peter's paper)
% =========================================================================

T_13_KB = -0.81261*      eye(16);
% T_13_KB = -0.098863*      eye(16);
T_1_KB = 0.171201.*     kron(i  ,kron(i,    kron(i  ,z)));
T_14_KB = 0.171201.*     kron(i  ,kron(i,    kron(z  ,z)));
T_3_KB = -0.2227965.*    kron(i  ,kron(z,    kron(i  ,i)));
T_5_KB = -0.2227965.*    kron(z  ,kron(z,    kron(z  ,i)));
T_12_KB = 0.16862325.*   kron(i  ,kron(i,    kron(z  ,i)));
T_7_KB = 0.12054625.*    kron(i  ,kron(z,    kron(i  ,z)));
T_10_KB = 0.165868.*     kron(i  ,kron(z,    kron(z  ,z)));
T_11_KB = 0.165868.*     kron(z  ,kron(z,    kron(z  ,z)));
T_9_KB =  0.12054625.*   kron(z  ,kron(z,    kron(i  ,z)));
T_15_KB = 0.17434925.*   kron(z  ,kron(i,    kron(z  ,i)));


T_2_KB = 0.04532175.*    kron(i  ,kron(x,    kron(z  ,x)));
T_4_KB = 0.04532175.*    kron(i  ,kron(y,    kron(z  ,y)));
T_6_KB = 0.04532175.*    kron(z  ,kron(x,    kron(z  ,x)));
T_8_KB = 0.04532175.*    kron(z  ,kron(y,    kron(z  ,y)));


else
end

%
% =========================================================================
%% Order for 2nd order trotterization (Bravyi-Kitaev)
% =========================================================================
% 
% T_5_KB = -0.81261*      eye(16);
% T_6_KB = 0.171201.*     kron(i  ,kron(i,    kron(i  ,z)));
% T_7_KB = 0.171201.*     kron(i  ,kron(i,    kron(z  ,z)));
% T_8_KB = -0.2227965.*    kron(i  ,kron(z,    kron(i  ,i)));
% T_9_KB = -0.2227965.*    kron(z  ,kron(z,    kron(z  ,i)));
% T_10_KB = 0.16862325.*   kron(i  ,kron(i,    kron(z  ,i)));
% T_11_KB = 0.12054625.*    kron(i  ,kron(z,    kron(i  ,z)));
% T_12_KB = 0.165868.*     kron(i  ,kron(z,    kron(z  ,z)));
% T_13_KB = 0.165868.*     kron(z  ,kron(z,    kron(z  ,z)));
% T_14_KB =  0.12054625.*   kron(z  ,kron(z,    kron(i  ,z)));
% T_15_KB = 0.17434925.*   kron(z  ,kron(i,    kron(z  ,i)));
% 
% 
% T_1_KB = 0.04532175.*    kron(i  ,kron(x,    kron(z  ,x)));
% T_2_KB = 0.04532175.*    kron(i  ,kron(y,    kron(z  ,y)));
% T_3_KB = 0.04532175.*    kron(z  ,kron(x,    kron(z  ,x)));
% T_4_KB = 0.04532175.*    kron(z  ,kron(y,    kron(z  ,y)));







ordering=1:15;
string_H_KB='';
for iter=1:15
    string_H_KB=strcat(strcat(strcat(string_H_KB,'(T_'),num2str(ordering(iter))),'_KB)');
    if iter~=15
        string_H_KB=strcat(string_H_KB,'+');
    end
end
H_KB=eval(string_H_KB);

%
% H_KB = 0.16586702396410952.* kron(z  ,kron(z,    kron(z  ,z)))+ ...
%        0.17434844170557137.* kron(z  ,kron(i,    kron(z  ,i)))+ ...
%        0.045322202098565384.*kron(i  ,kron(x,    kron(z  ,x)))+ ...
%        0.045322202098565384.*kron(z  ,kron(x,    kron(z  ,x)))+ ...
%        0.17119774853325864.*kron(i  ,kron(i,    kron(i  ,z)))+ ...
%        0.12054482186554413.*kron(z  ,kron(z,    kron(i  ,z))) ...
%        -0.22278592890107018.*kron(z  ,kron(z,    kron(z  ,i))) ...
%        -0.22278592890107018.*kron(i  ,kron(z,    kron(i  ,i)))+ ...
%        0.16586702396410952.*kron(i  ,kron(z,    kron(z  ,z))) + ...
%        0.045322202098565384.*kron(z  ,kron(y,    kron(z  ,y)))+ ...
%        0.12054482186554413.*kron(i  ,kron(z,    kron(i  ,z)))+ ...
%        0.1686221914334755.*kron(i  ,kron(i,    kron(z  ,i)))+ ...
%        -0.81261.*eye(16) + ...
%        0.045322202098565384.*kron(i  ,kron(y,    kron(z  ,y)))+ ...
%        0.1711977485332587.*kron(i  ,kron(i,    kron(z  ,z)));




[H_KB_evec,H_KB_eval]=eig(H_KB);
psi_gs_KB=H_KB_evec(:,1);

% keyboard


%  stab = -A12*A14*A23*A34;
stab= 1i*A23*A34*A24;
stab1=1i*A12*A23*A13;
stab2=1i*A12*A24*A14;
stab3=1i*A13*A34*A14;
stab4=-A12*A23*A34*A14;


[y1,z1]=eig(stab);
stabsub=y1(:,diag(z1)==1);

stab1projection= stabsub'*stab1*stabsub;
[y2,z2]=eig(stab1projection);
stab1sub=y2(:,round(1e5*diag(z2))/1e5==1);


stab2projection=stab1sub'*stabsub'*stab2*stabsub*stab1sub;
[y3,z3]=eig(stab2projection);
stab2sub=y3(:,round(1e5*diag(z3))/1e5==1);

stab3projection=stab2sub'*stab1sub'*stabsub'*stab3*stabsub*stab1sub*stab2sub;
[y4,z4]=eig(stab3projection);
stab3sub=y4(:,round(1e5*diag(z4))/1e5==1);

stab4projection=stab3sub'*stab2sub'*stab1sub'*stabsub'*stab4...
    *stabsub*stab1sub*stab2sub*stab3sub;
[y5,z5]=eig(stab4projection);
stab4sub=y5(:,round(1e5*diag(z5))/1e5==1);

stab5projection=stab4sub'*stab3sub'*stab2sub'*stab1sub'*stabsub'*stab4...
    *stabsub*stab1sub*stab2sub*stab3sub*stab4sub;
[y6,z6]=eig(stab5projection);
stab5sub=y6(:,round(1e5*diag(z6))/1e5==1);

vac=[];
psi0=[];
psi1=[];
np=[];
for d7=1:size(stab3sub,2)
    v=stabsub*stab1sub*stab2sub*stab3sub(:,d7);
    np(d7)=v'*N*v;
    if round(1e5*np(d7))/1e5==0
        psi0=v; %State Initialized to the Vacuum State
        vac=v;
    end
end

% =========================================================================
%% Using Bj's to define the initial Vacuum State.
% =========================================================================

% [y1,z1]=eig(B1);
% B1sub=y1(:,diag(z1)==1);
%
%
% B2projection= B1sub'*B2*B1sub;
% [y2,z2]=eig(B2projection);
% B2sub=y2(:,round(1e5*diag(z2))/1e5==1);
%
% B3projection= B2sub'*B1sub'*B3*B1sub*B2sub;
% [y3,z3]=eig(B3projection);
% B3sub=y3(:,round(1e5*diag(z3))/1e5==1);
%
% B4projection= B3sub'*B2sub'*B1sub'*B4*B1sub*B2sub*B3sub;
% [y4,z4]=eig(B4projection);
% B4sub=y4(:,round(1e5*diag(z4))/1e5==1);
%
%
%
% k=B1sub*B2sub*B3sub*B4sub;
% psi0=k(:,6)/norm(k(:,6));
% keyboard
% 12
% 13
% 14
% 24
% 23
% 34
% coef=[1;0;0;0;0;0];
%
% psi0=(-1i/2*coef(5)*((A23*B3)...
%     -(B2*A23) )-1i/2*coef(3)*((A14*B4)...
%     -(B1*A14) )-1i/2*coef(4)*((A24*B4)...
%     -(B2*A24) )-1i/2*coef(6)*((A34*B4)...
%     -(B3*A34) )-1i/2*coef(1)*((A12*B2)...
%     -(B1*A12) )-1i/2*coef(2)*((A13*B3)...
%     -(B1*A13) ))* psi0;
%
% psi0=-1i/2*((A23*B3)...
%     -(B2*A23) )* psi0;

% =========================================================================
%% Initializing ground state
% The ground state is calculated from Jordan Wigner diagonalization and the
% state will be constructed by initializing individual components and then
% taking superposition of those components.
% =========================================================================

psi_gs_1=-1i/2*((A12*B2)...
    -(B1*A12) )* vac;
psi_gs_1=psi_gs_1/sqrt(norm(psi_gs_1));
psi_gs_2=-1i/2*((A34*B4)...
    -(B3*A34) )* vac;
psi_gs_2=psi_gs_2/sqrt(norm(psi_gs_2));

psi_gs=0.9936.*psi_gs_1-0.1128.*psi_gs_2;


% nu=psi0'*(V)*psi0
% tau=psi0'*(T)*psi0

tau=[];
nu=[];
prob=[];
num=[];
t=linspace(0,100,2000);
%  psit=psi0;
psit=psi_gs;
data=[];
data_JW=[];
data_KB=[];
data_BK=[];
steps=8;
for n=1:steps
    string_H_JW_exp='';
    for iter=1:15
        string_H_JW_exp=strcat(strcat(strcat(string_H_JW_exp,'expm(-1j.*T_'),num2str(iter)),'_JW*1/n)');
        if iter~=15
            string_H_JW_exp=strcat(string_H_JW_exp,'*');
        end
    end
    %     'Jordan Wigner Eval'
    data_JW(n,1)=angle(psi_gs_JW'*(eval(string_H_JW_exp))^n*psi_gs_JW);
    data(n,1)=angle(psi_gs_JW'*expm(-1j.*H_JW)*psi_gs_JW)-angle(psi_gs_JW'*(eval(string_H_JW_exp))^n*psi_gs_JW);
    
    
    
    string_H_KB_exp='';
    for iter=1:15
        string_H_KB_exp=strcat(strcat(strcat(string_H_KB_exp,'expm(-1j.*T_'),num2str(iter)),'_KB*1/n)');
        if iter~=15
            string_H_KB_exp=strcat(string_H_KB_exp,'*');
        end
    end
    data(n,2)=angle(psi_gs_KB'*expm(-1j.*H_KB)*psi_gs_KB)-angle(psi_gs_KB'*(eval(string_H_KB_exp))^n*psi_gs_KB);
    data_KB(n,2)=angle(psi_gs_KB'*(eval(string_H_KB_exp))^n*psi_gs_KB);
    
    
    
    string_H_BK_exp='';
    for iter=1:14
        string_H_BK_exp=strcat(strcat(strcat(string_H_BK_exp,'expm(-1j.*T_'),num2str(iter)),'_BK*1/n)');
        if iter~=14
            string_H_BK_exp=strcat(string_H_BK_exp,'*');
        end
    end
    %     'BKSF Eval'
    data(n,3)=angle(psi_gs'*expm(-1j.*H)*psi_gs)-angle(psi_gs'*(eval(string_H_BK_exp))^n*psi_gs);
    data_BK(n,3)=angle(psi_gs'*(eval(string_H_BK_exp))^n*psi_gs);
end
data1=[];
data1_JW=[];
data1_KB=[];
data1_BK=[];
for n=1:8
    string_H_JW_exp='';
    for iter=1:15
        string_H_JW_exp=strcat(strcat(strcat(string_H_JW_exp,'expm(-1j.*T_'),num2str(16-iter)),'_JW*1/n)');
        if (16-iter)~=1
            string_H_JW_exp=strcat(string_H_JW_exp,'*');
        end
    end
    %     'Jordan Wigner Eval'
    data1_JW(n,1)=angle(psi_gs_JW'*(eval(string_H_JW_exp))^n*psi_gs_JW);
    data1(n,1)=angle(psi_gs_JW'*expm(-1j.*H_JW)*psi_gs_JW)-angle(psi_gs_JW'*(eval(string_H_JW_exp))^n*psi_gs_JW);
    
    
    
    string_H_KB_exp='';
    for iter=1:15
        string_H_KB_exp=strcat(strcat(strcat(string_H_KB_exp,'expm(-1j.*T_'),num2str(16-iter)),'_KB*1/n)');
        if (16-iter)~=1
            string_H_KB_exp=strcat(string_H_KB_exp,'*');
        end
    end
    data1(n,2)=angle(psi_gs_KB'*expm(-1j.*H_KB)*psi_gs_KB)-angle(psi_gs_KB'*(eval(string_H_KB_exp))^n*psi_gs_KB);
    data1_KB(n,1)=angle(psi_gs_KB'*(eval(string_H_KB_exp))^n*psi_gs_KB);
    
    
    
    string_H_BK_exp='';
    for iter=1:14
        string_H_BK_exp=strcat(strcat(strcat(string_H_BK_exp,'expm(-1j.*T_'),num2str(15-iter)),'_BK*1/n)');
        if (15-iter)~=1
            string_H_BK_exp=strcat(string_H_BK_exp,'*');
        end
    end
    %     'BKSF Eval'
    data1(n,3)=angle(psi_gs'*expm(-1j.*H)*psi_gs)-angle(psi_gs'*(eval(string_H_BK_exp))^n*psi_gs);
    data1_BK(n,1)=angle(psi_gs'*(eval(string_H_BK_exp))^n*psi_gs);
end
xaxis=linspace(1,steps,steps);
if (FLAG_PLOT)
    s_JW = plot(xaxis, -1.*data_JW(:,1), 'k','MarkerSize',7, 'LineWidth',2);
    hold on
    s_JW.Marker = 'd';
%     s_JW1 = plot(xaxis, -1.*data_JW(:,1), 'k','MarkerSize',7, 'LineWidth',2);
    hold on
    s_JW1.Marker = '^';
    s_KB = plot(xaxis, -1.*data_KB(:,2), '-b','MarkerSize',2, 'LineWidth',2);
    s_KB.Marker = '*';
    s_BK = plot(xaxis, -1.*data_BK(:,3), ':r','MarkerSize',7, 'LineWidth',2);
    s_BK.Marker = 's';
    ax=gca;
    ax.XLim=[0 9];
    ax.YLim=[-1.8513 -1.84945];
    xlabel('Trotter steps');
    ylabel('Energy (atomic units)');
    plot([0,9],-1.*[angle(psi_gs_JW'*expm(-1j.*H_JW)*psi_gs_JW),angle(psi_gs_JW'*expm(-1j.*H_JW)*psi_gs_JW)],'-.k');
    plot([0,9],-1.*[angle(psi_gs_JW'*expm(-1j.*H_JW)*psi_gs_JW)-0.0001,angle(psi_gs_JW'*expm(-1j.*H_JW)*psi_gs_JW)-0.0001],':r');
    plot([0,9],-1.*[angle(psi_gs_JW'*expm(-1j.*H_JW)*psi_gs_JW)+0.0001,angle(psi_gs_JW'*expm(-1j.*H_JW)*psi_gs_JW)+0.0001],':r');
    legend('Jordan-Wigner','Bravyi-Kitaev','BKSF');
%     legend([s_JW,s_BK,s_JW1],'JW, BK','BKSF','JW1, BK1 ', 'Orientation', 'horizontal')
end

% % =========================================================================
%% Second-order Trotterization
% % =========================================================================
% 
% 
% data_T2=[];
% data_JW_T2=[];
% data_KB_T2=[];
% data_BK_T2=[];
% for n=1:7
%     string_H_JW_exp_T2='';
%     for iter=1:19
%         if iter<5
%             string_H_JW_exp_T2=strcat(strcat(strcat(string_H_JW_exp_T2,'expm(-1j.*T_'),num2str(iter)),'_JW*1/(2*n))');
%             string_H_JW_exp_T2=strcat(string_H_JW_exp_T2,'*');
%         elseif iter<16
%             string_H_JW_exp_T2=strcat(strcat(strcat(string_H_JW_exp_T2,'expm(-1j.*T_'),num2str(iter)),'_JW*1/(n))');
%             string_H_JW_exp_T2=strcat(string_H_JW_exp_T2,'*');
%         else
%             string_H_JW_exp_T2=strcat(strcat(strcat(string_H_JW_exp_T2,'expm(-1j.*T_'),num2str(20-iter)),'_JW*1/(2*n))');
%             if iter~=19
%             string_H_JW_exp_T2=strcat(string_H_JW_exp_T2,'*');
%             end
%         end
%         
%         
%     end
%     %     'Jordan Wigner Eval'
%     data_JW_T2(n,1)=angle(psi_gs_JW'*(eval(string_H_JW_exp_T2))^n*psi_gs_JW);
%     data_T2(n,1)=angle(psi_gs_JW'*expm(-1j.*H_JW)*psi_gs_JW)-angle(psi_gs_JW'*(eval(string_H_JW_exp_T2))^n*psi_gs_JW);
%     
%     
%     
%     string_H_KB_exp_T2='';
%     for iter=1:19        
%         if iter<5
%             string_H_KB_exp_T2=strcat(strcat(strcat(string_H_KB_exp_T2,'expm(-1j.*T_'),num2str(iter)),'_KB*1/(2*n))');
%             string_H_KB_exp_T2=strcat(string_H_KB_exp_T2,'*');
%         elseif iter<16
%             string_H_KB_exp_T2=strcat(strcat(strcat(string_H_KB_exp_T2,'expm(-1j.*T_'),num2str(iter)),'_KB*1/n)');
%             string_H_KB_exp_T2=strcat(string_H_KB_exp_T2,'*');
%         else
%             string_H_KB_exp_T2=strcat(strcat(strcat(string_H_KB_exp_T2,'expm(-1j.*T_'),num2str(20-iter)),'_KB*1/(2*n))');
%             if iter~=19
%             string_H_KB_exp_T2=strcat(string_H_KB_exp_T2,'*');
%             end        
%         end
%     end
%     data_T2(n,2)=angle(psi_gs_KB'*expm(-1j.*H_KB)*psi_gs_KB)-angle(psi_gs_KB'*(eval(string_H_KB_exp_T2))^n*psi_gs_KB);
%     data_KB_T2(n,2)=angle(psi_gs_KB'*(eval(string_H_KB_exp_T2))^n*psi_gs_KB);
%     
%     
%     
%     string_H_BK_exp_T2='';
%     for iter=1:20
%         if iter<7
%             string_H_BK_exp_T2=strcat(strcat(strcat(string_H_BK_exp_T2,'expm(-1j.*T_'),num2str(iter)),'_BK*1/(2*n))');
%             string_H_BK_exp_T2=strcat(string_H_BK_exp_T2,'*');
%         elseif iter<15
%             string_H_BK_exp_T2=strcat(strcat(strcat(string_H_BK_exp_T2,'expm(-1j.*T_'),num2str(iter)),'_BK*1/n)');
%             string_H_BK_exp_T2=strcat(string_H_BK_exp_T2,'*');
%         else
%             string_H_BK_exp_T2=strcat(strcat(strcat(string_H_BK_exp_T2,'expm(-1j.*T_'),num2str(21-iter)),'_BK*1/(2*n))');
%             if iter~=20
%             string_H_BK_exp_T2=strcat(string_H_BK_exp_T2,'*');
%             end        
%         end
%     end
%     %     'BKSF Eval'
%     data_T2(n,3)=angle(psi_gs'*expm(-1j.*H)*psi_gs)-angle(psi_gs'*(eval(string_H_BK_exp_T2))^n*psi_gs);
%     data_BK_T2(n,3)=angle(psi_gs'*(eval(string_H_BK_exp_T2))^n*psi_gs);
% end
% 
% xaxis=linspace(1,7,7);
% if (FLAG_PLOT)
%     figure(2);
%     s_JW_T2 = plot(xaxis, -1.*data_JW_T2(:,1), 'r');
%     hold on
%     s_JW_T2.Marker = 'o';
%     s_KB_T2 = plot(xaxis, -1.*data_KB_T2(:,2), 'g');
%     s_KB_T2.Marker = 'o';
%     s_BK_T2 = plot(xaxis, -1.*data_BK_T2(:,3), 'b');
%     s_BK_T2.Marker = 'o';
%     ax=gca;
% %     ax.XLim=[0 9];
% %     ax.YLim=[-1.8513 -1.84945];
%     xlabel('Trotter steps');
%     ylabel('Energy (atomic units)');
%     plot([0,9],-1.*[angle(psi_gs_JW'*expm(-1j.*H_JW)*psi_gs_JW),angle(psi_gs_JW'*expm(-1j.*H_JW)*psi_gs_JW)],'--r');
%     plot([0,9],-1.*[angle(psi_gs_JW'*expm(-1j.*H_JW)*psi_gs_JW)-0.0001,angle(psi_gs_JW'*expm(-1j.*H_JW)*psi_gs_JW)-0.0001],'--b');
%     plot([0,9],-1.*[angle(psi_gs_JW'*expm(-1j.*H_JW)*psi_gs_JW)+0.0001,angle(psi_gs_JW'*expm(-1j.*H_JW)*psi_gs_JW)+0.0001],'--b');
%     legend('Jordan-Wigner','Bravyi-Kitaev','BKSF','Correct Energy', 'upper bound', 'lower bound');
% end







    for d4=1:length(t)
        psit=expm(-1i*H*(t(2)-t(1)))*psit;
        tau(d4)=psit'*(T)*psit;
        nu(d4)=psit'*V*psit;
        
        prob(d4,1)=psit'*((eye(length(B1))-B1)/2)*psit;
        prob(d4,2)=psit'*((eye(length(B2))-B2)/2)*psit;
        prob(d4,3)=psit'*((eye(length(B3))-B3)/2)*psit;
        prob(d4,4)=psit'*((eye(length(B4))-B4)/2)*psit;
        
        
        num(d4) =psit'*N*psit;
    end
    keyboard;
    
    
    y=prob;
end

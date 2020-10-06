x=[0  1 ; 1   0];
y=[0 -1j; 1j  0];
z=[1  0 ; 0  -1];
i=[1  0 ; 0   1];
 
A=(x+1j*y)/2;
%JW raising lowering ops
a1=kron(kron(A,i),i);
a2=kron(kron(z,A),i);
a3=kron(kron(z,z),A);
 
H_full=a1'*a2+a2'*a3+a3'*a1+...
  a2'*a1+a3'*a2+a1'*a3;
 
eig(H_full)
 
 
%A=12
%B=23
%C=13
%             A B  C
Y12=kron(kron(y,i),i);
Y23=kron(kron(i,y),i);
Y13=kron(kron(i,i),y);
%             A B  C
Z12=kron(kron(z,i),i);
Z23=kron(kron(i,z),i);
Z13=kron(kron(i,i),z);
 
H_BK2=.5*Y12*(Z23-Z13)+.5*Y23*(Z13-Z12)+.5*Y13*Z12*(Z12-Z23);
eig(H_BK2)

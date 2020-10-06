id=eye(2);
x=[0  1 ; 1   0];
y=[0 -1j; 1j  0];
z=[1  0 ; 0  -1];
psi=zeros(4,1);
psi(1)=1;

ak1=annhilator(2,1);
ak2=annhilator(2,2);
akd1=annhilator(2,1)';
akd2=annhilator(2,2)';
c0=ak1+akd1;
c1=1/1i*(ak1-akd1);
c2=ak2+akd2;
c3=1/1i*(ak2-akd2);

% c0=kron(eye(8),x);
% c1=kron(kron(eye(4),x),eye(2));
% c2=kron(kron(eye(2),x),eye(4));
% c3=kron(x,eye(8));



B0=-1i*c0*c1;
B1=-1i*c2*c3;
A01=-1i*c0*c2;


% B0=kron(eye(4),z);
% B1=kron(kron(id,z),id);
% B2=kron(z,eye(4));
% A01=-kron(eye(2),kron(x,y));
% A12=-1i*kron(eye(2),kron(x,z))*kron(x,kron(z,z));
% A02=-1i*kron(eye(4),x)*kron(x,eye(4));
% A01*B0*A01*psi
function y=Ajk(G,j,k)
x=[0  1 ; 1   0];
y=[0 -1j; 1j  0];
z=[1  0 ; 0  -1];
epsjk=0;
if j>k && H(j,k)~=0
    epsjk=1;
elseif k<j && H(j,k)~=0
    epsjk=-1;
end

% edge numbering goes like (34)-->3 (41)-->4
% (23)-->2 (12)-->1



cre_dum=eye;
if quj-1>0
    for i=1:quj-1
        cre_dum=kron(cre_dum,z);
%         cre_dum=kron(z,cre_dum);
    end
    cre_dum=kron(cre_dum,q);
    if tot_qub-quj>0
        for i=quj+1:tot_qub
            cre_dum=kron(cre_dum,eye(2));
        end
    end
else
    cre_dum=kron(cre_dum,q);
    if tot_qub>1
        for i=2 :tot_qub
            cre_dum=kron(cre_dum,eye(2));
        end
    end
end

y=cre_dum;
end

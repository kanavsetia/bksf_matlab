function y=bj_general(edges,j)
x=[0  1 ; 1   0];
y=[0 -1j; 1j  0];
z=[1  0 ; 0  -1];
I=eye(2);
A=eye;
S=[];
for d1=1:length(edges)
    if edges(d1,1)==j|edges(d1,2)==j
        A=kron(z,A);
        edstr= strcat('Z',strcat(num2str(edges(d1,1)),num2str(edges(d1,2))));
        edstr=strcat(strcat('(',edstr),')');
        S=strcat(edstr,S);
    else        
        A=kron(I,A);
        edstr=strcat('I',strcat(num2str(edges(d1,1)),num2str(edges(d1,2))));
        edstr=strcat(strcat('(',edstr),')');
        S=strcat(edstr,S);
    end
end

y=S;

end

        
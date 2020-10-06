function y=Cdisp(h)

x=size(h)
C={};
for d1=1:x(1)
    for d2=1:x(2)
        if 1/1e10*round(1e10*real(h(d1,d2)))==0 & 1/1e10*round(1e10*imag(h(d1,d2)))==0
            C{d1,d2}=0;
        elseif 1/1e10*round(1e10*real(h(d1,d2)))==0
            C{d1,d2}=strcat(num2str(1/1e10*round(1e10*imag(h(d1,d2)))),'i');
        elseif 1/1e10*round(1e10*imag(h(d1,d2)))==0
            C{d1,d2}=1/1e10*round(1e10*real(h(d1,d2)));
        else
            C{d1,d2}=strcat(num2str(1/1e10*round(1e10*real(h(d1,d2)))),'+',num2str(1/1e10*round(1e10*imag(h(d1,d2)))),'i');
        end
    end
end
y=C;
end

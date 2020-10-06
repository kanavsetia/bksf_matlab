function y=animation(p)

np=size(p);

r=linspace(0,0.5,15);
X=[];
Y=[];
Z=[];
X1=[];
Y1=[];
Z1=[];
C=[];
h1=figure(2);
h1.Position=[1 1 2500 1000];
for d1=1:np(1)
    if d1==1
        for d2=1:np(2)
            [val ind]=min(abs(p(d1,d2)/2-r));
            for d3=1:ind
                points = round(r(d3)*15);
                if points ~=0
                    [X1,Y1,Z1] = sphere(points);
                    X=[X+d2-1;r(d3).*X1(:)];
                    Y=[Y;r(d3).*Y1(:)];
                    Z=[Z;r(d3).*Z1(:)];
                    %     keyboard
                    
                    C=[C;d3.*ones(length(X1(:)),1)];
                end
            end
            h=scatter3(X,Y,Z,60,'C','.');
            h.XDataSource='X';
            h.YDataSource='Y';
            h.ZDataSource='Z';
            h.CDataSource='C';
            
        
            
        end
    else
        for d2=1:np(2)
            [val ind]=min(abs(p(d1,d2)/2-r));
            for d3=1:ind
                points = round(r(d3)*15);
                if points ~=0
                    
                    [X1,Y1,Z1] = sphere(points);
                    X=[X+d2-1;r(d3).*X1(:)];
                    Y=[Y;r(d3).*Y1(:)];
                    Z=[Z;r(d3).*Z1(:)];
                    C=[C;d3.*ones(length(X1(:)),1)];
                else
                    
                end
            end
        
            
            
        end
        
        h.XDataSource='X';
        h.YDataSource='Y';
        h.ZDataSource='Z';
        h.CDataSource='C';
        
        refreshdata;
        drawnow;
        pause(0.1)
        
    end
    
end

% map=colormap;
% for d1=1:length(r)
%     points = round(r(d1)*25);
%     if points ~=0
%         [X1,Y1,Z1] = sphere(points);
%         X=[X;r(d1).*X1(:)];
%         Y=[Y;r(d1).*Y1(:)];
%         Z=[Z;r(d1).*Z1(:)];
%         %     keyboard
%
%         C=[C;d1.*ones(length(X1(:)),1)];
%     end
% end
% % cm=linspace(1,length(r),length(r));
% % C = repmat([1 2 3],numel(X),1);
%
% scatter3(X,Y,Z,60,C,'.')
%
% hold on
% % for d1=1:length(r)
% %     points = round(r(d1)*15);
% %     [X1,Y1,Z1] = sphere(points);
% %     X=[X;r(d1).*X1(:)];
%     Y=[Y;r(d1).*Y1(:)];
%     Z=[Z;r(d1).*Z1(:)];
%     C=[C;d1.*ones(length(X1(:)),1)];
% end
% X=X+3;
% scatter3(X,Y,Z,100,C,'.')
%
% for d1=1:length(r)
%     points = round(r(d1)*15);
%     [X1,Y1,Z1] = sphere(points);
%     X=[X;r(d1).*X1(:)];
%     Y=[Y;r(d1).*Y1(:)];
%     Z=[Z;r(d1).*Z1(:)];
%     C=[C;d1.*ones(length(X1(:)),1)];
% end
% X=X+3;
% scatter3(X,Y,Z,100,C,'.')
% for d1=1:length(r)
%     points = round(r(d1)*15);
%     [X1,Y1,Z1] = sphere(points);
%     X=[X;r(d1).*X1(:)];
%     Y=[Y;r(d1).*Y1(:)];
%     Z=[Z;r(d1).*Z1(:)];
%     C=[C;d1.*ones(length(X1(:)),1)];
% end
% X=X+3;
% scatter3(X,Y,Z,100,C,'.')
%
end

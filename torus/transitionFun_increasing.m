function y = transitionFun_increasing(x,L,flag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This routine provides a smooth increasing transition function from 0 to 1 over the interval
% x in [0,L]. Depending on the select flag the shape of the function will
% change:
%
% flag=1: standard mollifier
% flag=2: unique cubic polynomial p(x) on [0,L] with p(0)=p'(0)=p'(L)=0 and p(L)=1
% flag=3: transition function built out of the prototypical smooth but non
%         analytic function, see https://en.wikipedia.org/wiki/Non-analytic_smooth_function
% flag=4: linear function on [0,L]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Running example:

% x=linspace(-1,5,1000);
% L=3;
% figure(1)
% col=jet(4);
% for j=1:4
%     y=transitionFun_increasing(x,L,j);
%     plot(x,y,'LineWidth',2,'Color',col(j,:))
%     if j==1 hold on, end
%     grid on 
%     grid minor
% end
% hold off
% tit=title(['different transition functions over the length scale $L=',num2str(L),'$']);
% tit.Interpreter='latex';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    y=zeros(size(x));
    idx_center=(x>0 & x<L);
    idx_right=(x>=L); y(idx_right)=1;
    
    switch flag
        case 1
            y(idx_center)=exp(-1./(1-(x(idx_center)/L-1).^2)+1);    
        case 2
            y(idx_center)=(-2/L^3)*x(idx_center).^3+(3/L^2)*x(idx_center).^2;
            
        case 3
            gFun1=gFun(x(idx_center)/L);
            gFun2=gFun(1-x(idx_center)/L);
            y(idx_center)=gFun1./(gFun1+gFun2);
            
        case 4
            y(idx_center)=x(idx_center)/L;     
    end
end

function y2=gFun(x)
    y2=exp(-1./x);
end

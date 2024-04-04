function y = transitionFun(x,L,flag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This routine provides a continuous/smooth decaying transition function from 1 to 0 over the interval
% x in [0,L]. Depending on the selected flag the shape of the function will
% change:
%
% flag=1: shifted and scaled standard mollifier (smooth)
% flag=2: unique cubic polynomial on [0,L] (smooth)
% flag=3: transition function built out of the prototypical smooth but non
%         analytic function (smooth), see https://en.wikipedia.org/wiki/Non-analytic_smooth_function
% flag=4: linear function on [0,L] (continuous)
%
% Running example:
%
% x=linspace(-1,5,1000);
% L=3;
% figure(1)
% col=jet(4);
% for j=1:4
%     y=transitionFun(x,L,j);
%     plot(x,y,'LineWidth',2,'Color',col(j,:))
%     if j==1 hold on, end
%     grid on 
%     grid minor
% end
% hold off
% tit=title(['different transition functions over the length scale $L=',num2str(L),'$']);
% tit.Interpreter='latex';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    y=ones(size(x));
    idx_center=(x>0 & x<L);
    idx_right=(x>=L); y(idx_right)=0;
    
    switch flag
        case 1
            y(idx_center)=exp(-1./(1-(x(idx_center)/L).^2)+1);    
        case 2
             y(idx_center)=-2*(1-(x(idx_center)/L).^3)+3*(1-(x(idx_center)/L).^2);    
        case 3
            gFun1=gFun(1-x(idx_center)/L);
            gFun2=gFun(x(idx_center)/L);
            y(idx_center)=gFun1./(gFun1+gFun2);
            
        case 4
            y(idx_center)=1-x(idx_center)/L;     
    end
end

function y2=gFun(x)
    y2=exp(-1./x);
end

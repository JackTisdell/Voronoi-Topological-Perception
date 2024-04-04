%% Parameters for displaying and saving animation frames

iterNum=10; % iteration number to be displatyed, if you want all to be ploted set it to 'all'
% iterNum='all';
iterInterval=0; % this gives an intrval of iterations to be displayed, precisely (iterNum-iterInterval:iterNum+iterInterval)

output_type=0;   % To produce ``non-MatLab-looking'' figures for professional use. set to 1 for only the torus with no margins, labels, titles, etc. and pdf output. Set to 1 for png output with data plots (intended for video frames)
show_plots=0;

b_0 = 0;      % minimum brightness value for HSB colormap (image of |V_i| = infty)
A_0 = 3/4;      % Values of n|V_i|/|\Omega| less than this map to brightness 1

% set the four quantities below larger to make heavier strokes for smaller intended printing. 
sz =300;        % weight of position markers
qw=3;         % weight of velocity markers (from quiver)
ql=qw/90;       % length of velocity markers (proportional to weight)
ve=1;           % weight of Voronoi edges

save_frames=0; % 1 if you want to save, 0 for otherwise.
pic_count=0;
if iterNum~='all', pic_count=max(0,iterNum-iterInterval); end
currDir=pwd; %current directory, to save pictures.
path='/output';

if save_frames==1
  answerQ= questdlg('Do you want to save figure frames?','Watch Out','Yes','No','No');
  if strcmp(answerQ,'No')
    save_frames=0;
  end
end

disp_plot=1; %display plot frames of the animation (Yes/No)
disp_a_colormap=1; %display colored Voronoi cells (Yes/No)

%% Loads the data produced by Script_Dynamics.
load("data.mat")

%% Useful variables to declare
Omega=[[0;Lx;Lx;0],[0;0;Ly;Ly]]; % hall way domain
N=Nr+Nl; %total number of agents;

dt=1; %time scale
kMax=size(DT_t,1); %maximal number of iterations

ee=exp(1);
rE=N^2*18*sqrt(3)/(5*(Lx*Ly)^2); %normalization factor for Voronoi energy
sectorx_basis=[0;0;2*Lx;2*Lx]; % To help define \overline(Wi)
sectory_basis=[Lx;-Lx;-Lx;Lx];
if q==Inf, chara='\infty'; else, chara=num2str(q); end %for the title of the figure.
c_del=1.5; % delta-vecinity factor for points to be consider by ID_copy (for coloring purposes), i.e we will use the color-delta neighborhood of the domain Omega

%% Plot Animation
if iterNum=='all'
    iterDisplay=1:kMax;
else
    iterDisplay=max(1,iterNum-iterInterval):min(kMax,iterNum+iterInterval);
end

fig=figure(1); 
fig.Units='inches';  
fig.Position=[0 0 8+output_type 9];     % only relevant for non-pdf output. If different from PaperSize, pdf output will look different but will fill page regardless.
fig.PaperSize=[8+output_type 9];

clf
for k=iterDisplay
    
    %% Extract information from the saved data
    DT=DT_t{k};
    [V,C]=voronoiDiagram(DT);
    x4=DT.Points(:,1); y4=DT.Points(:,2);
    x=x4(1:N);y=y4(1:N);
    idx_copy=idx_copy_t{k};
    
    angle_vi=Angle_vi_t(k,:)';
    area=Area_t(k,:)';
    vx_old=Vx_old_t(k,:)';
    vy_old=Vy_old_t(k,:)';
    vx=Vx_t(k,:)';
    vy=Vy_t(k,:)';
    vmag=sqrt(vx.*vx + vy.*vy);
    
    %% Plot
    if disp_plot==1
        
        sep=.015;
        marg=.04;
        pos1=[.106 .8-marg .2 .2];
        pos2=[.38 .9+sep-marg .514 .1-sep];
        pos3=[.38 .8-marg .514 .1-sep];
        pos4=[0 .02 1 .7];
        
        wt=1.2;
        if output_type==0
            subplot('Position',pos1)
            polarscatter(angle_vi,vmag,'k.')
            thetaticks([]);
            rlim([0 2]);
            rticks([.5 1 1.5]);
    %         title('$v_i$','Interpreter','Latex','FontSize',14);
            title('velocities');

            subplot('Position',pos2)
            plot(0:k-1,Polarization_t(1:k),'r','LineWidth',wt)
            xlim([-5 kMax+5])
            ylim([0 1.05])
            xticklabels({});
            yticks([0 1]);
            ylabel('$P$','Interpreter','Latex');
            title('polarization');

            subplot('Position',pos3)
            plot(0:k-1,Energy_t(1:k),'b','LineWidth',wt)
            xlim([-5 kMax+5])
            ylim([.95 3])
            ylabel('$E$','Interpreter','Latex');
            title('clustering');



            subplot('Position',pos4)
        
        end
        %0) Plot the Domain
        scatter(0,0,'.w','MarkerEdgeColor','white') % DON'T REMOVE THIS LINE!!! (we use it to open the graphic handle and doesn't affect the graphs)
        hold on
        rectangle('Position',[0 0 Lx Ly])
        
        %idx_ghosts=(N+5:length(x4))';
        %scatter(x4(idx_ghosts),y4(idx_ghosts),'k.')
        
        % 2) Voronoi Diagram of the cylinder
        %     vor=voronoi(DT,'k-');
        %     vor(1).MarkerEdgeColor='none';
        %     vor(1).MarkerFaceColor='none';
        
        % 3) Colouring of the Voronoi regions based on the norm of \tilda{a} (without the scaling factor \phi_i)
        if disp_a_colormap==1
            
            % GrayScale: we delcare a grayscale in HSV channels and cap it
            % at some "Oxford Gray", i.e. we remove a part of the spectrum
            % closest to black (so that Voronoi regions are never fully
            % black...it looks very ugly).
            % Finally, we transform the grayscale to RBG format
            hsv_grayScale=zeros(64,3);
            hsv_grayScale(:,3)=linspace(0,1,size(hsv_grayScale,1));
            hsv_grayScale=hsv_grayScale(20:end,:);
            rgb_grayScale=hsv2rgb(hsv_grayScale); % we 
            
            % Computing the HSV pallette
            hsv_colors=ones(N,3);
%             hsv_colors(:,1)=(angle_vi+pi)/(2*pi); % transform the values angle_vi to [0 1] to use them as Hue values
            % Next we compute the Valor channel (third one). We make it in
            % a way that Voronoi regions with area equal to 2*|Omega|/N are
            % matched at 1 (i.e. full black). But then we cap that value
            % 2*|Omega|/N to the max value in the GrayScale of hsv_grayScale
%             hsv_colors(:,3)=max(ones(N,1)-hsv_grayScale(1,3),0.5*min(2*ones(N,1),(Lx*Ly)./(N*(area)))); %the saturation given by the Voronoi area scaled w.r.t |Omega|/N, however under this scaling the saturation can exceed 1...if that's the case we simply bring the saturation to 1.
            hsv_colors(:,3) = rescale(1./max(ones(N,1), 1/(A_0*Lx*Ly/N).*area), b_0, 1);
%             hsv_colors(:,3) = ones(N,1)-tanh(1/2*N/(Lx*Ly).*area);
            hsv_colors(:,2)=ones(N,1).*0;
            rgb_colors=hsv2rgb(hsv_colors);
            %NOTE, the second color of hsv_colors corresponding to "Saturation" is left to be identically one, i.e. all colors in our palette will have max saturation!
            
   
            
            % Coordinates of the vertices of the rectangle forming the c_del
            % neigborhood of D0
            d1=Omega(:,1)*c_del+0.5*Lx*(1-c_del);
            d2=Omega(:,2)*c_del+0.5*Ly*(1-c_del);
            
            % ID_Copy will contain the index of all original+ghost generators whose
            % Voronoi cell need to be colred so that there no Voro cell remains white
            % when raping around the boundaries of the torus.
            ID_copy=num2cell((1:N)');
            for j=1:length(idx_copy)
                jj=idx_copy(j);
                idx=N+4+j+(0:7)*length(idx_copy);
                in=inpolygon(x4(idx),y4(idx),d1,d2); % find points inside the c_del neighborhood
                ID_copy{jj}=[ID_copy{jj}  idx(in)]; %adding to the ID_copy list the indices of the points in the c_del neighborhood
            end
            
            % Color the Voronoi cell of each generator as well as the pertinent
            % copies.
            for i=1:N
                for j=1:length(ID_copy{i})
                    jj=ID_copy{i}(j);
                    patch(V(C{jj},1),V(C{jj},2),rgb_colors(i,:),'linewidth',ve,'EdgeAlpha',1)
                end
            end
            
            % Plot the colorbar
            if (false)
                c=colorbar;
                colormap(c,'hsv')
                c.Location='southoutside';
                c.Label.String='Angle of $v_i$';
                c.Label.Interpreter='latex';
                c.Label.FontSize=15;
                c.Ticks=[0 .5 1];
                c.TickLabels={'$0$' '$\pi$' '$2\pi$'};
                c.TickLabelInterpreter='latex';
                c.FontSize=15;


                c2=colorbar;
                colormap(c2,flipud(rgb_grayScale)) %flip upside down the grayscale so that small Voronoi cells are brighter and larger ones are darker.
                labels=c2.TickLabels;
                labels=linspace(0,2,size(labels,1))';
                c2.TickLabels=labels;
                c2.Location='eastoutside';
                c2.Label.String='$n|V_i|/|\Omega|$';
                c2.Label.Interpreter='latex';
                c2.Label.FontSize=15;
            end
            
        end
        

        % 4) Unitary velocity vector field (simply to showcast heading of agents)
        quiv=quiver(x(1:Nr),y(1:Nr),vx_old(1:Nr)*Ly*ql,vy_old(1:Nr)*Ly*ql,0);
        quiv.ShowArrowHead='off';
        quiv.AutoScale='off';
        quiv.Color=[1 .42 0];
        quiv.LineWidth=qw;
        
        quiv=quiver(x(Nr+1:end),y(Nr+1:end),vx_old(Nr+1:end)*Ly*ql,vy_old(Nr+1:end)*Ly*ql,0);
        quiv.ShowArrowHead='off';
        quiv.AutoScale='off';
        quiv.Color=[1 .42 0];
        quiv.LineWidth=qw;
        
        % 5) plot all agents
        scatter(x(1:Nr),y(1:Nr),sz,[1 .42 0],'.')
        scatter(x(Nr+1:N),y(Nr+1:N),sz,[1 .42 0],'.')
        
        % 6) Labels of agents
        %     txtlabels = arrayfun(@(n) {sprintf('%d', n)}, idx_agent);
        %     labels=text(x_plot,y_plot,txtlabels);
        %
        %     txtlabels = arrayfun(@(n) {sprintf('%d', n)}, idx_ghosts);
        %     labels=text(x4(idx_ghosts),y4(idx_ghosts),txtlabels,'Color','red');
        
        axis equal
        axis([0 Lx 0 Ly ])
        hold off
        xticks([])
        if (false)
            yticks([0 L])
            yticklabels({'' 'L'})
            title({['Torus with $n=',num2str(N),': L=',num2str(L),', p=',num2str(p),', q=',chara,', \nu=',num2str(nu),'$'];['$\#$iter $=',num2str(k-1),'$']},'Interpreter','latex')
        else
            ax = gca;
            ax.LooseInset=[0 0 0 0];    % zero all margins
            yticks([])
        end

        pause(0.000001)
        
        if save_frames==1
            %tit=['L=',num2str(L),'_iter=',num2str(k-1)];
            tit=num2str(pic_count);
            if pic_count<10 tit=strcat('000',tit);
            elseif pic_count<100 tit=strcat('00',tit);
            elseif pic_count<1000 tit=strcat('0',tit);
            end
%             print(tit,'-dpng')
            ext='.png';
            if output_type==0
                print(tit,'-dpng','-r200') % To regulate the resolution (it becomes heavier)
            else
                print(tit,'-dpdf','-r600','-fillpage') % in .pdf it looks amazing (vectorial
                ext='.pdf';
            % image) but it's very heavy
            end
                        
            ss1=strcat(tit,ext);
            ss2=strcat(currDir,path);
            movefile(ss1,ss2)
        end
        
        pic_count=pic_count+1;
        
    end
    disp(num2str(k))
    
end
%% Plot order parameters (metrics)

if (show_plots==1)
    fig=figure(2); fig.Units='centimeters';  fig.Position=[0 5.9972 28.3633 22.3661];
    clf
    %Polarization
    subplot(2,2,[1 2])
    plot(0:kMax-1,Polarization_t,'k')
    xlabel('$\#iter$','Interpreter','latex','FontSize',15)
    xlim([-5 kMax+5])
    ylim([-0.01 1.01])
    %ylabel('$$','Interpreter','latex','FontSize',15)
    title('Polarization over time ','Interpreter','latex');

    subplot(2,2,3)
    plot(0:kMax-1,Energy_t,'b')
    xlabel('$\#iter$','Interpreter','latex','FontSize',15)
    xlim([-5 kMax+5])
    %ylabel('$$','Interpreter','latex','FontSize',15)
    title('Scaled Voronoi Energy over time ','Interpreter','latex');


    subplot(2,2,4)
    plot(0:kMax-1,AvgVel_t,'r')
    hold on
    plot(0:kMax-1,AvgVel_t-StdVel_t,'r:')
    plot(0:kMax-1,AvgVel_t+StdVel_t,'r:')
    xlabel('$\#iter$','Interpreter','latex','FontSize',15)
    xlim([-5 kMax+5])
    %ylabel('$$','Interpreter','latex','FontSize',15)
    title('Average velocity with standard deviation  ','Interpreter','latex');


    suptitle({['Torus with $N=',num2str(N),': L=',num2str(L),', p=',num2str(p),', q=',chara,', \nu=',num2str(nu),'$']})
end

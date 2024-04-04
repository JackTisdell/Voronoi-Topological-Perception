%% Loads the data produced by Script_Dynamics.
load("data.mat")

%% Parameters for displaying and saving animation frames
iterNum=10; % iteration number to be displatyed, if you want all to be ploted set it to 'all'
iterInterval=0; % this gives an intrval of iterations to be displayed, precisely (iterNum-iterInterval:iterNum+iterInterval)
%iterNum='all';

disp_a_colormap=1; % display colored Voronoi cells (Yes) or display the Earth (No). The latter is less expensive

save_frames=0; % 1 if you want to save, 0 for otherwise.

tri_layer=5; % number of layers in the triangulation (meshing) used to smooth out 
             % the surface of the sphrical Voronoi polygons. 
             % the larger the number, the smoother the plotting.

angle_view=[-97.6 1.2]; %Point of view for the sphere [azimut elevation]

currDir=pwd; %current directory, to save pictures.
%path='/runs/Jack/color_10x10_N300_L1_nu1_p1_qInf';
%path='/movies/demo1';
path='/runs/test';

if save_frames==1
  answerQ= questdlg('Do you want to save figure frames?','Watch Out','Yes','No','No');
  if strcmp(answerQ,'No')
    save_frames=0;
  end
end

%% Useful variables to declare
N=Nr+Nl; %total number of agents;
kMax=size(DT_t,1); %maximal number of iterations
% kMax=801;
ee=exp(1);
if q==Inf, chara='\infty'; else, chara=num2str(q); end %for the title of the figure.

%% Preparing the graphic axis handle
pic_count=0;
if strcmp(iterNum,'all')
    iterDisplay=1:kMax;
else
    iterDisplay=max(1,iterNum-iterInterval):min(kMax,iterNum+iterInterval);
    pic_count=max(0,iterNum-iterInterval);
end

fig=figure(2); 
fig.Units='inches';  
fig.Position=[0 0 8 9];     % only relevant for non-pdf output. If different from PaperSize, pdf output will look different but will fill page regardless.
fig.PaperSize=[8 9];
clf(fig);

color1=[1 .42 0]; %orange 
color2=color1;
%color2=[0 .6 0]; % forest green

sep=.015; %separation between subplots
marg=.04; %margins separation
pos1=[.106 .8-marg .2 .2]; %position of subplots
pos2=[.38 .9+sep-marg .514 .1-sep];
pos3=[.38 .8-marg .514 .1-sep];
pos4=[0 .02 1 .7];
wt=1.2; % Line Width for plots of Order Parameters.

b_0 = 0;      % minimum brightness value for HSB colormap (image of |V_i| = infty)
A_0 = 3/4;      % Values of n|V_i|/|\Omega| less than this map to brightness 1

% set the four quantities below larger to make heavier strokes for smaller intended printing. 
sz =300;        % weight of position markers
qw=3;         % weight of velocity markers (from quiver)
ql=qw/90;       % length of velocity markers (proportional to weight)
ve=1.5;           % weight of Voronoi edges

flag1=0;
flag2=0;

addpath(genpath('PlaneEarth')) %add path to the PlotEarth folder have the Earth as a graphic primitive
%%
for k=iterDisplay
    %% Extract information from the saved data
    DT=DT_t{k}; V=V_t{k}; C=C_t{k};
    XYZ=XYZ_t(:,:,k);
    Vel=Vel_t(:,:,k); % velocity vectors;
    Vel_old=Vel_old_t(:,:,k); % unitary velocity vectors;
    VelMag=sqrt(sum(Vel.^2,2)); %magnitude of velocity vectors.
    
    AngleVel=AngleVel_t(k,:)'; %Angles of Velocity vectors with respect to an arbitrary direcion (All vectors where brought to the same plane to compute these angles)
    Area=Area_t(k,:)'; %area of each Voronoi cell

    %% Plot
    
    subplot('Position',pos1); %plotting endpoints of velocity vectors at iteration k in a polar plot
    if flag1==0
        
        polscat=polarscatter(AngleVel,VelMag,'k.');
        hold on
        polax=polscat.Parent;
        polax.ThetaDir='clockwise';
        polax.ThetaZeroLocation='top';
        polax.RTickMode='auto';
        polax.RTick=[.5 1];
        polax.RLim=[0 1.25];
        polax.ThetaTick=[0 45 90 135 180];
        polax.ThetaTickLabel={'$0$','$\pi/4$','$\pi/2$','$3\pi/4$','$\pi$'};
        polax.TickLabelInterpreter='latex';
        polax.ThetaLim=[0 180];
        hold off
        flag1=1;
    end
    
    polscat.ThetaData=AngleVel;
    polscat.RData=VelMag;
    drawnow limitrate % display updates
    title({'deviation from $\langle\omega\rangle_j$'},'Interpreter','Latex','FontSize',13);
    
    subplot('Position',pos2) % plotting the polarization up to iteration k
    plot(0:k-1,Polarization_t(1:k),'r','LineWidth',wt)
    xlim([-5 kMax+5])
    ylim([0 1.05])
    xticklabels({});
    yticks([0 1]);
    ylabel('$P$','Interpreter','Latex');
    title('polarization','Interpreter','Latex','FontSize',13);
    
    subplot('Position',pos3) % plotting the Voronoi Energy (clustering measure) up to iteration k
    plot(0:k-1,Energy_t(1:k),'b','LineWidth',wt)
    xlim([-5 kMax+5])
    ylim([.95 max(Energy_t)])
    ylabel('$E$','Interpreter','Latex');
    title('clustering','Interpreter','Latex','FontSize',13);
    
    sbplt4=subplot('Position',pos4); % Plot the configuration at iteration k

    if disp_a_colormap==1

        % Computing the GrayScale, first in HSV channels and then to RGB
        hsv_colors=ones(N,3); hsv_colors(:,2)=zeros(N,1);
        hsv_colors(:,3) = rescale(1./max(ones(N,1), N*Area./(A_0*4*pi)), b_0, 1);
        % hsv_colors(:,3) = ones(N,1)-tanh(1/2*N/(Lx*Ly).*area);
        rgb_colors=hsv2rgb(hsv_colors);

        % Agents and unit length orientation
        quiv=quiver3(XYZ(1:Nr,1),XYZ(1:Nr,2),XYZ(1:Nr,3),Vel_old(1:Nr,1),Vel_old(1:Nr,2),Vel_old(1:Nr,3),0.15);
        quiv.ShowArrowHead='off'; quiv.AutoScale='on'; quiv.Color=color1; quiv.LineWidth=1.8;
        hold on
        
        quiv=quiver3(XYZ(Nr+1:end,1),XYZ(Nr+1:end,2),XYZ(Nr+1:end,3),Vel_old(Nr+1:end,1),Vel_old(Nr+1:end,2),Vel_old(Nr+1:end,3),0.15);
        quiv.ShowArrowHead='off'; quiv.AutoScale='on'; quiv.Color=color2; quiv.LineWidth=1.8;
     
        % Recursive triangulation of Voronoi Polygons to make the spherical
        % polygons appear smoother and well round.
        inds_tri=(1:N)'; vertices_tri=V; tri=C;
        for tl=1:tri_layer
            [vertices_tri,tri,inds_tri]=triangulateVoronoiSphere(vertices_tri,tri,inds_tri);
        end
        COL=rgb_colors(inds_tri,:); % attribute the color of the original face to each of the final triangles.
        patch('Vertices',vertices_tri,'Faces', tri, 'EdgeColor','none','FaceVertexCData', COL,'FaceColor','flat')
        
        % Discretize the spherical edges so that they do not appear as
        % straight lines. We also plot the Voronoi cells a second time...it
        % gives a smoother plot.
        [Vx_discrete,Vy_discrete,Vz_discrete]=VoronoiEdges_discretized_sphere(V,C,2*pi/180);
        for i=1:N
            patch(Vx_discrete{i},Vy_discrete{i},Vz_discrete{i},rgb_colors(i,:),'LineWidth',ve)
        end

        scatter3(XYZ(1:Nr,1),XYZ(1:Nr,2),XYZ(1:Nr,3),10,'MarkerFaceColor',color1,'MarkerEdgeColor',color1) %objectives
        scatter3(XYZ(Nr+1:end,1),XYZ(Nr+1:end,2),XYZ(Nr+1:end,3),10,'MarkerFaceColor',color2,'MarkerEdgeColor',color2) %objectives
        
        sbplt4.View=angle_view;
        hold off
        axis equal
        axis manual
        axis off
        pause(0.000001)
        
    else
        if flag2==0
            hold on
            Earth=plotearth('MapType','bluemarble');
            quiv1=quiver3(XYZ(1:Nr,1),XYZ(1:Nr,2),XYZ(1:Nr,3),Vel_old(1:Nr,1),Vel_old(1:Nr,2),Vel_old(1:Nr,3),0.15);
            quiv1.ShowArrowHead='off'; quiv1.AutoScale='on'; quiv1.Color=color1; quiv1.LineWidth=1.8;
            quiv2=quiver3(XYZ(Nr+1:end,1),XYZ(Nr+1:end,2),XYZ(Nr+1:end,3),Vel_old(Nr+1:end,1),Vel_old(Nr+1:end,2),Vel_old(Nr+1:end,3),0.15);
            quiv2.ShowArrowHead='off'; quiv2.AutoScale='on'; quiv2.Color=color2; quiv2.LineWidth=1.8;
            sbplt4.View=angle_view;%[-216.3 2];
            hold off
            axis equal
            axis manual
            axis off
            flag2=1;
        end

        quiv1.XData=XYZ(1:Nr,1); quiv1.YData=XYZ(1:Nr,2); quiv1.ZData=XYZ(1:Nr,3); %update quiver's arrow starting point
        quiv1.UData=Vel_old(1:Nr,1); quiv1.VData=Vel_old(1:Nr,2); quiv1.WData=Vel_old(1:Nr,3); %update quiver's arrows directions
        quiv2.XData=XYZ(Nr+1:end,1); quiv2.YData=XYZ(Nr+1:end,2); quiv2.ZData=XYZ(Nr+1:end,3); %update quiver's arrow starting point
        quiv2.UData=Vel_old(Nr+1:end,1); quiv2.VData=Vel_old(Nr+1:end,2); quiv2.WData=Vel_old(Nr+1:end,3); %update quiver's arrows directions
        drawnow limitrate % display updates
    end

    if save_frames==1
        tit=num2str(pic_count); tit=strcat('0000',tit); tit=tit(end-3:end); % up to iteration 9,999.
        ext='.png';
        print(tit,'-dpng','-r150') % To regulate the resolution (it becomes heavier)
        ss1=strcat(tit,ext);
        ss2=strcat(currDir,path);
        movefile(ss1,ss2)
    end
    
    pic_count=pic_count+1;
    
    disp(num2str(k))
    
end
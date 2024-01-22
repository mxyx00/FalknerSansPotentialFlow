function [separation] = panelCode(airfoil,aoa,velocity, resolution)

%Inputs: Airfoil Data, Angle of Attack (degrees), Flow velocity (m/s),
%Resolution (n)
%Outputs: Potential Flow Vector Figure

%% Errors
% Angle of Attack and Flow Velocity are optional parameters. If not
% specified then default is:
%   aoa = 0º
%   flow velocity = 2m/s
%   resolution = 50 pts.

if nargin < 2
    velocity = 2;
    aoa = 0;
    fprintf("Default aoa of 0º.\n")
end

if nargin <3
    velocity = 2;
    fprintf("Default velocity of 2m/s used.\n")
end

if nargin <4
    resolution = 50;
    fprintf("Default resolution used.\n")
end

% Setup
% Opening an airfoil data file. Reading bonyname and creating vector.

fid = fopen(airfoil);
xycell = textscan(fid, '%f %f','headerlines', 1); % Skipping titles and text
xy=cell2mat(xycell); % 
fclose(fid);
bodyname=textread(airfoil,'%s',1,'delimiter','\n');


% Preliminary Calculations
% Converting aoa from degrees to radians. Removing double-used points such
% as LE and TE, re-sorting direction of vector.

U = velocity; % free stream velocity
aoad = aoa
aoa=aoa*pi/180; %Converting angle of attack to radians. 

% Separating coordinates into x and y columns.
x=xy(:,1);
y=xy(:,2);


% Delete last data point (Redundant leading edge data point)
x(end)=[];
y(end)=[];

 
LE=find(x==min(abs(x)));

x=[flipud(x(1:LE)); flipud(x(LE:end))]*100;
y=[flipud(y(1:LE)); flipud(y(LE:end))]*100;

x(end)=[];
y(end)=[];

xy = [x,y];

% Parameterization
% Using spline function and cummulative distance in order to create an
% airfoil with collocation points incrementally further apart from each
% other.

% Splitting the airfoil points so they are symmetrical.
dist = vecnorm(diff(xy, 1, 1), 2, 2);

% Parameterize x and y coordinates
abscissa = [0; cumsum(dist)];

% Compute spline coefficients

pp = spline(abscissa, flip(xy',2));
suction = flip(coord(0, abscissa(LE),resolution)); % Spacing
pressure = coord(abscissa(end),abscissa(LE),resolution); % Spacing 
abscissa = [flip(suction), flip(pressure(1:end-1))];
fitSpline = flip(ppval(pp,abscissa),2); % coord 

xPlot = fitSpline(1,:)';
yPlot = fitSpline(2,:)';

xPlot = x;
yPlot = y;

n = length(xPlot);
t=find(x==max(x));
c=x(t);

xL=[xPlot(n) ; xPlot(1:(n-1))];	% nx1
yL=[yPlot(n) ; yPlot(1:(n-1))];	% nx1

xc=(xL+xPlot)./2;		% nx1
yc=(yL+yPlot)./2;

% Panel Plot
% Creates a plot with the airfoil and it's panels. Collocation points at
% the center of each panel. 
subplot(2,1,1)
AF_F1=zeros(n,1);
hold on
box on
count_even=0;
count_odd=0;

for i=1:n
    if i==n
        ip1=1;
    else
        ip1=i+1;
    end
    
    if mod(i,2)==0
        AF_F1(i)=plot([xPlot(i) xPlot(ip1)],[yPlot(i), yPlot(ip1)],'Color','b','LineWidth',2);
        count_even=count_even+1;
    else
        AF_F1(i)=plot([xPlot(i) xPlot(ip1)],[yPlot(i), yPlot(ip1)],'Color',[1 1 1],'LineWidth',4);
        count_odd=count_odd+1;
    end
end


AF_F2=plot(xc,yc,'rO','LineWidth',2,'MarkerEdgeColor','r','MarkerFaceColor','g','MarkerSize',5);
set(gca,'color','#d3d3d3');


% Plot Formatting and Etc. 

leg3=['Collocation Points ( ' num2str(length(xc)) ' )'];
legend([AF_F1(1),AF_F1(2),AF_F2],'Panels','Panels',leg3,'Location','NorthEast');

 
% RESAMPLE

x = xPlot(1:end-1);
y = yPlot(1:end-1);

xy = [x,y];

n = length(x);
i = n;
% Find index of trailing edge
t=find(x'==max(x'));
% Chord Length
c=x(t);

xL=[x(n) ; x(1:(n-1))];	% nx1
yL=[y(n) ; y(1:(n-1))];	% nx1

xc=(xL+x)./2;		% nx1
yc=(yL+y)./2;		% nx1

s=((x-xL).^2+(y-yL).^2).^0.5;	% nx1
tx=(x-xL)./s;	% nx1
ty=(y-yL)./s;	% nx1

nx=-ty;	% nx1
ny=tx;	% nx1
  
subplot(2,1,2)
    xT=x.*cos(aoa)+y.*sin(aoa);
    yT=-x.*sin(aoa)+y.*cos(aoa);

    maxT=max(xT);
    
    left=-round(2.5*n);
    right=100+round(2.5*n);
    bottom=min(yT)-round(1.2*n);
    top=max(yT)+round(1.2*n);
    
    Outer=[left right+4 bottom top];

    if max(yT)>40 | min(yT)<-40
        Focussed=[-80  150  min(yT)-40   max(yT)+30];
    else
        Focussed=[-32.5679  136.5679  min(yT)-40   max(yT)+30];
    end
    hold on
    box on
    axis equal
    
    VP1(1)=title(['Vector Plot of ' bodyname]);
    VP1(2)=xlabel('X \rightarrow');
    P1(3)=ylabel('Y \rightarrow');

    VP1(4)=plot([left right+4 right+4 left ,left], [bottom bottom  top top ,bottom],'r','LineWidth',3);

    axis(Outer);
    zoom reset

    Ntheta_i=atan2(ny,nx);
    Ntheta_j=[Ntheta_i(2:n);Ntheta_i(1)];

    Ntheta=(Ntheta_i+Ntheta_j)./2;
    Ntheta(1)=pi+Ntheta(1);

    xx=linspace(left,right,n);
    yy=linspace(bottom,top,n);

    xxc=repmat(xx,n,1);
    yyc=repmat(yy',1,n);

    xxcT=xxc*cos(-aoa)+yyc*sin(-aoa);
    yycT=-xxc*sin(-aoa)+yyc*cos(-aoa);

    xxc=xxcT;
    yyc=yycT;
    
    xlowT=[xT(find(xT==max(xT)):n);xT(1)]; % 26x1
    xupT=xT(1:find(xT==max(xT))); % 26x1

    ylowT=[yT(find(xT==max(xT)):n);yT(1)]; % 26x1
    yupT=yT(1:find(xT==max(xT))); % 26x1
    
    dxij_j=repmat(xc,1,n)-repmat(xc',n,1);	% nxn
    dyij_j=repmat(yc,1,n)-repmat(yc',n,1);  % nxn

    xij_j=dxij_j.*(repmat(tx',n,1))+dyij_j.*(repmat(ty',n,1));	% nxn
    yij_j=dxij_j.*(repmat(nx',n,1))+dyij_j.*(repmat(ny',n,1));	% nxn

    ss=repmat(s',n,1);	% nxn

    vxij_j=0.5.*(log((xij_j+0.5.*ss).^2+yij_j.^2)-log((xij_j-0.5.*ss).^2+yij_j.^2));    % nxn
    vyij_j=atan((xij_j+0.5.*ss)./yij_j)-atan((xij_j-0.5.*ss)./yij_j);   % nxn

    nitj=(nx*tx')+(ny*ty'); % nxn
    ninj=(nx*nx')+(ny*ny'); % nxn
    titj=(tx*tx')+(ty*ty'); % nxn
    tinj=(tx*nx')+(ty*ny'); % nxn

    sNij=vxij_j.*nitj+vyij_j.*ninj; % nxn
    sTij=vxij_j.*titj+vyij_j.*tinj; % nxn

    vNij=vyij_j.*nitj-vxij_j.*ninj; % nxn
    vTij=vyij_j.*titj-vxij_j.*tinj; % nxn

    sNij=sNij-diag(diag(sNij))+diag(repmat(pi,1,n));
    vTij=vTij-diag(diag(vTij))+diag(repmat(pi,1,n));

    vN=sum(vNij')'; % nx1
    vT=sum(vTij')'; % nx1

    N=[sNij vN];    % nx(n+1)
    T=[sTij vT];   % nx(n+1)

    Mnp1_row= T(t,:)+T((t+1),:); % 1x(n+1)

    M=[N;Mnp1_row];   % (n+1)x(n+1)

    Ux=U.*cos(aoa);   % 1x1
    Uy=U.*sin(aoa);   % 1x1

    Bi=-(Ux.*nx+Uy.*ny);    % nx1
    Bnp1=-(Ux.*(tx(t)+tx(t+1))+Uy.*(ty(t)+ty(t+1)));    % 1x1

    B=[Bi;Bnp1];    % (n+1)x1

    A=M^-1*B;   % (n+1)x1

    vni=N*A+(Ux.*nx+Uy.*ny);    % nx1 (=0)
    vti=T*A+(Ux.*tx+Uy.*ty);    % nx1 (!=0)

    TotalGamma=2*pi*A(n+1)*sum(s);
    L=1.225*U*TotalGamma; 

    Cpi=1-(vti./U).^2;  % nx1

    cL=L/(0.5*1.225*(U^2)*c);

    theta_i=atan(ty./tx);

    th1=[theta_i(n);theta_i(1:(n-1))];
    th2=[theta_i(2:n);theta_i(1)];

    Psi_i=abs((th1-th2)./2);

    R_i=s./(2.*sin(Psi_i./2));

    xxcT=xxc*cos(aoa)+yyc*sin(aoa);
    yycT=-xxc*sin(aoa)+yyc*cos(aoa);

    VP2=zeros(1,n);
    Vaoa=aoa;

    for j=1:n        
        dxij_jH=repmat(xxc(:,j),1,n)-repmat(xc',n,1);	% nxn
        dyij_jH=repmat(yyc(:,j),1,n)-repmat(yc',n,1);  % nxn

        xij_jH=dxij_jH.*(repmat(tx',n,1))+dyij_jH.*(repmat(ty',n,1));	% nxn
        yij_jH=dxij_jH.*(repmat(nx',n,1))+dyij_jH.*(repmat(ny',n,1));	% nxn


        vxij_jH=0.5.*(log((xij_jH+0.5.*ss).^2+yij_jH.^2)-log((xij_jH-0.5.*ss).^2+yij_jH.^2));    % nxn
        vyij_jH=atan((xij_jH+0.5.*ss)./yij_jH)-atan((xij_jH-0.5.*ss)./yij_jH);   % nxn


        sNijH=vxij_jH.*nitj+vyij_jH.*ninj; % nxn
        sTijH=vxij_jH.*titj+vyij_jH.*tinj; % nxn

        vNijH=vyij_jH.*nitj-vxij_jH.*ninj; % nxn
        vTijH=vyij_jH.*titj-vxij_jH.*tinj; % nxn


        sTijH=sTijH-diag(diag(sTijH));
        vNijH=vNijH-diag(diag(vNijH));

        vNH=sum(vNijH')'; % nx1
        vTH=sum(vTijH')'; % nx1

        NH=[sNijH vNH];    % nx(n+1)
        TH=[sTijH vTH];    % nx(n+1)

        vniH=NH*A+Ux.*nx+Uy.*ny; % nx1
        vtiH=TH*A+Ux.*tx+Uy.*ty; % nx1

        vx(:,j)=vniH.*nx+vtiH.*tx;%+Ux;
        vy(:,j)=vniH.*ny+vtiH.*ty;%+Uy;
     
            if max(yT)>40 | min(yT)<-40
                tolerance_over_surface=4;
            else
                tolerance_over_surface=2;
            end

            
            
       %------------------------------ Emptying the velocities inside the airfoil
        if xxcT(2,j)>=min(xT) & xxcT(2,j)<=max(xT) % coloums have same values

            lowdiff=abs( round(xlowT*100)-round(xxcT(2,j)*100) );
            updiff=abs( round(xupT*100)-round(xxcT(2,j)*100) )  ;
            
            icl=find(lowdiff==min(lowdiff));
            icu=find(updiff==min(updiff));
            
             
            if xxcT(i,j)<0

                for i=1:n % length(yy)

                    if yycT(i,j)-ylowT(icl)<=-tolerance_over_surface & yycT(i,j)-yupT(icu)>=tolerance_over_surface
                        vx(i,j)=0;
                        vy(i,j)=0;
                    end
                end

            else
        
                for i=1:n % length(yy)

                    if yycT(i,j)-ylowT(icl)>=-tolerance_over_surface & yycT(i,j)-yupT(icu)<=tolerance_over_surface
                        vx(i,j)=0;
                        vy(i,j)=0;
                    end
                end
                
            end

       end
%         ------------------------------ Emptying the velocities inside the airfoil

        vxT(:,j)=vx(:,j)*cos(aoa)+vy(:,j)*sin(aoa);
        vyT(:,j)=-vx(:,j)*sin(aoa)+vy(:,j)*cos(aoa);

        if xxcT(2,j)>min(xT) && xxcT(2,j)<max(xT)
            VP2(j)=quiver(xxcT(:,j),yycT(:,j),vxT(:,j),vyT(:,j),14/n,'b'); % vertical line vectors
        else
            VP2(j)=quiver(xxcT(:,j),yycT(:,j),vxT(:,j),vyT(:,j),8/n,'b'); % vertical line vectors
        end
            
    end
    
    VP4=quiver(xxcT,yycT,vxT,vyT,0.8,'b'); %

    
    VP3(:,:)=fill([xT;xT(1)],[yT;yT(1)],[0.5020         0    0.2510]);

    axis(Focussed)
    zoom reset




% Pressure Gradient Plot

figure()
s=((x-xL).^2+(y-yL).^2).^0.5;
psi1 = gradient(Cpi(2:t,1),s(2:t));
psi2 = gradient([ Cpi(t:end,1); Cpi(1,1)],[s(t:end);s(1)]);
plot(s(2:t),psi1,'r','LineWidth',1.5)
hold on
plot([s(t:end);s(1)],psi2,'b','LineWidth',1.5)
xlabel('X \rightarrow');
title(['Pressure gradient with respect to distance of ' bodyname]);
legend('Upper Surface', 'Lower Surface','Location','NorthEast' )
grid on
pbaspect([5 2 1])

%% Pressure through Source Vortex Method

% Resetting variables 
fid = fopen(airfoil);
xycell = textscan(fid, '%f %f','headerlines', 1); % Skipping titles and text
airfoil=cell2mat(xycell); % 
fclose(fid);
numPanels = 131;
alpha = aoad
numXGrid = 140;
numYGrid = 140;

dist = vecnorm(diff(airfoil, 1, 1), 2, 2);
abscissa = [0; cumsum(dist)];
pp = spline(abscissa, flip(airfoil', 2));
abscissaNew = linspace(abscissa(1), abscissa(end), numPanels + 1);
airfoilRefined = ppval(pp, abscissaNew);
airfoilRefined = airfoilRefined'; 

XB = airfoilRefined(:,1);
YB = airfoilRefined(:,2);
edge = zeros(numPanels,1);                                                     

for index = 1:numPanels                                                          
    edge(index) = ( XB(index+1)-XB(index) ) * ( YB(index+1)+YB(index) ); 
end

sumEdge = sum(edge);                                                        

if sumEdge < 0                                                          
    XB = flipud(XB);                                                        
    YB = flipud(YB);                                                       
end

XC   = zeros(numPanels,1);                                       
YC   = zeros(numPanels,1);                                      
SVec    = zeros(numPanels,1);                                              
phiDeg = zeros(numPanels,1);                                                

for j = 1:numPanels                                                          
    XC(j)   = XB(j) + 0.5 * (XB(j+1) - XB(j));                                         
    YC(j)   = YB(j) + 0.5 * (YB(j+1) - YB(j)); 
    deltaX      = XB(j+1)-XB(j);                                                
    deltaY      = YB(j+1)-YB(j); 
    SVec(j)    = (deltaX^2 + deltaY^2)^0.5;                                          
    phiDeg(j) = atan2d(deltaY, deltaX); 
    if (phiDeg(j) < 0)                                                           
        phiDeg(j) = phiDeg(j) + 360; 
    end
end

deltaDeg = phiDeg + 90; 
betaDeg = deltaDeg - alpha;                                        

betaDeg(betaDeg > 360) = betaDeg(betaDeg > 360) - 360;                              
phi  = phiDeg.*(pi/180); 
delta = deltaDeg * (pi/180);
beta = betaDeg.*(pi/180);                                                 

Vinf = U;

% Find the normal (I) and tangential (J) geometric integrals for the source
% panel method. 
[I,J] = computeIJ_SPM(XC,YC,XB,YB,phi,SVec);  

% Find the normal (K) and tangential (L) geometric integrals for the vortex
% panel method. 
[K,L] = computeKL_VPM(XC,YC,XB,YB,phi,SVec);              

% Construct the primary A matrix. 
A = I + pi * eye(numPanels,numPanels);

% Populate the right column of the A matrix. 
for j = 1:numPanels                                                  
    A(j,numPanels+1) = -sum(K(j,:));      
end

% Enforce the Kutta condition in the bottom row of the A matrix.
for j = 1:numPanels 
    
    % Account for the contribution of the sources to the Kutta condition.
   
    A(numPanels+1,j) = J(1,j) + J(numPanels,j) ;                                

end

%Account for the contribution of the vortex sheets to the Kutta condition.

A(numPanels+1,numPanels+1) = -sum(L(1,:) + L(numPanels,:)) + 2*pi;                  

% Define the b vector 
b = -Vinf * 2 * pi * cos(beta);

% Modify the last element of the b array to satisfy the Kutta condition.
b(numPanels+1) = -2 * pi * Vinf *( sin(beta(1)) + sin(beta(numPanels)) );          

% Compute the result array that contains the vortex strength and
% the source strengths for each of the panels. 
resultArray = A\b;                                                               

% Extract source strengths and vortex strength from the result array.
lambda = resultArray(1:end-1);                                          
gamma  = resultArray(end);                         

%% Determine the velocity and pressure coefficient along each panel. 

%Initialize tangential velocity array and pressure coefficient array.
Vtan = zeros(numPanels,1);                                               
CpAirfoil = zeros(numPanels,1);                                                  

for j = 1:numPanels
    %Determine tangential velocity and pressure coefficient for panel j.

    %Contribution from uniform flow. 
    term1 = Vinf*sin(beta(j));

    %Contribution from the source sheets on the other panels.
    term2 = (1/(2*pi))*sum(lambda .* J(j,:)');                              
    
    %Contribution from vortex sheet on panel j.
    term3 = gamma/2;   

    %Contribution from vortex sheets on the other panels. 
    term4 = -(gamma/(2*pi)) * sum(L(j,:));                                   
    
    %Define tangential velocity on panel j.
    Vtan(j) = term1 + term2 + term3 + term4;                                 
    
    %Define pressure coefficient on panel j. 
    CpAirfoil(j) = 1 - ( Vtan(j)/Vinf )^2;                                              
end

%% Compute the lift and moment coefficients.

% Determine the normal force coefficients on all the panels. Recall that
% the chord line is parallel to the x-axis and that delta is the angle 
% between a panel normal vector and the +x-axis.
CN = -CpAirfoil .* SVec .* sin(delta);                                                    

% Determine the axial force coefficients on all the panels. 
CA = -CpAirfoil .* SVec .* cos(delta);                                                     

% Compute lift, drag, and moment coefficients. 
Cl = sum(CN .* cosd(alpha)) - sum(CA .* sind(alpha));  

Cd = sum(-CN .* sind(alpha)) + sum(CA .* cosd(alpha));  

Cm = sum(CpAirfoil .* (XC - 0.25) .*SVec .* cos(phi));    


%% Calculate the streamlines away from the airfoil. 

% Define grid parameters.                                                       
xVals  = [min(XB)-0.5, max(XB)+0.5]; % x range                                
yVals  = [min(YB)-0.4, max(YB)+0.4]; % y range                          
    

% Define streamline parameters:

% Step size for propagation of streamlines. 
stepsize = 0.01;                                                       

% Maximum number of vertices. 
maxVert  = numXGrid*numYGrid*100; 

% Percentage of streamlines to be plotted on the grid.  
streamlinePercentage = 25;                                                       

% Create an array of streamline starting points parallel to the y-axis. 
yStreamLineStarts = linspace(yVals(1), yVals(2), ...
    floor( (streamlinePercentage/100)*numYGrid) )';     
    
% Generate the evenly-spaced grid points.
Xgrid   = linspace(xVals(1), xVals(2), numXGrid)'; 

Ygrid   = linspace(yVals(1), yVals(2), numYGrid)';                         

[xGrid, yGrid] = meshgrid(Xgrid, Ygrid);                                       
    
% Initialize matrices containing the normalized velocities. 
uMat = zeros(numXGrid, numYGrid);                                              
vMat = zeros(numXGrid, numYGrid);                                          
    
% Solve for grid point x and y velocities.
for m = 1:1:numXGrid
    for n = 1:1:numYGrid

        % Extract current grid point location. 
        XP      = xGrid(m,n);                                             
        YP      = yGrid(m,n);  

        % Compute Source Panel Method streamline geometric integrals. 
        [Mx,My] = streamline_SPM(XP,YP,XB,YB,phi,SVec); 

        % Compute Vortex Panel Method streamline geometric integrals. 
        [Nx,Ny] = streamline_VPM(XP,YP,XB,YB,phi,SVec);    


        % If the current grid point is on the airfoil, then set the velocity
        % equal to zero here. 
        [in,on] = inpolygon(XP,YP,XB,YB);

        if (in == 1 || on == 1)                                       
            uMat(m,n) = 0;                                               
            vMat(m,n) = 0; 
        
        % If the grid point is off the airfoil, then compute the true x-
        % and y-components of velocity.
       
        else                                                           
            uMat(m,n) = Vinf * cosd(alpha) + sum(lambda.*Mx./ (2*pi)) + ...    
                            sum(-gamma .*Nx ./ (2*pi));
                
            vMat(m,n) = Vinf * sind(alpha) + sum(lambda.*My./ (2*pi)) + ...    
                            sum(-gamma .* Ny ./ (2*pi));
        end
    end
end
    
% Determine the magnitude of velocity at the grid point and use that
% magnitude to compute the pressure coefficient. 
VMag  = sqrt(uMat.^2 + vMat.^2);                                      

CpMat = 1-(VMag ./ Vinf).^2;                                       

    % Now, plot the Cp distribution around the airfoil. 

    figure()
    cla
    hold on
    grid on

    set(gcf,'Color','White')                                               
    set(gca,'FontSize',15)   

    
    % Determine index of airfoil LE, and call it "midIndS," because it will
    % most likely be the index of the middle value in an array that
    % contains the cumulative arclength around the airfoil surface. 
    minXControlPoint = min(XC);

    midIndS = find(XC == minXControlPoint);

    % Plot Cp vs. x for the suction surface. 
    CpUpper = plot(XC(midIndS+1:end), CpAirfoil(midIndS+1 : end),'b-',...
        'LineWidth', 1.5);
    
    % Plot Cp vs. x for the pressure surface. 
    CpLower = plot(XC(1:midIndS),CpAirfoil(1:midIndS),'r-','LineWidth', 1.5);


    % Create legend and format the rest of the plot.
    legend([CpUpper,CpLower], {'Suction Surface', 'Pressure Surface'}, ...
        'FontSize', 12);
    
    title('$C_{p}$ vs. $x/c$', 'Interpreter','latex', 'FontSize', 20)
    
    xlabel('$x/c$', 'Interpreter','latex', 'FontSize', 20);                                              
    ylabel('$C_{p}$', 'Interpreter','latex', 'FontSize', 20);                                                          
    xlim(['auto']);                                                           
    ylim('auto');                                                          
    set(gca,'Ydir','reverse')  

    zoom reset; 

%% CALCULATING M + BETA
% Beta = angle
figure()
% for leading edge 
len = length(airfoil)-2;
beta = zeros(1,len);
m = zeros(1,len);

% A=[x(2),y(2)]-[x(1),y(1)];
% B=[x(end),y(end)]-[x(1),y(1)];
% 
% beta(1)=acos(sum(A.*B)/(norm(A)*norm(B)));
% 
% m(1) = beta/((2*pi)-beta);
% 
% % every other panel 
% for i = 3:(len)
%     if y(i) > y(i-1)
%         A=[x(i),y(i)]-[x(i-1),y(i-1)];
%         B=[100,0]-[0,0];
%         beta(i-1)=acos(sum(A.*B)/(norm(A)*norm(B)));
%     elseif y(i) < y(i-1)
%         A=[x(i),y(i)]-[x(i-1),y(i-1)];
%         B=[100,0]-[0,0];
%         beta(i-1)=-acos(sum(A.*B)/(norm(A)*norm(B)));
%     end
%     m(i-1) = beta(i-1)./((2*pi)-beta(i-1));
% end
% 
% %trailing edge index
te = find(x==100)

% 
% 
% line1 = plot(x(1:te),beta(1:te));
% hold on
% line2 = plot(x(1:te),m(1:te));
% 



% Beta as function of deriv. pressure/distance(s)

beta = zeros(1,te);
m = zeros(1,te);

CpUpper = CpAirfoil(midIndS+1 : end);
s2 = [0,0]
for i = 2:te
    distance = (sqrt((y(i)-y(i-1))^2+(x(i)-x(i-1))^2));
    s2(i) = s2(i-1) + distance;
end



dCpdx = gradient(CpUpper(1:te-1),x(1:te-1));
dCpds = gradient(CpUpper(1:te-1),s2(1:te-1));


for i = 1:te-1
     m(i) = -(0.5*x(i)*dCpdx(i));
     beta(i) = 2*m(i)/(m(i)+1);
end

for i = 1:te-1
     m2(i) = -(0.5*s2(i)*dCpds(i));
     beta2(i) = 2*m2(i)/(m2(i)+1);
end


for i = 1:te-1
    if beta(i) > -0.199
        sepIndx = beta(i)
    elseif beta(i) < -0.199
        break
    end
end

pt = find(beta==sepIndx)
separation = x(pt)

plot(x(1:te-1),beta(1:te-1));
hold on
plot(x(1:te-1),m(1:te-1));
yline(-0.199)
separationOccurs = (xline(separation,'LineWidth',1));
separationOccurs.LineStyle = '--';
dim = [.38 .4 .3 .3];
str = ['Separation Occurs at x=' int2str(separation)];
annotation('textbox',dim,'String',str,'FitBoxToText','on');
t = ["Beta & M vs. X (Upper side), Angle =" aoad]
title(t)
legend('Beta', 'M (dCp/dx)')

%% Ext Functions

function point = coord(a,b,N)
    point = a + (b-a)*fStretch(N);
end

function sy = fStretch(my)
    yStretch = 0.5;
    yOffset = 2.95;
    iy = linspace(0,1,my);
    sy = scale01(exp(yStretch*(yOffset+iy)).^2);
end

function out = scale01(z)

    out = (z-min(z))/(max(z)-min(z));
end

function [Mx,My] = streamline_SPM(XP,YP,XB,YB,phi,SVec)

% STREAMLINE_SPM is a function that will compute the geometric integrals at
% a point (XP, YP) away from an airfoil due to source panels on the
% airfoil. 
%
% Credit for this algorithm goes to jte0419 on GitHub, whose original 
% code can be found here: 
%   
% https://github.com/jte0419/Panel_Methods
%
% A YouTube video explaining the derivations behind the following formulae
% can be found at: 
%
% https://www.youtube.com/watch?v=BnPZjGCatcg
%
% INPUTS:
% - XP: x-coordinate of computation point.
%
% - YP: y-coordinate of computation point.
%
% - XB: x-coordinate of boundary points.
%
% - YB: y-coordinate of boundary points.
%
% - phi: [numPanels x 1] column vector containing the angles between the
%        positive x-axis and interiors of the panels [radians].
%
% - SVec: [numPanels x 1] column vector containing the lengths of all the
%         panels. 
% 
% OUTPUTS:
% - Mx: [numPanels x 1] vector of x-direction source panel geometric integrals. 
%
% - My: [numPanels x 1] vector of y-direction source panel geometric integrals. 

% Extract the number of panels.
numPanels = length(XB)-1;                                                      

% Initialize Mx and My geometric integral arrays.
Mx = zeros(numPanels,1);                                                      
My = zeros(numPanels,1);                                                       

% Compute Mx and My
for j = 1:numPanels                                                      
    % Compute the appropriate terms that appear in the derivations of Mx
    % and My. 
    A  = -(XP-XB(j)) * cos(phi(j)) - (YP-YB(j))*sin(phi(j));                    
    B  = (XP-XB(j))^2 + (YP-YB(j))^2;                                        
    
    % Define C and D terms for x-direction geometric integral.
    Cx = -cos(phi(j));                                                     
    Dx = XP - XB(j);

    % Define C and D terms for y-direction geometric integral. 
    Cy = -sin(phi(j));                                                      
    Dy = YP - YB(j);                                                   
    
    E  = sqrt(B-A^2); 

    % Set E to zero if E is not real. 
    if (~isreal(E))
        E = 0;
    end
    
    % Compute Mx
    term1 = 0.5 * Cx * log((SVec(j)^2 + 2*A*SVec(j) + B)/B); 

    term2 = ((Dx-A*Cx)/E)*(atan2((SVec(j)+A),E) - atan2(A,E));                 
    
    % Compose geometric integral for x-direction. 
    Mx(j) = term1 + term2;                                                  
    
    % Compute My
    term1 = 0.5*Cy*log((SVec(j)^2 + 2*A*SVec(j) + B) / B); 

    term2 = ((Dy-A*Cy)/E) * ( atan2( (SVec(j)+A),E ) - atan2(A,E) );                 
    
    % Compose geometric integral for y-direction.
    My(j) = term1 + term2;                                                  
    
    % Set NANs, INFs, or imaginary numbers to zero in Mx and My arrays 
    if (isnan(Mx(j)) || isinf(Mx(j)) || ~isreal(Mx(j)))
        Mx(j) = 0;
    end
    
    if (isnan(My(j)) || isinf(My(j)) || ~isreal(My(j)))
        My(j) = 0;
    end
end

end 

function [Nx,Ny] = streamline_VPM(XP,YP,XB,YB,phi,SVec)

% STREAMLINE_VPM is a function that will compute the geometric integrals at
% a point (XP, YP) away from an airfoil due to vortex panels on the
% airfoil. 
%
% Credit for this algorithm goes to jte0419 on GitHub, whose original 
% code can be found here: 
%   
% https://github.com/jte0419/Panel_Methods
%
% A YouTube video explaining the derivations behind the following formulae
% can be found at: 
%
% https://www.youtube.com/watch?v=TBwBnW87hso
%
% INPUTS
% - XP: x-coordinate of computation point.
%
% - YP: y-coordinate of computation point.
%
% - XB: x-coordinate of boundary points.
%
% - YB: y-coordinate of boundary points.
%
% - phi: A [numPanels x 1] column vector containing the angles between the 
%        positive x-axis and interior of panel [radians].
%
% - SVec: A [numPanels x 1] column vector containing the lengths of all the
%         panels. 
% 
% OUTPUTS
% - Nx: [numPanels x 1] array of x-direction vortex panel geometric integrals. 
%
% - Ny: [numPanels x 1] array of y-direction vortex panel geometric integrals. 


% Extract the number of panels.
numPanels = length(XB)-1;                                             

% Initialize Nx and Ny geometric integral arrays. 
Nx = zeros(numPanels,1);                                                   
Ny = zeros(numPanels,1);                                                  

% Calculate Nx and Ny using the derived formulae. 
for j = 1:1:numPanels                                                         
    % Compute intermediate values.
    A  = -(XP-XB(j)) * cos(phi(j)) - (YP-YB(j))*sin(phi(j));                   
    B  = (XP-XB(j))^2 + (YP-YB(j))^2;                                        
    
    % Define C and D terms for x-direction integrals.
    Cx = sin(phi(j));                                                       
    Dx = -(YP-YB(j));                                                      

    % Define C and D terms for y-direction integrals. 
    Cy = -cos(phi(j));                                                     
   
    Dy = XP-XB(j);                                                          
   
    E  = sqrt(B-A^2);                                                      
    
    % Make E = 0 if E is not real. 
    if (~isreal(E))
        E = 0;
    end
    
    % Compute Nx
    term1 = 0.5*Cx*log( (SVec(j)^2 + 2*A*SVec(j)+B)/B );                              
    term2 = ((Dx-A*Cx)/E) * ( atan2((SVec(j)+A),E) - atan2(A,E) );                 
   
    % Formulate Nx geometric integral. 
    Nx(j) = term1 + term2;                                                  
    
    % Compute Ny
    term1 = 0.5*Cy*log((SVec(j)^2+2*A*SVec(j)+B)/B);                             
    term2 = ( (Dy-A*Cy)/E )*( atan2((SVec(j)+A),E) - atan2(A,E) );                 
    
    % Formulate Ny geometric integral. 
    Ny(j) = term1 + term2;                                                  
    
	% Turn NANs, INFs, or imaginary numbers to zeros. 
    if (isnan(Nx(j)) || isinf(Nx(j)) || ~isreal(Nx(j)))
        Nx(j) = 0;
    end
    if (isnan(Ny(j)) || isinf(Ny(j)) || ~isreal(Ny(j)))
        Ny(j) = 0;
    end
end

end 

function [I,J] = computeIJ_SPM(XC,YC,XB,YB,phi,SVec)
numPanels = length(XC);                                                        
I = zeros(numPanels,numPanels);                                                   
J = zeros(numPanels,numPanels);                                                  
for i = 1:1:numPanels                                                        
    for j = 1:1:numPanels                                                
        if (j ~= i)                                            
            % Compute intermediate values
            A  = -(XC(i)-XB(j)) * cos(phi(j))-(YC(i)-YB(j)) * sin(phi(j));   
            B  = (XC(i)-XB(j))^2+(YC(i)-YB(j))^2;

            % Define C and D coefficients for the panel-normal direction. 
            Cn = sin(phi(i)-phi(j));                                     
            Dn = -(XC(i)-XB(j))*sin(phi(i))+(YC(i)-YB(j)) * cos(phi(i)); 

            % Define C and D coefficients for the panel-tangential
            % direction.
            Ct = -cos(phi(i)-phi(j));                                    
            Dt = (XC(i)-XB(j))*cos(phi(i))+(YC(i)-YB(j)) * sin(phi(i));     
            
            E  = sqrt(B-A^2);                                               
            
            % Set E = 0 if E is not real. 
            if (~isreal(E))
                E = 0;
            end
            
            % Determine I (needed for normal velocity).
            term1  = 0.5*Cn*log((SVec(j)^2+2*A*SVec(j)+B)/B);                     
            term2  = ( (Dn-A*Cn)/E )*( atan2((SVec(j)+A),E) - atan2(A,E) );        
            
            I(i,j) = term1 + term2;                                        
            
            % Determine J (needed for tangential velocity).
            term1  = 0.5*Ct*log((SVec(j)^2+2*A*SVec(j)+B)/B);                     
            term2  = ( (Dt-A*Ct)/E )*( atan2((SVec(j)+A),E) - atan2(A,E) );        
            
            J(i,j) = term1 + term2;                                       
        end
        
        % Change NANs, INFs, or imaginary numbers to zeros.
        if (isnan(I(i,j)) || isinf(I(i,j)) || ~isreal(I(i,j)))
            I(i,j) = 0;
        end

        if (isnan(J(i,j)) || isinf(J(i,j)) || ~isreal(J(i,j)))
            J(i,j) = 0;
        end
    end
end
end

function [K,L] = computeKL_VPM(XC,YC,XB,YB,phi,SVec)

numPanels = length(XC);                                                      

% Initialize arrays
K = zeros(numPanels, numPanels);                                                   
L = zeros(numPanels, numPanels); 


% Compute integral
for i = 1:1:numPanels                                                          
    for j = 1:1:numPanels                                                     
        if (j ~= i)                                                         
            A  = -(XC(i)-XB(j))*cos(phi(j))-(YC(i)-YB(j))*sin(phi(j));      
            B  = (XC(i)-XB(j))^2+(YC(i)-YB(j))^2;                           
           
            % Determine normal C and D terms. 
            Cn = -cos(phi(i)-phi(j));                                    
            Dn = (XC(i)-XB(j))*cos(phi(i))+(YC(i)-YB(j))*sin(phi(i));      
            
            % Determine tangential C and D terms. 
            Ct = sin(phi(j)-phi(i));                                       
            Dt = (XC(i)-XB(j))*sin(phi(i))-(YC(i)-YB(j))*cos(phi(i));     
            
            E  = sqrt(B-A^2);                                             
            
            % If E is not real, then set E = 0.
            if (~isreal(E))
                E = 0;
            end
            
            % Compute K matrix element
            term1  = 0.5*Cn * log((SVec(j)^2+2*A*SVec(j)+B)/B);                     
            term2  = ((Dn-A*Cn)/E) * (atan2( (SVec(j)+A),E) - atan2(A,E) );         
            K(i,j) = term1 + term2;                                        
            
            % Compute L matrix element
            term1  = 0.5*Ct*log((SVec(j)^2+2*A*SVec(j)+B)/B);                     
            term2  = ((Dt-A*Ct)/E)*(atan2((SVec(j)+A),E)-atan2(A,E));         
            L(i,j) = term1 + term2;                                   
        end
        
        % Set NANs, INFs, or imaginary numbers to zero.
        if (isnan(K(i,j)) || isinf(K(i,j)) || ~isreal(K(i,j)))
            K(i,j) = 0;
        end

        if (isnan(L(i,j)) || isinf(L(i,j)) || ~isreal(L(i,j)))
            L(i,j) = 0;
        end
    end
end
end

end

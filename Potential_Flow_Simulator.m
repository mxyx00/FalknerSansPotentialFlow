function potentialFlowPlot = panelCode(airfoil,aoa,velocity, resolution)

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

%% Setup
% Opening an airfoil data file. Reading bonyname and creating vector.

fid = fopen(airfoil);
xycell = textscan(fid, '%f %f','headerlines', 1); % Skipping titles and text
xy=cell2mat(xycell); % 
fclose(fid);
bodyname=textread(airfoil,'%s',1,'delimiter','\n');

%% Preliminary Calculations
% Converting aoa from degrees to radians. Removing double-used points such
% as LE and TE, re-sorting direction of vector.

U = velocity; % free stream velocity
aoa=aoa*pi/180; %Converting angle of attack to radians. 

%Separating coordinates into x and y columns.
x=xy(:,1);
y=xy(:,2);


%Delete last data point (Redundant leading edge data point)
x(end)=[];
y(end)=[];

 
LE=find(x==min(abs(x)));

x=[flipud(x(1:LE)); flipud(x(LE:end))]*100;
y=[flipud(y(1:LE)); flipud(y(LE:end))]*100;

% x(end)=[];
% y(end)=[];

xy = [x,y]

% Parameterization
% Using spline function and cummulative distance in order to create an
% airfoil with collocation points incrementally further apart from each
% other.

% Splitting the airfoil points so they are symmetrical.
dist = vecnorm(diff(xy,1,1), 2, 2);

%Parameterize x and y coordinates
abscissa = [0; cumsum(dist)];

%Compute spline coefficients
pp = spline(abscissa, flip(xy',2));
suction = flip(coord(0, abscissa(LE),resolution/2)); % Spacing
pressure = coord(abscissa(end),abscissa(LE),resolution/2); % Spacing


abscissa = [flip(suction), flip(pressure(1:end-1))]

fitSpline = flip(ppval(pp,abscissa),2);

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

%% Panel Plot
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


%% Plot Formatting and Etc. 

leg3=['Collocation Points ( ' num2str(length(xc)) ' )'];
legend([AF_F1(1),AF_F1(2),AF_F2],'Panels','Panels',leg3,'Location','NorthEast');

 
%% RESAMPLE

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

    maxT=max(xT)
    
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
%         figure(figure1)
        
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

            
            
       % %------------------------------ Emptying the velocities inside the airfoil
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

%% Pressure

figure()
subplot(2,1,1)

% Testing
aoa=0;
xy=cell2mat(xycell)
x=xy(:,1);
y=xy(:,2);
x(end)=[];
y(end)=[];
LE=find(x==min(abs(x)));
x=[flipud(x(1:LE)); flipud(x(LE:end))]*100;
y=[flipud(y(1:LE)); flipud(y(LE:end))]*100;
x(end)=[];
y(end)=[];
n=length(x);		% 1x1
t=find(x==max(x));
c=x(t);
xL=[x(n) ; x(1:(n-1))];	% nx1
yL=[y(n) ; y(1:(n-1))];	% nx1
xc=(xL+x)./2;		% nx1
yc=(yL+y)./2;
s=((x-xL).^2+(y-yL).^2).^0.5;	% nx1
tx=(x-xL)./s;	% nx1
ty=(y-yL)./s;	% nx1
nx=-ty;	% nx1
ny=tx;
% End of testing

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
T=[sTij vT];    % nx(n+1)

[xSorted,indx]=sort(x);
t=indx(length(indx));   % 1x1

Mnp1_row=T(t,:)+T((t+1),:); % 1x(n+1)

M=[N;Mnp1_row];   % (n+1)x(n+1)

Ux=U.*cos(aoa);   % 1x1
Uy=U.*sin(aoa);   % 1x1

bi=-(Ux.*nx+Uy.*ny);    % nx1
bnp1=-(Ux.*(tx(t)+tx(t+1))+Uy.*(ty(t)+ty(t+1)));    % 1x1

b=[bi;bnp1];    % (n+1)x1

a=M^-1*b;   % (n+1)x1

vni=N*a+(Ux.*nx+Uy.*ny);    % nx1 (=0)
vti=T*a+(Ux.*tx+Uy.*ty);    % nx1 (!=0)

Cpi(:,1)=(vti./U).^2;  % nx1

TotalGamma=2*pi*a(n+1)*sum(s);
c=x(t);
L1=1.225*U*TotalGamma; 
cL=L1/(0.5*1.225*(U^2)*c);


hold on
box on

Cpi=fliplr(Cpi)
cL=fliplr(cL);

cpiplot1(1)=plot(xc(2:t),Cpi(2:t,1),'r','LineWidth',1.5);
cpiplot1(2)=plot([ xc(t:end); xc(1)],[ Cpi(t:end,1); Cpi(1,1)],'b','LineWidth',1.5);
cpiplot2(4)=plot([0,100],[0,0],'Color',[0.5020         0         0],'LineWidth',1.5);

cpiplot2(2)=ylabel('C_p');


v=axis;
set(gca,'XTick',0:5:100)
d=abs(v(4)-v(3))/10;
d=d-mod(d,0.01);
set(gca,'YTick',v(3):d:v(4))


cpiplot2(1)=xlabel('X \rightarrow');
cpiplot2(3)=title(['C_p Distribution of ' bodyname]);
legend('Upper Surface', 'Lower Surface','Location','NorthEast' )
pbaspect([7 2 1])


%% Pressure Gradient Plot

subplot(2,1,2)
s = sqrt(xc.^2+yc.^2);
psi1 = gradient(Cpi(2:t,1),s(2:t));
psi2 = gradient([ Cpi(t:end,1); Cpi(1,1)],[s(t:end);s(1)]);
plot(s(2:t),psi1,'r','LineWidth',1.5)
hold on
plot([s(t:end);s(1)],psi2,'b','LineWidth',1.5)
xlabel('X \rightarrow');
title(['Pressure gradient with respect to distance of ' bodyname]);
legend('Upper Surface', 'Lower Surface','Location','NorthEast' )
grid on
pbaspect([7 2 1])

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
end    

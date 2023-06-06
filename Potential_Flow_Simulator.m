function potentialFlowPlot = panelCode(angleofA, airfoil)

% inputs: angle of attack, airfoil, (flow velocity)
% outputs: potential flow figure

%% Setup

fid = fopen(airfoil);
xycell = textscan(fid, '%f %f','headerlines', 1);
xy=cell2mat(xycell);
fclose(fid);

%% Ignore

handles.ExtraFile=1;
handles.isccm=0;

handles.ExtraFile_1=1;
handles.isccm_1=0;
handles.isrun=0;
handles.isnote=0;

%% RUN

% creating for loop, running through iterations of the angle of attack
U = 2;

for i = 0:0.1:angleofA
    aoa=i*pi/180;
end


x=xy(:,1);
y=xy(:,2);

x(end)=[];
y(end)=[];


LE=find(x==min(abs(x)));

x=[flipud(x(1:LE)); flipud(x(LE:end))]*100; %flipud --> reverses array (eg. 123 --> 321)
y=[flipud(y(1:LE)); flipud(y(LE:end))]*100;

x(end)=[];
y(end)=[];

handles.x=x;
handles.y=y;

n=length(x);		% 1x1
t=find(x==max(x));
handles.t=t;
c=x(t);

xL=[x(n) ; x(1:(n-1))];	% nx1
yL=[y(n) ; y(1:(n-1))];	% nx1

xc=(xL+x)./2;		% nx1
yc=(yL+y)./2;



%% Features of Given Airfoil


%% Airflow Vector

xL=[x(n) ; x(1:(n-1))];	% nx1
yL=[y(n) ; y(1:(n-1))];	% nx1

xc=(xL+x)./2;		% nx1
yc=(yL+y)./2;		% nx1

handles.xc=xc;
handles.yc=yc;

s=((x-xL).^2+(y-yL).^2).^0.5;	% nx1
tx=(x-xL)./s;	% nx1
ty=(y-yL)./s;	% nx1

nx=-ty;	% nx1
ny=tx;	% nx1


xT=x*cos(aoa)+y*sin(aoa);
yT=-x*sin(aoa)+y*cos(aoa);

handles.maxT=max(xT);

        hold on
box on
axis equal
    
% Figure Boundaries
left=-round(2.5*n);
right=100+round(2.5*n);
bottom=min(yT)-round(1.2*n);
top=max(yT)+round(1.2*n);


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
T=[sTij vT];    % nx(n+1)

Mnp1_row=T(t,:)+T((t+1),:); % 1x(n+1)

M=[N;Mnp1_row];   % (n+1)x(n+1)

Ux=U.*cos(aoa);   % 1x1
Uy=U.*sin(aoa);   % 1x1

Bi=-(Ux.*nx+Uy.*ny);    % nx1
Bnp1=-(Ux.*(tx(t)+tx(t+1))+Uy.*(ty(t)+ty(t+1)));    % 1x1

B=[Bi;Bnp1];    % (n+1)x1

A=M^-1.*B;   % (n+1)x1

vni=N*A+(Ux.*nx+Uy.*ny);    % nx1 (=0)
vti=T*A+(Ux.*tx+Uy.*ty);    % nx1 (!=0)

TotalGamma=2*pi*A(n+1)*sum(s);
L=1.225*U*TotalGamma; 

Cpi=1-(vti./U).^2;  % nx1

cL=L/(0.5*1.225*(U^2)*c);
set(handles.cl_placeholder,'String',['cL = ' num2str( round(cL*1e5)/1e5 )])

theta_i=atan(ty./tx);

th1=[theta_i(n);theta_i(1:(n-1))];
th2=[theta_i(2:n);theta_i(1)];

Psi_i=abs((th1-th2)./2);

R_i=s./(2.*sin(Psi_i./2));

xxcT=xxc*cos(aoa)+yyc*sin(aoa);
yycT=-xxc*sin(aoa)+yyc*cos(aoa);

handles.VP2=zeros(1,n);
handles.Vaoa=aoa;

    for j=1:n
        figure(finished)
        figure(handles.figure1)
        
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

        if xxcT(2,j)>min(xT) & xxcT(2,j)<max(xT)
            handles.VP2(j)=quiver(xxcT(:,j),yycT(:,j),vxT(:,j),vyT(:,j),14/n,'b'); % vertical line vectors
        else
            handles.VP2(j)=quiver(xxcT(:,j),yycT(:,j),vxT(:,j),vyT(:,j),8/n,'b'); % vertical line vectors
        end
            
        set(handles.VP2(j),'Visible','off')
    
    
    if isempty( find( strcmp( get(handles.ccm,'Checked'),'on')==1 ) )==1
        handles.VP4=quiver(xxcT,yycT,vxT,vyT,1.25,'b'); %
    else
        handles.VP4=quiver(xxcT,yycT,vxT,vyT,0.8,'b'); %
    end


%% Coefficient of Pressure Distribution

S=-16;
D=1;
F=16;

alpha=(S:D:F)*pi/180;

for i=1:length(alpha)
       
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

    Ux=U.*cos(alpha(i));   % 1x1
    Uy=U.*sin(alpha(i));   % 1x1

    bi=-(Ux.*nx+Uy.*ny);    % nx1
    bnp1=-(Ux.*(tx(t)+tx(t+1))+Uy.*(ty(t)+ty(t+1)));    % 1x1

    b=[bi;bnp1];    % (n+1)x1

    a=M^-1*b;   % (n+1)x1

    vni=N*a+(Ux.*nx+Uy.*ny);    % nx1 (=0)
    vti=T*a+(Ux.*tx+Uy.*ty);    % nx1 (!=0)

   handles.Cpi(:,i)=1-(vti./U).^2;  % nx1

   TotalGamma=2*pi*a(n+1)*sum(s);
   c=x(t);
   L1=1.225*U*TotalGamma; 
   handles.cL(i)=L1/(0.5*1.225*(U^2)*c);

%    clc

end

Cpi=fliplr(handles.Cpi);
cL=fliplr(handles.cL);

v=axis;
set(gca,'XTick',0:5:100)
d=abs(v(4)-v(3))/20;
d=d-mod(d,0.01);
set(gca,'YTick',v(3):d:v(4))

handles.cpiplot2(1)=xlabel('X \rightarrow');
handles.cpiplot2(3)=title(['C_p Distribution of ' bodyname]);

set(handles.cpiplot1,'Visible','off')
set(handles.cpiplot2,'Visible','off')
grid on



close(finished)
figure(handles.figure1)


axes(handles.axes4)
handles.cLalpha(1)=plot(alpha*180/pi,handles.cL','r','LineWidth',1.5);
hold on
box on
set(gca,'XTick',-16:2:16)
v=axis;
axis([-16 16 v(3) v(4)])
set(gca,'YTick',v(3):0.25:v(4))
v3=min(handles.cL);
v4=max(handles.cL);

if mod(abs(v3),1)<=0.25
    v3=ceil(v3)-0.25;
elseif mod(abs(v3),1)<=0.5
    v3=ceil(v3)-0.5;
elseif mod(abs(v3),1)<=0.75
    v3=ceil(v3)-0.75;
else
    v3=floor(v3);
end
    
if mod(v4,1)<=0.25
    v4=floor(v4)+0.25;
elseif mod(v4,1)<=0.5
    v4=floor(v4)+0.5;
elseif mod(v4,1)<=0.75
    v4=floor(v4)+0.75;
else
    v4=ceil(v4);
end

axis([-16 16 v3 v4])

handles.cLalpha(3)=plot([-16 16],[0 0],'Color',[0.5020         0         0],'LineWidth',1.5);
handles.cLalpha(4)=plot([0 0],[v3 v4],'Color',[0.5020         0         0],'LineWidth',1.5);
grid on
handles.Result=[round(alpha*180/pi)' handles.cL'];



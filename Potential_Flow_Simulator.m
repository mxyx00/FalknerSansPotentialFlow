function varargout = Potential_Flow_Simulator(varargin)
% POTENTIAL_FLOW_SIMULATOR M-file for Potential_Flow_Simulator.fig
%      POTENTIAL_FLOW_SIMULATOR, by itself, creates a new POTENTIAL_FLOW_SIMULATOR or raises the existing
%      singleton*.
%
%      H = POTENTIAL_FLOW_SIMULATOR returns the handle to a new POTENTIAL_FLOW_SIMULATOR or the handle to
%      the existing singleton*.
%
%      POTENTIAL_FLOW_SIMULATOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in POTENTIAL_FLOW_SIMULATOR.M with the given input arguments.
%
%      POTENTIAL_FLOW_SIMULATOR('Property','Value',...) creates a new POTENTIAL_FLOW_SIMULATOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Potential_Flow_Simulator_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Potential_Flow_Simulator_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help Potential_Flow_Simulator

% Last Modified by GUIDE v2.5 31-Oct-2006 07:22:58

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Potential_Flow_Simulator_OpeningFcn, ...
                   'gui_OutputFcn',  @Potential_Flow_Simulator_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

% End initialization code - DO NOT EDIT


% --- Executes just before Potential_Flow_Simulator is made visible.
function Potential_Flow_Simulator_OpeningFcn(hObject, eventdata, handles)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Potential_Flow_Simulator (see VARARGIN)

% Choose default command line output for Potential_Flow_Simulator
handles.output = hObject;

handles.Files_old = dir(['Airfoil_Coordinate_Files' filesep 'Current_Files' filesep '*.dat']);
fid = fopen(['Airfoil_Coordinate_Files' filesep 'Current_Files' filesep 'Extra.dat'],'w');
fprintf(fid,'This Just an additional text file');
fclose(fid);
delete(handles.ccm)

handles.Files_old_1 = dir(['Arbitrary_Body_Coordinate_Files' filesep '*.dat']);
fid = fopen(['Arbitrary_Body_Coordinate_Files' filesep 'Extra.dat'],'w');
fprintf(fid,'This Just an additional text file');
fclose(fid);
delete(handles.ccm_1)

handles.ExtraFile=1;
handles.isccm=0;

handles.ExtraFile_1=1;
handles.isccm_1=0;
handles.isrun=0;
handles.isnote=0;

% Update handles structure
guidata(hObject, handles);


% UIWAIT makes Potential_Flow_Simulator wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = Potential_Flow_Simulator_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% ####################### --------------- >>>>  Menu (File) 

function file_Callback(hObject, eventdata, handles)
% hObject    handle to file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = guidata(gcbo);

%.........................Listing the airfoils in menu


if handles.isrun==0;

    handles.Files_new = dir(['Airfoil_Coordinate_Files' filesep 'Current_Files' filesep '*.dat']);

    set(handles.Empty,'Visible','off')

    if handles.ExtraFile==1 & length(handles.Files_new)==1
        set(handles.Empty,'Visible','on')
    end


    if isempty(handles.Files_new)==1
        set(handles.Empty,'Visible','on')
    end

    if isequal(handles.Files_new,handles.Files_old)==0
        % Changing Menu if files are added, deleted or renamed

        handles.ischanged=1;

        if handles.ExtraFile==1
            delete(['Airfoil_Coordinate_Files' filesep 'Current_Files' filesep 'Extra.dat'])
            handles.ExtraFile=0;
            handles.Files_new = dir(['Airfoil_Coordinate_Files' filesep 'Current_Files' filesep '*.dat']);
        end

        if handles.isccm==1
            set(handles.ccm,'Visible','off')
            handles.isccm=0;
        end

        for i=1:length(handles.Files_new)
            handles.ccm(i) = uimenu(handles.file_chooseaf,'Tag',['ccm' num2str(i)],'Label',handles.Files_new(i).name, ...
                'Callback','Potential_Flow_Simulator(''ccm_Callback'',gcbo,[],guidata(gcbo))','Checked','off');

            if isempty(handles.Files_new)==1
                handles.isccm=0;
            else
                handles.isccm=1;
            end

        end

    else
        handles.ischanged=0;
    end

    handles.Files_old=handles.Files_new;


    % Listing the Arbitrary Bodies in menu ####################################

    handles.Files_new_1 = dir(['Arbitrary_Body_Coordinate_Files' filesep '*.dat']);

    set(handles.Empty_1,'Visible','off')

    if handles.ExtraFile_1==1 & length(handles.Files_new_1)==1
        set(handles.Empty_1,'Visible','on')
    end


    if isempty(handles.Files_new_1)==1
        set(handles.Empty_1,'Visible','on')
    end

    if isequal(handles.Files_new_1,handles.Files_old_1)==0
        % Changing Menu if files are added, deleted or renamed

        handles.ischanged_1=1;

        if handles.ExtraFile_1==1
            delete(['Arbitrary_Body_Coordinate_Files' filesep 'Extra.dat'])
            handles.ExtraFile_1=0;
            handles.Files_new_1 = dir(['Arbitrary_Body_Coordinate_Files' filesep '*.dat']);
        end

        if handles.isccm_1==1
            set(handles.ccm_1,'Visible','off')
            handles.isccm_1=0;
        end

        for i=1:length(handles.Files_new_1)
            handles.ccm_1(i) = uimenu(handles.file_choosearbitrarybody,'Tag',['ccm_1' num2str(i)],'Label',handles.Files_new_1(i).name, ...
                'Callback','Potential_Flow_Simulator(''ccm_1_Callback'',gcbo,[],guidata(gcbo))','Checked','off');
            handles.Files_old_1=handles.Files_new_1;

            if isempty(handles.Files_new)==1
                handles.isccm_1=0;
            else
                handles.isccm_1=1;
            end

        end

    else
        handles.ischanged_1=0;
    end

    handles.Files_old_1=handles.Files_new_1;

end 
guidata(gcbo,handles)     
        
% ####################### --------------- >>>>  Menu (File ->Choose Airfoil Co-ordinate file )

function file_chooseaf_Callback(hObject, eventdata, handles)
% hObject    handle to file_chooseaf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% ####################### --------------- >>>>  Menu (File -> Choose Airfoil Co-ordinate file -> Airoil File Name )
% ................................Getting the Selected Airfoil file name
function ccm_Callback(hObject, eventdata, handles)
handles = guidata(gcbo);
set(handles.text7,'String','Airfoil :')

handles.FileName=get(hObject,'Label');
handles.Dir=['Airfoil_Coordinate_Files' filesep 'Current_Files' filesep];
bodyname=textread([handles.Dir handles.FileName],'%s',1,'delimiter','\n');


if strcmp(get(handles.Empty,'Visible'),'on' )~=1
    set(handles.file_placeholder,'String',['  ' char(bodyname)])
    
    if handles.ischanged==0 & handles.isccm==1
        set(handles.ccm,'Checked','off')
        set(handles.ccm_1,'Checked','off')
    end
    
    set(hObject,'Checked','on')
end

guidata(gcbo,handles)


% ####################### --------------- >>>>  Menu (File ->Choose Arbitrary Body Co-ordinate file )
function file_choosearbitrarybody_Callback(hObject, eventdata, handles)
% hObject    handle to file_choosearbitrarybody (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% ####################### --------------- >>>>  Menu (File -> Choose Arbitrary Body Co-ordinate file -> Arbitrary Body Name )
% ................................Getting the Selected Arbitrary Body file name
function ccm_1_Callback(hObject, eventdata, handles)
handles = guidata(gcbo);
set(handles.text7,'String','Body :')

handles.FileName=get(hObject,'Label');
handles.Dir=['Arbitrary_Body_Coordinate_Files' filesep];
bodyname=textread([handles.Dir handles.FileName],'%s',1,'delimiter','\n');

if strcmp(get(handles.Empty_1,'Visible'),'on' )~=1
    set(handles.file_placeholder,'String',['  ' char(bodyname) ])
    
    if handles.ischanged_1==0 & handles.isccm_1==1
        set(handles.ccm_1,'Checked','off')
        set(handles.ccm,'Checked','off')
    end
    
    set(hObject,'Checked','on')
end

guidata(gcbo,handles)


% ####################### --------------- >>>>  Menu (File ->Run)
% ............................................ CORE PANEL METHOD CODE 
function run_Callback(hObject, eventdata, handles)
% hObject    handle to run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(gcbo);
handles.whichplot=1;
handles.isrun=1;
load alltools
if isempty( find( strcmp( get(handles.ccm,'Checked'),'on')==1 ) )==1 & isempty( find( strcmp( get(handles.ccm_1,'Checked'),'on')==1 ) )==1
    msgbox('Please Choose Co-ordinate file of Airfoil or Arbitrary Body','Which Airfoil or Arbitrary Body....?','custom',Twhataf,Tcolormap)
    return
end

  
set(handles.figure1,'Pointer','watch')
AOA=16:-1:-16;
i=get(handles.GET_AOA,'Value');

aoa=AOA(i)*pi/180;
U=2;
fid = fopen([handles.Dir handles.FileName], 'r');
xycell = textscan(fid, '%f %f','headerlines', 1);
xy=cell2mat(xycell);
fclose(fid);

x=xy(:,1);
y=xy(:,2);

% x([1 end])=1;
% y([1 end])=0;


% Problem in Airfoil Coordinate file... !!!
if isempty( find( strcmp( get(handles.ccm_1,'Checked'),'on')==1 ) )==1
    if x(1)~=x(end) | y(1)~=y(end)
       
        x([1 end])=1;
        y([1 end])=( y(1)+y(end) )/2;
        
        handles.isnote=1;
    else
        handles.isnote=0;
    end
    
end
       

x(end)=[];
y(end)=[];


LE=find(x==min(abs(x)));

x=[flipud(x(1:LE)); flipud(x(LE:end))]*100;
y=[flipud(y(1:LE)); flipud(y(LE:end))]*100;

x(end)=[];
y(end)=[];

handles.x=x;
handles.y=y;

filename=handles.FileName(1:length(handles.FileName)-4);

n=length(x);		% 1x1
t=find(x==max(x));
handles.t=t;
c=x(t);

xL=[x(n) ; x(1:(n-1))];	% nx1
yL=[y(n) ; y(1:(n-1))];	% nx1

xc=(xL+x)./2;		% nx1
yc=(yL+y)./2;	

handles.Cpi=zeros(length(xc),33);
axes(handles.axes1)

set(handles.axes1,'color',[0.7490    0.7490    0.7569])
hold on
box on

%############################################## Airfoil Features

handles.AF_F1=zeros(n,1);

handles.count_even=0;
handles.count_odd=0;

for i=1:n
    if i==n
        ip1=1;
    else
        ip1=i+1;
    end
    
    if mod(i,2)==0
        handles.AF_F1(i)=plot([x(i) x(ip1)],[y(i), y(ip1)],'Color','b','LineWidth',2);
        handles.count_even=handles.count_even+1;
    else
        handles.AF_F1(i)=plot([x(i) x(ip1)],[y(i), y(ip1)],'Color',[1 1 1],'LineWidth',4);
        handles.count_odd=handles.count_odd+1;
    end
end

set(handles.AF_F1,'Visible','off')

handles.AF_F2=plot(xc,yc,'rO','LineWidth',2,'MarkerEdgeColor','r','MarkerFaceColor','g','MarkerSize',5);
set(handles.AF_F2,'Visible','off')

axis equal
grid on
bodyname=char( get(handles.file_placeholder,'String') );
handles.AF_F3(1)=title(['Profile and Mesh Details of ',bodyname ]);
handles.AF_F3(2)=xlabel('X \rightarrow');
handles.AF_F3(3)=ylabel('Y \rightarrow');
set(handles.AF_F3,'Visible','off')

%############################################## Vecotor Plot

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

if  strcmp( get(handles.options_run_ineedvectorplot,'Checked'), 'on' ) ==1
    
    finished = waitbar(0,'Calculating and plotting the Velocity vector Field...');

    xT=x*cos(aoa)+y*sin(aoa);
    yT=-x*sin(aoa)+y*cos(aoa);
    
    handles.maxT=max(xT);

    if isempty( find( strcmp( get(handles.ccm,'Checked'),'on')==1 ) )==1
        left=-round(2.5*n);
        right=100+round(2.5*n);
        bottom=min(yT)-round(1.2*n);
        top=max(yT)+round(1.2*n);
    else
        left=-round(3*n);
        right=100+round(3*n);
        bottom=min(yT)-round(2*n);
        top=max(yT)+round(2*n);
    end
    
    handles.Outer=[left right+4 bottom top];

    if max(yT)>40 | min(yT)<-40
        handles.Focussed=[-80  150  min(yT)-40   max(yT)+30];
    else
        handles.Focussed=[-32.5679  136.5679  min(yT)-40   max(yT)+30];
    end
        
    axes(handles.axes2)
    hold on
    box on
    axis equal
    
    handles.VP1(1)=title(['Vector Plot of ' bodyname]);
    handles.VP1(2)=xlabel('X \rightarrow');
    handles.VP1(3)=ylabel('Y \rightarrow');

    handles.VP1(4)=plot([left right+4 right+4 left ,left], [bottom bottom  top top ,bottom],'r','LineWidth',3);
    set(handles.VP1,'Visible','off')

    axis(handles.Outer)
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
    T=[sTij vT];    % nx(n+1)

    Mnp1_row=T(t,:)+T((t+1),:); % 1x(n+1)

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
        waitbar(j/n,finished)

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
     

        if isempty( find( strcmp( get(handles.ccm,'Checked'),'on')==1 ) )==1
            if max(yT)>40 | min(yT)<-40
                tolerance_over_surface=4;
            else
                tolerance_over_surface=2;
            end

        else
            tolerance_over_surface=1;
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
    end
    close(finished)
    
    if isempty( find( strcmp( get(handles.ccm,'Checked'),'on')==1 ) )==1
        handles.VP4=quiver(xxcT,yycT,vxT,vyT,1.25,'b'); %
    else
        handles.VP4=quiver(xxcT,yycT,vxT,vyT,0.8,'b'); %
    end
    
    handles.VP3(:,:)=fill([xT;xT(1)],[yT;yT(1)],[0.5020         0    0.2510]);

    set(handles.VP2,'Visible','off')
    set(handles.VP4,'Visible','off')
    set(handles.VP3,'Visible','off')
    axis(handles.Focussed)
    zoom reset

end
 
%--------------------- Cp Distribution

S=-16;
D=1;
F=16;

finished = waitbar(0,'Calculating the cL & Cp Distribution for Range of AOA...');

i=0;
alpha=(S:D:F)*pi/180;

for i=1:length(alpha)
    
    waitbar(i/length(alpha),finished)
   
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

axes(handles.axes3)
hold on
box on

Cpi=fliplr(handles.Cpi);
cL=fliplr(handles.cL);

i=get(handles.GET_AOA,'Value');

if strcmp( get(handles.plot_1_cp,'Checked') ,'on')==1
    handles.cpiplot1(1)=plot(handles.xc(2:handles.t),1-Cpi(2:handles.t,i),'r','LineWidth',1.5);
    handles.cpiplot1(2)=plot([ handles.xc(handles.t:end); handles.xc(1)],1-[ Cpi(handles.t:end,i); Cpi(1,i)],'b','LineWidth',1.5);
    handles.cpiplot2(4)=plot([0,100],[1,1],'Color',[0.5020         0         0],'LineWidth',1.5);
    
    handles.cpiplot2(2)=ylabel('1-C_p');
else
    handles.cpiplot1(1)=plot(handles.xc(2:handles.t),Cpi(2:handles.t,i),'r','LineWidth',1.5);
    handles.cpiplot1(2)=plot([ handles.xc(handles.t:end); handles.xc(1)],[ Cpi(handles.t:end,i); Cpi(1,i)],'b','LineWidth',1.5);
    handles.cpiplot2(4)=plot([0,100],[0,0],'Color',[0.5020         0         0],'LineWidth',1.5);

    handles.cpiplot2(2)=ylabel('C_p');
end

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

handles.cLalpha(5)=title(['C_L-\alpha Curve of ' bodyname]);
handles.cLalpha(6)=xlabel('\alpha \rightarrow');
handles.cLalpha(7)=ylabel('C_L \rightarrow');

set(handles.cLalpha(1),'Visible','off')
set(handles.cLalpha(3),'Visible','off')
set(handles.cLalpha(5),'Visible','off')
set(handles.cLalpha(6),'Visible','off')
set(handles.cLalpha(7),'Visible','off')

set(handles.axes4,'Visible','off')
legend off

set(handles.GET_AOA,'Visible','off')
set(handles.vectorplot_aoa,'Visible','off')
set(handles.AOA_Title,'Visible','off')
set(handles.cl_placeholder,'Visible','off')
set(handles.aoaforvectorplot,'Visible','off')
set(handles.options_run_ineedvectorplot,'Enable','off')


set(handles.run,'Enable','off')
set(handles.file_restart,'Enable','on')
set(handles.ccm,'Enable','off')
set(handles.ccm_1,'Enable','off')

pan off
zoom off
%panzoom off

set([handles.options_vector_plot_1 handles.options_vector_plot_2],'Enable','off')
set([handles.options_image_capture_1 handles.options_image_capture_2 handles.options_image_capture_3],'Enable','on')
set(handles.axes1,'Visible','on')
set(handles.AF_F1,'Visible','on')
set(handles.AF_F2,'Visible','on')
set(handles.AF_F3,'Visible','on')


set(handles.view_plotwhat_meshedairfoil,'Enable','on')

if  strcmp( get(handles.options_run_ineedvectorplot,'Checked'), 'on' ) ==1
    set(handles.view_plotwhat_vectorplot,'Enable','on')
end

set(handles.view_plotwhat_cpdistribution,'Enable','on')
set(handles.view_plotwhat_clalphacurve,'Enable','on')

set(handles.Tzoom,'Enable','on')
set(handles.Tpan,'Enable','on')
set(handles.Tpanzoom,'Enable','on')
set(handles.Tdatacursor,'Enable','on')
set(handles.Ttable,'Enable','on')
set(handles.Tcapture,'Enable','on')


set(handles.view_plotwhat_meshedairfoil,'Checked','on')
set(handles.view_plotwhat_vectorplot,'Checked','off')
set(handles.view_plotwhat_cpdistribution,'Checked','off')
set(handles.view_plotwhat_clalphacurve,'Checked','off')

set(handles.plot_cp,'Enable','off')
set(handles.plot_1_cp,'Enable','off')

axes(handles.axes1)
leg3=['Collocation Points ( ' num2str(length(handles.xc)) ' )'];
legend([handles.AF_F1(1),handles.AF_F1(2) handles.AF_F2],'Panels','Panels',leg3,'Location','NorthEast');

guidata(gcbo,handles)
set(handles.figure1,'Pointer','arrow')

if handles.isnote==1
    
    msgbox(['The selected Airfoil does not have an unique Trailing Edge                                ';...
            '                                                                                          ';...
            'For Simulation, aft most 2 points have been shrunk as Trailing Edge ( mid point of the 2 )'],...
            'Selected Airfoil does not have unique Trailing Edge','custom',Ttemsg,Tcolormap)
end


%  run over

% ####################### --------------- >>>>  Menu (File ->Restart)
function file_restart_Callback(hObject, eventdata, handles)
% hObject    handle to file_restart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(gcbo);
handles.isrun=0;

set(handles.ccm,'Checked','off')
set(handles.ccm_1,'Checked','off')

set(handles.options_vector_plot_2,'Checked','off')
set(handles.options_vector_plot_1,'Checked','on')

pan off
zoom off
zoom out
%panzoom off

set( [ handles.view_plotwhat_vectorplot handles.view_plotwhat_cpdistribution ...
    handles.view_plotwhat_clalphacurve] , 'Checked','off')

set( handles.view_plotwhat_meshedairfoil,'Checked','on')

set(handles.Tzoom,'Enable','off')
set(handles.Tpan,'Enable','off')
set(handles.Tpanzoom,'Enable','off')
set(handles.Tdatacursor,'Enable','off')
set(handles.Ttable,'Enable','off')
set(handles.Tcapture,'Enable','off')
set(handles.Tmovie,'Enable','off')

set(handles.file_placeholder,'String','   < Airfoil Name >')

set(handles.options_image_capture_1,'Enable','off')
set(handles.options_image_capture_2,'Enable','off')
set(handles.options_image_capture_3,'Enable','off')

set(handles.GET_AOA,'Visible','on','Enable','on')
set(handles.vectorplot_aoa,'Visible','off')
set(handles.text_AOA,'Visible','off')
set(handles.aoaforvectorplot,'Visible','on')
set(handles.options_run_ineedvectorplot,'enable','on')


set(handles.run,'Enable','on')
set(handles.ccm,'Enable','on')
set(handles.ccm_1,'Enable','on')
set(handles.file_restart,'Enable','off')

set( [handles.view_plotwhat_meshedairfoil, handles.view_plotwhat_vectorplot,...
    handles.view_plotwhat_cpdistribution, handles.view_plotwhat_clalphacurve], 'Enable','off')

axes(handles.axes1)
legend off
cla

axes(handles.axes2)
cla

axes(handles.axes3)
legend off
cla

axes(handles.axes4)
cla

set(handles.axes1,'Visible','off')
set(handles.axes2,'Visible','off')
set(handles.axes3,'Visible','off')
set(handles.axes4,'Visible','off')

guidata(gcbo,handles)

% ####################### --------------- >>>>  Menu (File ->Exit)
function file_exit_Callback(hObject, eventdata, handles)
% hObject    handle to file_exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(gcbo);
delete(handles.figure1);
clc
a=0;
save Extra.dat a
delete('Extra.dat')
clear
closereq
load alltools
msg=['  Please Review in mathworks.com or Send your Feedbacks to author   ';...
     '                                                                    ';...
     '     j.divahar@yahoo.com                                            ';...
     '     HomePage: jdivahar.com                                         '];
    

%button = msgbox(msg,'Thank you For Trying !','custom',Tvanakam,Tcolormap);

% ####################### --------------- >>>>  Menu (View)
function view_Callback(hObject, eventdata, handles)
% hObject    handle to view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% ####################### --------------- >>>>  Menu (View ->Plot What ?)
function view_plotwhat_Callback(hObject, eventdata, handles)
% hObject    handle to view_plotwhat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% ####################### --------------- >>>>  Menu (View ->Plot What ?->Profile and Mesh Details)
function view_plotwhat_meshedairfoil_Callback(hObject, eventdata, handles)
% hObject    handle to view_plotwhat_meshedairfoil (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(gcbo);


handles.whichplot=1;
set([handles.options_vector_plot_1 handles.options_vector_plot_2],'Enable','off')

set(handles.plot_cp,'Enable','off')
set(handles.plot_1_cp,'Enable','off')

set(handles.GET_AOA,'Visible','off')
set(handles.vectorplot_aoa,'Visible','off')
set(handles.AOA_Title,'Visible','off')
set(handles.cl_placeholder,'Visible','off')

set(handles.axes1,'Visible','on')
set(handles.axes2,'Visible','off')
set(handles.axes3,'Visible','off')
set(handles.axes4,'Visible','off')

if  strcmp( get(handles.options_run_ineedvectorplot,'Checked'), 'on' ) ==1
    set(handles.VP1,'Visible','off')
    set(handles.VP2,'Visible','off')
    set(handles.VP4,'Visible','off')
    set(handles.VP3,'Visible','off')
end
    
set(handles.Tmovie,'Enable','off')

set(handles.cpiplot1,'Visible','off')
set(handles.cpiplot2,'Visible','off')

set(handles.cLalpha(1),'Visible','off')
set(handles.cLalpha(3),'Visible','off')
set(handles.cLalpha(5),'Visible','off')
set(handles.cLalpha(6),'Visible','off')
set(handles.cLalpha(7),'Visible','off')

set(handles.view_plotwhat_meshedairfoil,'Checked','on')
set(handles.view_plotwhat_vectorplot,'Checked','off')
set(handles.view_plotwhat_cpdistribution,'Checked','off')
set(handles.view_plotwhat_clalphacurve,'Checked','off')

axes(handles.axes3)
legend off
axes(handles.axes4)
legend off
axes(handles.axes1)

set(handles.AF_F1,'Visible','on')
set(handles.AF_F2,'Visible','on')
set(handles.AF_F3,'Visible','on')

leg3=['Collocation Points ( ' num2str(length(handles.xc)) ' )'];

legend([handles.AF_F1(1),handles.AF_F1(2) handles.AF_F2],'Panels','Panels',leg3,'Location','NorthEast');
zoom off
pan off
%panzoom off
guidata(gcbo,handles)

% ####################### --------------- >>>>  Menu (View ->Plot What ?->Vector Plot)
function view_plotwhat_vectorplot_Callback(hObject, eventdata, handles)
% hObject    handle to view_plotwhat_vectorplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(gcbo);


handles.whichplot=2;

set([handles.options_vector_plot_1 handles.options_vector_plot_2],'Enable','on')

set(handles.GET_AOA,'Visible','off')
set(handles.vectorplot_aoa,'Visible','on')
set(handles.vectorplot_aoa,'String',[num2str(round(handles.Vaoa*180/pi)) ' deg'])

set(handles.plot_cp,'Enable','off')
set(handles.plot_1_cp,'Enable','off')

set(handles.AOA_Title,'Visible','on')
set(handles.cl_placeholder,'Visible','on')

set(handles.axes1,'Visible','off')
set(handles.axes2,'Visible','on')
set(handles.axes3,'Visible','off')
set(handles.axes4,'Visible','off')

set(handles.AF_F1,'Visible','off')
set(handles.AF_F2,'Visible','off')
set(handles.AF_F3,'Visible','off')

set(handles.VP1,'Visible','on')
set(handles.VP2,'Visible','off')
set(handles.VP4,'Visible','on')
set(handles.VP3,'Visible','on')
set(handles.Tmovie,'Enable','on')

set(handles.cpiplot1,'Visible','off')
set(handles.cpiplot2,'Visible','off')

set(handles.cLalpha(1),'Visible','off')
set(handles.cLalpha(3),'Visible','off')
set(handles.cLalpha(5),'Visible','off')
set(handles.cLalpha(6),'Visible','off')
set(handles.cLalpha(7),'Visible','off')

set(handles.view_plotwhat_meshedairfoil,'Checked','off')
set(handles.view_plotwhat_vectorplot,'Checked','on')
set(handles.view_plotwhat_cpdistribution,'Checked','off')
set(handles.view_plotwhat_clalphacurve,'Checked','off')

axes(handles.axes1)
legend off
axes(handles.axes3)
legend off
axes(handles.axes4)
legend off
axes(handles.axes2)
% zoom out
pan on

% axis(handles.Focussed)
% zoom out
% dfbdsf


guidata(gcbo,handles)

% ####################### --------------- >>>>  Menu (View ->Plot What ?->Cp Distribution)
function view_plotwhat_cpdistribution_Callback(hObject, eventdata, handles)
% hObject    handle to view_plotwhat_cpdistribution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(gcbo);

handles.whichplot=3;

set(handles.GET_AOA,'Enable','on')
set(handles.vectorplot_aoa,'Visible','off')

set([handles.options_vector_plot_1 handles.options_vector_plot_2],'Enable','off')

set(handles.plot_cp,'Enable','on')
set(handles.plot_1_cp,'Enable','on')

set(handles.GET_AOA,'Visible','on')
set(handles.AOA_Title,'Visible','on')
set(handles.cl_placeholder,'Visible','on')

set(handles.axes1,'Visible','off')
set(handles.axes2,'Visible','off')
set(handles.axes3,'Visible','on')
set(handles.axes4,'Visible','off')

set(handles.AF_F1,'Visible','off')
set(handles.AF_F2,'Visible','off')
set(handles.AF_F3,'Visible','off')

if strcmp( get(handles.options_run_ineedvectorplot,'Checked') , 'on')==1
    set(handles.VP1,'Visible','off')
    set(handles.VP2,'Visible','off')
    set(handles.VP4,'Visible','off')
    set(handles.VP3,'Visible','off')
    set(handles.Tmovie,'Enable','off')
end

set(handles.cpiplot1,'Visible','on')
set(handles.cpiplot2,'Visible','on')

set(handles.cLalpha(1),'Visible','off')
set(handles.cLalpha(3),'Visible','off')
set(handles.cLalpha(5),'Visible','off')
set(handles.cLalpha(6),'Visible','off')
set(handles.cLalpha(7),'Visible','off')

set(handles.view_plotwhat_meshedairfoil,'Checked','off')
set(handles.view_plotwhat_vectorplot,'Checked','off')
set(handles.view_plotwhat_cpdistribution,'Checked','on')
set(handles.view_plotwhat_clalphacurve,'Checked','off')

axes(handles.axes1)
legend off
axes(handles.axes4)
legend off
axes(handles.axes3)
legend('Upper Surface', 'Lower Surface','Location','NorthEast' )

zoom off
pan off
%panzoom off
guidata(gcbo,handles)

% ####################### --------------- >>>>  Menu (View ->Plot What ?->CL alpha Curve)
function view_plotwhat_clalphacurve_Callback(hObject, eventdata, handles)
% hObject    handle to view_plotwhat_clalphacurve (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(gcbo);

handles.whichplot=4;

set([handles.options_vector_plot_1 handles.options_vector_plot_2],'Enable','off')

axes(handles.axes4)

set(handles.GET_AOA,'Visible','off')
set(handles.vectorplot_aoa,'Visible','off')

set(handles.AOA_Title,'Visible','off')
set(handles.cl_placeholder,'Visible','off')

set(handles.plot_cp,'Enable','off')
set(handles.plot_1_cp,'Enable','off')

set(handles.axes1,'Visible','off')
set(handles.axes2,'Visible','off')
set(handles.axes3,'Visible','off')
set(handles.axes4,'Visible','on')

set(handles.AF_F1,'Visible','off')
set(handles.AF_F2,'Visible','off')
set(handles.AF_F3,'Visible','off')

if strcmp( get(handles.options_run_ineedvectorplot,'Checked'),'on')
    set(handles.VP1,'Visible','off')
    set(handles.VP2,'Visible','off')
    set(handles.VP4,'Visible','off')
    set(handles.VP3,'Visible','off')
end

set(handles.Tmovie,'Enable','off')

set(handles.cpiplot1,'Visible','off')
set(handles.cpiplot2,'Visible','off')

set(handles.cLalpha(1),'Visible','on')
set(handles.cLalpha(3),'Visible','on')
set(handles.cLalpha(5),'Visible','on')
set(handles.cLalpha(6),'Visible','on')
set(handles.cLalpha(7),'Visible','on')

set(handles.view_plotwhat_meshedairfoil,'Checked','off')
set(handles.view_plotwhat_vectorplot,'Checked','off')
set(handles.view_plotwhat_cpdistribution,'Checked','off')
set(handles.view_plotwhat_clalphacurve,'Checked','on')

axes(handles.axes1)
legend off
axes(handles.axes3)
legend off

axes(handles.axes4)


zoom off
pan off
%panzoom off

guidata(gcbo,handles)

% ####################### --------------- >>>>  Menu (View ->Author's Information)
function view_plotwhat_author_Callback(hObject, eventdata, handles)
% hObject    handle to view_plotwhat_author (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(gcbo);
if get(handles.file_placeholder,'Value')==0
    set(hObject,'Checked','on')
    set(handles.Author,'Visible','on')
    set(handles.file_placeholder,'Value',1)
else
    set(hObject,'Checked','off')
    set(handles.Author,'Visible','off')
    set(handles.file_placeholder,'Value',0)
end
    
guidata(gcbo,handles)


% ####################### --------------- >>>>  Menu (Options)
function view_options_Callback(hObject, eventdata, handles)
% hObject    handle to view_options (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% ####################### --------------- >>>>  Menu (Options ->Image Cature Settings)
function options_image_capture_Callback(hObject, eventdata, handles)
% hObject    handle to options_image_capture (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(gcbo);

guidata(gcbo,handles)


% ####################### --------------- >>>>  Menu (Options ->Image Cature Settings ->Include Inputs, Title & Axis)
function options_image_capture_1_Callback(hObject, eventdata, handles)
% hObject    handle to options_image_capture_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(gcbo);
set(handles.options_image_capture_1,'Checked','on')
set([ handles.options_image_capture_2 handles.options_image_capture_3],'Checked','off')
guidata(gcbo,handles)

% ####################### --------------- >>>>  Menu (Options ->Image Cature Settings ->Include Title & Axis)
function options_image_capture_2_Callback(hObject, eventdata, handles)
% hObject    handle to options_image_capture_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(gcbo);
set(handles.options_image_capture_2,'Checked','on')
set([ handles.options_image_capture_3 handles.options_image_capture_1],'Checked','off')
guidata(gcbo,handles)

% ####################### --------------- >>>>  Menu (Options ->Image Cature Settings ->Exclude Inputs, Title & Axis)
function options_image_capture_3_Callback(hObject, eventdata, handles)
% hObject    handle to options_image_capture_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = guidata(gcbo);
set(handles.options_image_capture_3,'Checked','on')
set([ handles.options_image_capture_1 handles.options_image_capture_2],'Checked','off')
guidata(gcbo,handles)


% ####################### --------------- >>>>  Menu (Options ->Vector Plot View Settings)
function options_vector_plot_Callback(hObject, eventdata, handles)
% hObject    handle to options_vector_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% ####################### --------------- >>>>  Menu (Options ->Vector Plot View Settings ->Focussed View)
function options_vector_plot_1_Callback(hObject, eventdata, handles)
% hObject    handle to options_vector_plot_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(gcbo);
set(handles.options_vector_plot_1,'Checked','on')
set(handles.options_vector_plot_2,'Checked','off')
guidata(gcbo,handles)

% zoom out
axis(handles.Focussed)
zoom reset
pan on

% ####################### --------------- >>>>  Menu (Options ->Vector Plot View Settings ->Full field View)
function options_vector_plot_2_Callback(hObject, eventdata, handles)
% hObject    handle to options_vector_plot_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(gcbo);
set(handles.options_vector_plot_2,'Checked','on')
set(handles.options_vector_plot_1,'Checked','off')
guidata(gcbo,handles)

axis(handles.Outer)
% zoom out
zoom reset
pan off

% ####################### --------------- >>>>  Menu (Options ->Run Settings)
function options_run_Callback(hObject, eventdata, handles)
% hObject    handle to options_run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% ####################### --------------- >>>>  Menu (Options ->Run Settings ->I need Vector Plot also)
function options_run_ineedvectorplot_Callback(hObject, eventdata, handles)
% hObject    handle to options_run_ineedvectorplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if  strcmp( get(handles.options_run_ineedvectorplot,'Checked'), 'on' ) ==1
    set(handles.options_run_ineedvectorplot,'Checked', 'off')
    set(handles.GET_AOA,'Enable', 'off')
    
else
    set(handles.options_run_ineedvectorplot,'Checked', 'on')
    set(handles.GET_AOA,'Enable', 'on')
end
    

% ####################### --------------- >>>>  Menu (Options ->Export Settings)
function export_options_Callback(hObject, eventdata, handles)
% hObject    handle to export_options (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% ####################### --------------- >>>>  Menu (Options ->Export Settings ->Export cL-Alpha Data)
function export_cl_alpha_Callback(hObject, eventdata, handles)
% hObject    handle to export_cl_alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(gcbo);
set(handles.export_cl_alpha,'Checked','on')
set(handles.export_cp_distribution,'Checked','off')
guidata(gcbo,handles)

% ####################### --------------- >>>>  Menu (Options ->Export Settings ->Export Cp Distribution)
function export_cp_distribution_Callback(hObject, eventdata, handles)
% hObject    handle to export_cp_distribution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(gcbo);
set(handles.export_cl_alpha,'Checked','off')
set(handles.export_cp_distribution,'Checked','on')
guidata(gcbo,handles)


% ####################### --------------- >>>>  Menu (Options ->Cp Distribution Settings)
function cpdistribution_settings_Callback(hObject, eventdata, handles)
% hObject    handle to cpdistribution_settings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% ####################### --------------- >>>>  Menu (Options ->Cp Distribution Settings ->Plot 1-Cp)
function plot_1_cp_Callback(hObject, eventdata, handles)
% hObject    handle to plot_1_cp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(gcbo);
set(handles.plot_cp,'Checked','off')
set(handles.plot_1_cp,'Checked','on')


if  strcmp( get(handles.axes3,'Visible'), 'on' ) ==1


    Cpi=fliplr(handles.Cpi);
    cL=fliplr(handles.cL);
    i=get(handles.GET_AOA,'Value');

    axes(handles.axes3)
    delete(handles.cpiplot1)
    delete(handles.cpiplot2)

    handles.cpiplot1(1)=plot(handles.xc(2:handles.t),1-Cpi(2:handles.t,i),'r','LineWidth',1.5,'LineWidth',1.5);
    handles.cpiplot1(2)=plot([ handles.xc(handles.t:end); handles.xc(1)],1-[ Cpi(handles.t:end,i); Cpi(1,i)],'b','LineWidth',1.5);
    handles.cpiplot2(4)=plot([0,100],[1,1],'Color',[0.5020         0         0],'LineWidth',1.5);
    handles.cpiplot2(2)=ylabel('1-C_p');
  
    bodyname=char( get(handles.file_placeholder,'String') );
    handles.cpiplot2(1)=xlabel('X \rightarrow');
    handles.cpiplot2(3)=title(['C_p Distribution of ' bodyname]);
    
    v=axis;
    set(gca,'XTick',0:5:100)
    d=abs(v(4)-v(3))/20;
    d=d-mod(d,0.01);
    set(gca,'YTick',v(3):d:v(4))
    legend('Upper Surface', 'Lower Surface','Location','NorthEast' )
    grid on
    
end

guidata(gcbo,handles)

% ####################### --------------- >>>>  Menu (Options ->Cp Distribution Settings ->Plot Cp)
function plot_cp_Callback(hObject, eventdata, handles)
% hObject    handle to plot_cp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(gcbo);
set(handles.plot_cp,'Checked','on')
set(handles.plot_1_cp,'Checked','off')

if  strcmp( get(handles.axes3,'Visible'), 'on' ) ==1

    xc=(handles.xc);
    Cpi=fliplr(handles.Cpi);
    cL=fliplr(handles.cL);

    i=get(handles.GET_AOA,'Value');

    axes(handles.axes3)
    delete(handles.cpiplot1)
    delete(handles.cpiplot2)
    
    handles.cpiplot1(1)=plot(handles.xc(2:handles.t),Cpi(2:handles.t,i),'r','LineWidth',1.5);
    handles.cpiplot1(2)=plot([ handles.xc(handles.t:end); handles.xc(1)],[ Cpi(handles.t:end,i); Cpi(1,i)],'b','LineWidth',1.5);
    handles.cpiplot2(4)=plot([0,100],[0,0],'Color',[0.5020         0         0],'LineWidth',1.5);

    handles.cpiplot2(2)=ylabel('C_p');

    bodyname=char( get(handles.file_placeholder,'String') );
    handles.cpiplot2(1)=xlabel('X \rightarrow');
    handles.cpiplot2(3)=title(['C_p Distribution of ' bodyname]);
    
    v=axis;
    set(gca,'XTick',0:5:100)
    d=abs(v(4)-v(3))/20;
    d=d-mod(d,0.01);
    set(gca,'YTick',v(3):d:v(4))
    legend('Upper Surface', 'Lower Surface','Location','NorthEast' )

    grid on
    
end

guidata(gcbo,handles)


% ####################### --------------- >>>>  Menu (Help)
function help_Callback(hObject, eventdata, handles)
% hObject    handle to help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% ####################### --------------- >>>>  Menu (Help ->About...)
function help_about_Callback(hObject, eventdata, handles)
% hObject    handle to help_about (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cd Help_files
web('about.htm','-helpbrowser')
cd ..

% ####################### --------------- >>>>  Menu (Help ->Intro to Panel Methods)
function help_intro_panel_Callback(hObject, eventdata, handles)
% hObject    handle to help_intro_panel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cd Help_files
web('intro_panel.htm','-helpbrowser')
cd ..

% ####################### --------------- >>>>  Menu (Help ->Manual)
function help_manual_Callback(hObject, eventdata, handles)
% hObject    handle to help_manual (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cd Help_files
web('manual.htm','-helpbrowser')
cd ..

% ####################### --------------- >>>>  Menu (Help ->References)
function help_references_Callback(hObject, eventdata, handles)
% hObject    handle to help_references (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cd Help_files
web('references.htm','-helpbrowser')
cd ..


% #################### All the Menus are OVER

%  TTTTTTTTTTTTTTTTTT Tool bar (GET_AOA)
function GET_AOA_Callback(hObject, eventdata, handles)
% hObject    handle to GET_AOA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns GET_AOA contents as cell array
%        contents{get(hObject,'Value')} returns selected item from GET_AOA

handles = guidata(gcbo);
bodyname=char( textread([handles.Dir handles.FileName],'%s',1,'delimiter','\n') );

if  strcmp( get(handles.axes3,'Visible'), 'on' ) ==1

    xc=(handles.xc);
    Cpi=fliplr(handles.Cpi);

    i=34-get(hObject,'Value');

    axes(handles.axes3)
    cla
    
    if strcmp( get(handles.plot_1_cp,'Checked') ,'on')==1
        handles.cpiplot1(1)=plot(handles.xc(2:handles.t),1-Cpi(2:handles.t,i),'r','LineWidth',1.5);
        handles.cpiplot1(2)=plot([ handles.xc(handles.t:end); handles.xc(1)],1-[ Cpi(handles.t:end,i); Cpi(1,i)],'b','LineWidth',1.5);
        handles.cpiplot2(4)=plot([0,100],[1,1],'Color',[0.5020         0         0],'LineWidth',1.5);
        handles.cpiplot2(2)=ylabel('1-C_p');
    else
        handles.cpiplot1(1)=plot(handles.xc(2:handles.t),Cpi(2:handles.t,i),'r','LineWidth',1.5);
        handles.cpiplot1(2)=plot([ handles.xc(handles.t:end); handles.xc(1)],[ Cpi(handles.t:end,i); Cpi(1,i)],'b','LineWidth',1.5);
        handles.cpiplot2(4)=plot([0,100],[0,0],'Color',[0.5020         0         0],'LineWidth',1.5);

        handles.cpiplot2(2)=ylabel('C_p');
    end

    set(handles.cl_placeholder,'String',['cL = ' num2str( num2str(round(handles.cL(i)*1e5)/1e5 ) )])

    handles.cpiplot2(1)=xlabel('X \rightarrow');
    handles.cpiplot2(3)=title(['C_p Distribution of ' bodyname]);
    
    v=axis;
    set(gca,'XTick',0:5:100)
    d=abs(v(4)-v(3))/20;
    d=d-mod(d,0.01);
    set(gca,'YTick',v(3):d:v(4))
    legend('Upper Surface', 'Lower Surface','Location','NorthEast' )

    grid on
    
end
guidata(gcbo,handles)

% TTTTTTTTTTTTTTTTTT Tool bar -- Executes during object creation, after setting all properties.
function GET_AOA_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GET_AOA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    ccm - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% TTTTTTTTTTTTTTTTTT Tool bar --- Executes on button press in Tmovie.
function Tmovie_Callback(hObject, eventdata, handles)
% hObject    handle to Tmovie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(gcbo);
n=length(handles.xc);

if get(hObject,'Value')==get(hObject,'Max')
    set(handles.Tdatacursor,'Enable','off')
else
    set(handles.Tdatacursor,'Enable','on')
end

set(handles.view_plotwhat_meshedairfoil,'Enable','off')
set(handles.view_plotwhat_vectorplot,'Enable','off')
set(handles.view_plotwhat_cpdistribution,'Enable','off')
set(handles.view_plotwhat_clalphacurve,'Enable','off')

i=0;
set(handles.VP4,'Visible','off')
set(handles.VP2,'Visible','on')

ii=1:3:n;
jj=1:6:n;
while get(handles.Tmovie,'Value')==get(handles.Tmovie,'Max')
    set(handles.VP2,'Visible','off')
    set(handles.VP2,'Color','b')
    i=i+1;
    ii=ii+1;
    jj=jj+1;
    
    ii(find(ii>n))=ii(find(ii>n))-n;
    jj(find(jj>n))=jj(find(jj>n))-n;

    set(handles.VP2(ii),'Visible','on')
    set(handles.VP2(jj),'Color',[0         0    0.5020])
    drawnow
end

set(handles.VP4,'Visible','on')
set(handles.VP2,'Visible','off')

set(handles.axes2,'Color','w')
set(handles.view_plotwhat_meshedairfoil,'Enable','on')
set(handles.view_plotwhat_vectorplot,'Enable','on')
set(handles.view_plotwhat_cpdistribution,'Enable','on')
set(handles.view_plotwhat_clalphacurve,'Enable','on')

guidata(gcbo,handles)


% TTTTTTTTTTTTTTTTTT Tool bar --- Executes on button press in Ttable.
function Ttable_Callback(hObject, eventdata, handles)
% hObject    handle to Ttable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(gcbo);

if strcmp( get(handles.export_cl_alpha,'Checked'),'on')==1
    exportdata='_cL_Vs_Alpha.xls';
    exportmessage='Save cL-Alpha Data as Excel File';
elseif strcmp( get(handles.export_cp_distribution,'Checked'),'on')==1
    exportdata='_Cp_Distribution.xls';
    exportmessage='Save Cp Distribution Data as Excel File';
end
    
filename=handles.FileName(1:length(handles.FileName)-4);
filename=['Results_files' filesep filename exportdata];
[file,path] = uiputfile(filename,exportmessage);

if file==0 & path==0

else
    filename=[path file];
    if strcmp( get(handles.export_cl_alpha,'Checked'),'on')==1
        d1={'Alpha','cL'}
        xlswrite(filename, d1,'Sheet1')
        xlswrite(filename, handles.Result,'Sheet1','A2')
    elseif strcmp( get(handles.export_cp_distribution,'Checked'),'on')==1
        
        first={'AOA-Row/x-Col'};
        xlswrite(filename, first,'Sheet1')
        xlswrite(filename, handles.xc,'Sheet1','A2')
        xlswrite(filename, (-16:1:16),'Sheet1','B1')        
        
        if strcmp( get(handles.plot_1_cp,'Checked') ,'on')==1
            xlswrite(filename, 1-handles.Cpi,'Sheet1','B2')
        else
            xlswrite(filename, handles.Cpi,'Sheet1','B2')
        end
            
    end
end
 
guidata(gcbo,handles)


% TTTTTTTTTTTTTTTTTT Tool bar --- Executes on button press in Tcapture.
function Tcapture_Callback(hObject, eventdata, handles)
% hObject    handle to Tcapture (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(gcbo);

if strcmp( get(handles.options_image_capture_1,'Checked'),'on')==1
    option=1;
elseif strcmp( get(handles.options_image_capture_2,'Checked'), 'on')==1
    option=2;
else
    option=3;
end


AOA=16:-1:-16;

switch handles.whichplot
    case 1
        plottype=' Meshed Airfoil ';
    case 2
        plottype=[ ' Vector Plot @ AOA ' get(handles.vectorplot_aoa,'String') ];
    case 3
        plottype=[ ' Cp Distribution @ AOA ' num2str( AOA(get(handles.GET_AOA,'Value')) ) ' deg' ];
    case 4
        plottype=[ ' cL-Alpha Curve ' ];
end

filename=handles.FileName(1:length(handles.FileName)-4);
filename=['Captured_images' filesep filename plottype '.jpg'];
[file,path] = uiputfile(filename,['Capture' plottype 'As Image']);
filename=[path,file];

set(handles.GET_AOA,'Visible','off')
aoa=16:-1:-16;

i=get(handles.GET_AOA,'Value');

switch option
    case 1
        set([ handles.Tzoom,handles.Tpan,handles.Tpanzoom,handles.Tdatacursor,handles.Ttable,handles.Tcapture,handles.Tmovie],'Visible','off')
        
        
        if handles.whichplot==2
            set(handles.vectorplot_aoa,'Visible','on')
            set(handles.text_AOA,'Visible','off')
        elseif handles.whichplot==3
            set([handles.GET_AOA,handles.vectorplot_aoa],'Visible','off')
            set(handles.text_AOA,'Visible','on','String',[num2str(aoa(i)) ' deg'])
        end

    case 2
        set([handles.Tzoom,handles.Tpan,handles.Tpanzoom,handles.Tdatacursor,handles.Ttable,handles.Tcapture,handles.Tmovie],'Visible','off')
        set([handles.vectorplot_aoa,handles.text_AOA, handles.AOA_Title, handles.GET_AOA, handles.text7, handles.file_placeholder, handles.cl_placeholder],'Visible','off')
    case 3
        set([handles.Tzoom,handles.Tpan,handles.Tpanzoom,handles.Tdatacursor,handles.Ttable,handles.Tcapture,handles.Tmovie],'Visible','off')
        set([handles.vectorplot_aoa,handles.text_AOA,handles.AOA_Title, handles.GET_AOA, handles.text7, handles.file_placeholder, handles.cl_placeholder],'Visible','off')

        eval(['set(handles.axes' num2str(handles.whichplot) ',' '''Visible''' ',' '''off''' ')'])
        axes(handles.axes1)
        legend off
        axes(handles.axes3)
        legend off
        axes(handles.axes4)
        legend off
end

set(handles.toolbar_box,'Visible','off')

if file==0 & path==0
    
else
    saveas(gca,filename)
end

set(handles.toolbar_box,'Visible','on')

if handles.whichplot==2
    set(handles.vectorplot_aoa,'Visible','on')
    set(handles.text_AOA,'Visible','off')
elseif handles.whichplot==3
    set(handles.text_AOA,'Visible','on','String',[num2str(aoa(i)) ' deg'])
    set(handles.vectorplot_aoa,'Visible','off')
end

if handles.whichplot==1
      axes(handles.axes1)
      leg3=['Collocation Points ( ' num2str(length(handles.xc)) ' )'];
      legend([handles.AF_F1(1),handles.AF_F1(2) handles.AF_F2],'Panels','Panels',leg3,'Location','NorthEast');
end

set([ handles.Tzoom,handles.Tpan,handles.Tpanzoom,handles.Tdatacursor,handles.Ttable,handles.Tcapture,handles.Tmovie],'Visible','on')

if handles.whichplot==1
    set([handles.text7, handles.file_placeholder],'Visible','on')
elseif handles.whichplot==2
    set([handles.AOA_Title, handles.text7, handles.file_placeholder, handles.cl_placeholder],'Visible','on')
elseif handles.whichplot==3
    set([handles.AOA_Title, handles.text7, handles.file_placeholder, handles.cl_placeholder],'Visible','on')
else
    set([handles.text7, handles.file_placeholder],'Visible','on')
end

eval(['set(handles.axes' num2str(handles.whichplot) ',' '''Visible''' ',' '''on''' ')'])

if strcmp( get(handles.axes3,'Visible'),'on')==1
    set([ handles.text_AOA,handles.vectorplot_aoa ],'Visible','off')
    set(handles.GET_AOA,'Visible','on')
end
guidata(gcbo,handles)


% TTTTTTTTTTTTTTTTTT Tool bar --- Executes during object creation, after setting all properties.
function Author_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Author (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    ccm - handles not created until after all CreateFcns called
handles = guidata(gcbo);
load alltools
set(hObject,'CData',Banner)
guidata(gcbo,handles)


% TTTTTTTTTTTTTTTTTT Tool bar --- Executes on button press in Author.
function Author_Callback(hObject, eventdata, handles)
% hObject    handle to Author (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Author







% TTTTTTTTTTTTTTTTTTT  ------ All the Tool Bar Button's Create Fcn's (Mapping image on the Buttons)

% --- Executes during object creation, after setting all properties.
function Tpan_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Tpan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    ccm - handles not created until after all CreateFcns called

handles = guidata(gcbo);
load alltools
set(hObject,'CData',Tpan1)
guidata(gcbo,handles)

% --- Executes during object creation, after setting all properties.
function Tpanzoom_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Tpanzoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    ccm - handles not created until after all CreateFcns called
handles = guidata(gcbo);
load alltools
set(hObject,'CData',Tpanzoom)
guidata(gcbo,handles)

% --- Executes during object creation, after setting all properties.
function Tdatacursor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Tdatacursor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    ccm - handles not created until after all CreateFcns called
handles = guidata(gcbo);
load alltools
set(hObject,'CData',Tdatacursor)
guidata(gcbo,handles)

% --- Executes during object creation, after setting all properties.
function Tcapture_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Tcapture (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    ccm - handles not created until after all CreateFcns called
handles = guidata(gcbo);
load alltools
set(hObject,'CData',Tcapture)
guidata(gcbo,handles)

% --- Executes during object creation, after setting all properties.
function Ttable_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ttable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    ccm - handles not created until after all CreateFcns called
handles = guidata(gcbo);
load alltools
set(hObject,'CData',Ttable)
guidata(gcbo,handles)


% --- Executes during object creation, after setting all properties.
function Tmovie_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Tmovie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    ccm - handles not created until after all CreateFcns called
handles = guidata(gcbo);
load alltools
set(hObject,'CData',Tmovie)
guidata(gcbo,handles)

% --- Executes during object creation, after setting all properties.
function Tzoom_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Tzoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    ccm - handles not created until after all CreateFcns called
handles = guidata(gcbo);
load alltools
set(hObject,'CData',Tzoom)
guidata(gcbo,handles)

% --------------------------------------------------------------------
function Empty_Callback(hObject, eventdata, handles)
% hObject    handle to ccm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in Tzoom.
function Tzoom_Callback(hObject, eventdata, handles)
% hObject    handle to Tzoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pan off
%panzoom off
zoom;

% --- Executes on button press in Tpan.
function Tpan_Callback(hObject, eventdata, handles)
% hObject    handle to Tpan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
zoom off
%panzoom off
pan;

% --- Executes on button press in Tpanzoom.
function Tpanzoom_Callback(hObject, eventdata, handles)
% hObject    handle to Tpanzoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pan off
zoom off
%panzoom;

% --- Executes on button press in Tdatacursor.
function Tdatacursor_Callback(hObject, eventdata, handles)
% hObject    handle to Tdatacursor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
datacursormode;
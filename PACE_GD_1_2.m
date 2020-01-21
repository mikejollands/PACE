%% PACE GD - PACE extension to get dimensions of beam 

function varargout = PACE_GD_1_2(varargin)
% PACE_GD_1_2 MATLAB code for PACE_GD_1_2.fig
%      PACE_GD_1_2, by itself, creates a new PACE_GD_1_2 or raises the existing
%      singleton*.
%
%      H = PACE_GD_1_2 returns the handle to a new PACE_GD_1_2 or the handle to
%      the existing singleton*.
%
%      PACE_GD_1_2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PACE_GD_1_2.M with the given input arguments.
%
%      PACE_GD_1_2('Property','Value',...) creates a new PACE_GD_1_2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PACE_GD_1_2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PACE_GD_1_2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PACE_GD_1_2

% Last Modified by GUIDE v2.5 10-Jan-2020 20:51:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PACE_GD_1_2_OpeningFcn, ...
                   'gui_OutputFcn',  @PACE_GD_1_2_OutputFcn, ...
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
end

% --- Executes just before PACE_GD_1_2 is made visible.
function PACE_GD_1_2_OpeningFcn(hObject, ~, handles, varargin)
set( handles.figure1, 'Units', 'pixels' );
screenSize = get(0, 'ScreenSize');
position = get( handles.figure1, 'Position' );
position(3) = 490;
position(4) = 430;
position(1) = (screenSize(3)-position(3))/2; % put the GUI in the middle of the screen
position(2) = (screenSize(4)-position(4))/2; % put the GUI in the middle of the screen
set( handles.figure1,'Position', position );
ax1=handles.axes1;
cla(ax1,'reset');
clear global
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PACE_GD_1_2 (see VARARGIN)

% Choose default command line output for PACE_GD_1_2
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
end
% UIWAIT makes PACE_GD_1_2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = PACE_GD_1_2_OutputFcn(~, ~, handles) 
set(handles.pushbutton2,'Enable','off')
set(handles.pushbutton3,'Enable','off')
set(handles.pushbutton4,'Enable','off')
set(handles.pushbutton5,'Enable','off')
set(handles.pushbutton6,'Enable','off')
set(handles.pushbutton7,'Enable','off')
set(handles.pushbutton9,'Enable','off')


% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
end

%% data import
function pushbutton1_Callback(~, ~, handles) %#ok<DEFNU>

clear global
global ax1
global profile
global profile_length
global s

ax1=handles.axes1;
cla(ax1,'reset');
original_directory=cd;
[file,path] = uigetfile({'*.txt;*.csv;*.xls;*.xlsx;*.dat'});
cd (path);
data=readmatrix(file);
if size(data,2)>2
    prompt={'Distance column (A,B,C etc)','Concentration column (A,B,C etc)'};
    dlgtitle='Import extra input';
    dims=[1 35];
    definput={'',''};
    answer=inputdlg(prompt,dlgtitle,dims,definput);
    data_n(:,1)=readmatrix(file,'Range',strcat(answer{1},':',answer{1}));
    data_n(:,2)=readmatrix(file,'Range',strcat(answer{2},':',answer{2}));
    data=data_n;
end
cd(original_directory);
profile_length=data(:,1);
profile=data(:,2);
s=scatter(ax1,data(:,1),data(:,2), 'filled','MarkerEdgeColor','black','MarkerFaceColor',[0.5 0.5 0.5]);
ylim(ax1,[min([min(profile)*0.8 0]) max(profile)*1.05]);
hold(ax1,'on');
xlim(ax1,[min(data(:,1)) max(data(:,1))])
xlabel(ax1,'Distance')
ylabel(ax1,'Concentration')
set(handles.pushbutton2,'Enable','on')
set(handles.pushbutton3,'Enable','on')
set(handles.pushbutton4,'Enable','on')
set(handles.pushbutton5,'Enable','on')
set(handles.pushbutton6,'Enable','on')
set(handles.pushbutton9,'Enable','off')
end


% circle
function pushbutton2_Callback(~, ~, handles) %#ok<DEFNU>
global shape_type
shape_type=1;
set(handles.pushbutton2,'Enable','on')
set(handles.pushbutton3,'Enable','off')
set(handles.pushbutton4,'Enable','off')
set(handles.pushbutton5,'Enable','off')
set(handles.pushbutton6,'Enable','off')
set(handles.pushbutton7,'Enable','on')
set(handles.pushbutton9,'Enable','on')
end

% square
function pushbutton3_Callback(~, ~, handles) %#ok<DEFNU>
global shape_type
shape_type=2;
set(handles.pushbutton2,'Enable','off')
set(handles.pushbutton3,'Enable','on')
set(handles.pushbutton4,'Enable','off')
set(handles.pushbutton5,'Enable','off')
set(handles.pushbutton6,'Enable','off')
set(handles.pushbutton7,'Enable','on')
set(handles.pushbutton9,'Enable','on')
end

% gauss
function pushbutton4_Callback(~, ~, handles) %#ok<DEFNU>
global shape_type
shape_type=3;
set(handles.pushbutton2,'Enable','off')
set(handles.pushbutton3,'Enable','off')
set(handles.pushbutton4,'Enable','on')
set(handles.pushbutton5,'Enable','off')
set(handles.pushbutton6,'Enable','off')
set(handles.pushbutton7,'Enable','on')
set(handles.pushbutton9,'Enable','on')
end

% lorentz
function pushbutton5_Callback(~, ~, handles) %#ok<DEFNU>
global shape_type
shape_type=4;
set(handles.pushbutton2,'Enable','off')
set(handles.pushbutton3,'Enable','off')
set(handles.pushbutton4,'Enable','off')
set(handles.pushbutton5,'Enable','on')
set(handles.pushbutton6,'Enable','off')
set(handles.pushbutton7,'Enable','on')
set(handles.pushbutton9,'Enable','on')
end

% voigt
function pushbutton6_Callback(~, ~, handles) %#ok<DEFNU>
global shape_type
shape_type=5;
set(handles.pushbutton2,'Enable','off')
set(handles.pushbutton3,'Enable','off')
set(handles.pushbutton4,'Enable','off')
set(handles.pushbutton5,'Enable','off')
set(handles.pushbutton6,'Enable','on')
set(handles.pushbutton7,'Enable','on')
set(handles.pushbutton9,'Enable','on')
end

% reset curve choice
function pushbutton9_Callback(~, ~, handles) %#ok<DEFNU>
global profile
global profile_length
global ax1
global s
cla(ax1,'reset')
s=scatter(ax1,profile_length,profile, 'filled','MarkerEdgeColor','black','MarkerFaceColor',[0.5 0.5 0.5]);
xlabel(ax1,'Distance')
ylabel(ax1,'Concentration')
ylim(ax1,[min([min(profile)*0.8 0]) max(profile)*1.05]);
xlim(ax1,[min(profile_length) max(profile_length)])
hold(ax1,'on')
set(handles.pushbutton2,'Enable','on')
set(handles.pushbutton3,'Enable','on')
set(handles.pushbutton4,'Enable','on')
set(handles.pushbutton5,'Enable','on')
set(handles.pushbutton6,'Enable','on')
set(handles.pushbutton7,'Enable','on')
set(handles.pushbutton9,'Enable','off')
end

% fit
function pushbutton7_Callback(~, ~, ~) %#ok<DEFNU>
global profile
global profile_length
global shape_type
global width
global ax1
global l1
global s
warning off
cla(ax1,'reset')
s=scatter(ax1,profile_length,profile, 'filled','MarkerEdgeColor','black','MarkerFaceColor',[0.5 0.5 0.5]);
xlabel(ax1,'Distance')
ylabel(ax1,'Concentration')
ylim(ax1,[min([min(profile)*0.8 0]) max(profile)*1.05]);
xlim(ax1,[min(profile_length) max(profile_length)])
hold(ax1,'on')

if shape_type==1
fun=@(x,profile_length) x(2)+(x(1)-x(2))*0.5*erfc((profile_length-x(3))/(x(4)));
MSE_min=1e20;
for k=[0 10 20 50 100 200 500 1000 2000 5000]
    if mean(profile(1:10))>mean(profile(max(size(profile))-10:end))
        x0 = [max(profile),min(profile),mean(profile_length),0+k];
    else
        x0 = [min(profile),max(profile),mean(profile_length),0+k];
    end
    [x1,~,~,~,MSE]=nlinfit(profile_length,profile,fun,x0);
    if MSE<MSE_min
        MSE_min=MSE;
        x_min=x1;
    end
end
x1=x_min;
c1=x1(1);
c2=x1(2);
mp=x1(3);
width=2*x1(4);
dx=min([width/10 ((max(profile_length))-min(profile_length))/(max([25 4*(max(size(profile_length)))]))]);
model_x=min(profile_length):dx:max(profile_length);
options=optimset('Display','off','DiffMinChange',1);
fun=@(x) f_convolution_circ(x(1),x(2),x(3),x(4),model_x,profile_length,dx)-profile;
x0=[c1 c2 mp width];
x = lsqnonlin(fun,x0,[],[],options);
width=x(4);
[~,convoluted,model_x]=f_convolution_circ(x(1),x(2),x(3),x(4),model_x,profile_length,dx);
l1=plot(ax1,model_x,convoluted,'red','Linewidth',2);
hold on 
t=text(ax1,((max(profile_length)-min(profile_length))*0.7+min(profile_length)),(((max(profile)-min(profile))*0.5+min(profile))),sprintf('%s%s', 'Diameter = ',num2str(round(width,3,'significant'))));
end

if shape_type==2
fun=@(x) f_convolution_sq(x(1),x(2),x(3),profile_length,x(4))-profile;
options=optimset('Display','off','DiffMinChange',1);
rsn_min=1e100;
for k=[1 20 50 100 200 500 1000 2000 5000]
    if mean(profile(1:10))>mean(profile(max(size(profile))-10:end))
        x0 = [max(profile),min(profile),mean(profile_length),k];
    else
        x0 = [min(profile),max(profile),mean(profile_length),k];
    end
    [x,rsn]=lsqnonlin(fun,x0,[],[],options);
    if rsn<rsn_min
        rsn_min=rsn;
        x_min=x;
    end
end
x=x_min;
width=x(4);
[~,convoluted,model_x]=f_convolution_sq(x(1),x(2),x(3),profile_length,x(4));
[~,index1]=min(abs(model_x-min(profile_length)));
[~,index2]=min(abs(model_x-max(profile_length)));
model_x=model_x(index1:index2);
convoluted=convoluted(index1:index2);
l1=plot(ax1,model_x,convoluted,'red','Linewidth',2);
t=text(ax1,((max(profile_length)-min(profile_length))*0.7+min(profile_length)),(((max(profile)-min(profile))*0.5+min(profile))),sprintf('%s%s', 'Width = ',num2str(round(width,3,'significant'))));
end

if shape_type==3   
options=optimset('Display','off','DiffMinChange',1);
fun=@(x,profile_length) x(2)+(x(1)-x(2))*0.5*erfc((profile_length-x(3))/(sqrt(2)*x(4)));
MSE_min=1e20;
for k=[0.001 0.1 1 10 20 50 100 200 500 1000 2000 5000]
    if mean(profile(1:10))>mean(profile(max(size(profile))-10:end))
        x0 = [max(profile),min(profile),mean(profile_length),k];
    else
        x0 = [min(profile),max(profile),mean(profile_length),k];
    end
    [x,~,~,~,MSE]=nlinfit(profile_length,profile,fun,x0,options);
    if MSE<MSE_min
        MSE_min=MSE;
        x_min=x;
    end
end
x=x_min;
sigma=x(4);
FWHM=sigma*2*sqrt(2*log(2));
dx=min([FWHM/10 ((max(profile_length))-min(profile_length))/(max([10 (max(size(profile_length)))]))]);
model_x=min(profile_length):dx:max(profile_length);
model=x(2)+(x(1)-x(2))*0.5*erfc((model_x-x(3))/(sqrt(2)*x(4)));
l1=plot(model_x,model,'red','Linewidth',2);
t=text(ax1,((max(profile_length)-min(profile_length))*0.7+min(profile_length)),(((max(profile)-min(profile))*0.5+min(profile))),sprintf('%s%s\n%s%s', 'sigma = ',num2str(round(sigma,3,'significant')), 'FWHM = ',num2str(round(FWHM,3,'significant'))));
end

if shape_type==4  
    options=optimset('Display','off','DiffMinChange',1);

fun=@(x,profile_length) x(2)+(x(1)-x(2))*(1/(pi)*atan((profile_length-x(3))/x(4))+0.5);
MSE_min=1e20;
for k=[0.001 0.1 1 10 20 50 100 200 500 1000 2000 5000]
    if mean(profile(1:10))>mean(profile(max(size(profile))-10:end))
        x0 = [min(profile),max(profile),mean(profile_length),k];
    else
        x0 = [max(profile),min(profile),mean(profile_length),k];
    end
    [x,~,~,~,MSE]=nlinfit(profile_length,profile,fun,x0,options);
    if MSE<MSE_min
        MSE_min=MSE;
        x_min=x;
    end
end
x=x_min;
FWHM=x(4)*2;
% FWHM=sigma*2*sqrt(2*log(2));
dx=min([FWHM/10 ((max(profile_length))-min(profile_length))/(max([10 (max(size(profile_length)))]))]);
model_x=min(profile_length):dx:max(profile_length);
model=x(2)+(x(1)-x(2))*(1/(pi)*atan((model_x-x(3))/x(4))+0.5);
plot(model_x,model,'red','Linewidth',2)
t=text(ax1,((max(profile_length)-min(profile_length))*0.7+min(profile_length)),(((max(profile)-min(profile))*0.5+min(profile))),sprintf('%s%s','FWHM = ',num2str(round(FWHM,3,'significant'))));

end

if shape_type==5
    options=optimset('Display','off','DiffMinChange',1);
fun=@(x,profile_length) (x(2)+(x(1)-x(2))*((1.36603*(x(4)/(0.5346*x(4)+sqrt(0.2166*x(4)^2+x(5)^2)))-0.47719*(x(4)/(0.5346*x(4)+sqrt(0.2166*x(4)^2+x(5)^2)))^2+0.11116*(x(4)/(0.5346*x(4)+sqrt(0.2166*x(4)^2+x(5)^2)))^3)*((1/(pi)*atan((profile_length-x(3))/(x(4)/2))+0.5))+(1-(1.36603*(x(4)/(0.5346*x(4)+sqrt(0.2166*x(4)^2+x(5)^2)))-0.47719*(x(4)/(0.5346*x(4)+sqrt(0.2166*x(4)^2+x(5)^2)))^2+0.11116*(x(4)/(0.5346*x(4)+sqrt(0.2166*x(4)^2+x(5)^2)))^3))*((1+erf((profile_length-x(3))/(sqrt(2)*(x(5)/2.355))))*0.5)));
MSE_min=1e20;
for k=[0.001 0.1 1 10 20 50 100 200 500 1000 2000 5000]
    if mean(profile(1:10))>mean(profile(max(size(profile))-10:end))
        x0 = [min(profile),max(profile),mean(profile_length),k,k];
    else
        x0 = [max(profile),min(profile),mean(profile_length),k,k];
    end
    [x,~,~,~,MSE]=nlinfit(profile_length,profile,fun,x0);
    if MSE<MSE_min
        MSE_min=MSE;
        x_min=x;
    end
end
x=x_min;  
FWHM_L=x(4);
FWHM_G=x(5);
FWHM=(FWHM_G^5+2.69269*FWHM_G^4*FWHM_L+2.42843*FWHM_G^3*FWHM_L^2+4.47163*FWHM_G^2*FWHM_L^3+0.07842*FWHM_G*FWHM_L^4+FWHM_L^5)^(1/5);
dx=min([max([FWHM_L FWHM_G])/10 ((max(profile_length))-min(profile_length))/(max([10 (max(size(profile_length)))]))]);
model_x=min(profile_length):dx:max(profile_length);
model=(x(2)+(x(1)-x(2))*((1.36603*(x(4)/(0.5346*x(4)+sqrt(0.2166*x(4)^2+x(5)^2)))-0.47719*(x(4)/(0.5346*x(4)+sqrt(0.2166*x(4)^2+x(5)^2)))^2+0.11116*(x(4)/(0.5346*x(4)+sqrt(0.2166*x(4)^2+x(5)^2)))^3)*((1/(pi)*atan((model_x-x(3))/(x(4)/2))+0.5))+(1-(1.36603*(x(4)/(0.5346*x(4)+sqrt(0.2166*x(4)^2+x(5)^2)))-0.47719*(x(4)/(0.5346*x(4)+sqrt(0.2166*x(4)^2+x(5)^2)))^2+0.11116*(x(4)/(0.5346*x(4)+sqrt(0.2166*x(4)^2+x(5)^2)))^3))*((1+erf((model_x-x(3))/(sqrt(2)*(x(5)/2.355))))*0.5)));
plot(model_x,model,'red','Linewidth',2)       
t=text(ax1,((max(profile_length)-min(profile_length))*0.6+min(profile_length)),(((max(profile)-min(profile))*0.5+min(profile))),sprintf('%s%s\n%s%s\n%s%s', 'FWHM(Gauss) = ',num2str(round(FWHM_G,3,'significant')), 'FWHM(Lorentz) = ',num2str(round(FWHM_L,3,'significant')),'FWHM(Voigt) = ',num2str(round(FWHM,3,'significant'))));
end
t.FontUnits='pixels';
t.FontSize=10;
end

% square fitting function
function [convoluted_downsampled_sq,convoluted,model_x]=f_convolution_sq(c1,c2,mp,profile_length,width)
dx=min([width/10 ((max(profile_length))-min(profile_length))/(max([10 (max(size(profile_length)))]))]);
model_x=min(profile_length)-(max(profile_length)-min(profile_length))/2:dx:max(profile_length)+(max(profile_length)-min(profile_length))/2;
convoluted=movmean((c2+(c1-c2)*0.5*erfc(((model_x-mp))/(1e-50))),[floor((round(width/mean(diff(model_x)))-1)/2) ceil((round(width/mean(diff(model_x)))-1)/2)]);
convoluted_downsampled_sq=interp1(model_x,convoluted,profile_length,'nearest');
end

function [convoluted_downsampled,convoluted,model_x]=f_convolution_circ(c1,c2,mp,width,model_x,profile_length,dx)
original=(c2+(c1-c2)*0.5*erfc(((model_x-mp))/(1e-100)));
[X,Y] = meshgrid(0:dx:width,0:dx:width);
centre=[mean(mean(X,1)) mean(mean(Y,2))];
dist=sqrt(abs(X-centre(1)).^2+abs(Y-centre(2)).^2);
dist(dist>width/2)=NaN;
dist(~isnan(dist))=1;
dist(isnan(dist))=0;
proportions=sum(dist)./(sum(sum(dist)));
proportions=[zeros(1,floor((size(model_x,2)-size(proportions,2))/2)) proportions zeros(1,ceil((size(model_x,2)-size(proportions,2))/2))];
props_matrix=zeros(size(proportions,2)+400,size(model_x,2));
for k=1:size(props_matrix,2)
    props_matrix(:,(k))=circshift([zeros(1,200) proportions zeros(1,200)]',(k-round(size(props_matrix,2)/2)));
end
props_matrix(1:200,:)=[];
props_matrix=flipud(props_matrix);
props_matrix(1:200,:)=[];
props_matrix=flipud(props_matrix);
pmsum=1./(sum(props_matrix));
convoluted=original.*props_matrix';
convoluted=(sum(convoluted')); %#ok<*UDIM>
convoluted=convoluted.*pmsum;
convoluted_downsampled=interp1(model_x,convoluted,profile_length,'nearest');
end

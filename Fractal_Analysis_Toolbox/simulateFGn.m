function varargout = simulateFGn(varargin)
% SIMULATEFGN MATLAB code for simulateFGn.fig
%      SIMULATEFGN, by itself, creates a new SIMULATEFGN or raises the existing
%      singleton*.
%
%      H = SIMULATEFGN returns the handle to a new SIMULATEFGN or the handle to
%      the existing singleton*.
%
%      SIMULATEFGN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SIMULATEFGN.M with the given input arguments.
%
%      SIMULATEFGN('Property','Value',...) creates a new SIMULATEFGN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before simulateFGn_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to simulateFGn_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% Last Modified by Michael J. Richardson v2.5 29-May-2015 14:42:03

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @simulateFGn_OpeningFcn, ...
                   'gui_OutputFcn',  @simulateFGn_OutputFcn, ...
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

% --- Executes just before simulateFGn is made visible.
function simulateFGn_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to simulateFGn (see VARARGIN)

% Choose default command line output for simulateFGn
handles.output = hObject;
handles.fdata = 0;
handles.FBnDo = 0;

% Update handles structure
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = simulateFGn_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in pbSimulate.
function pbSimulate_Callback(hObject, eventdata, handles) %#ok<*DEFNU>
% hObject    handle to pbSimulate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
H = str2double(get(handles.edParamH,'String'));
datalength = str2double(get(handles.edParamLength,'String'));

if handles.FBnDo == 1
    fdata = simulate_fGn(H, datalength);
    fdata = cumsum(fdata);
else
    fdata = simulate_fGn(H, datalength);
end

handles.fdata = fdata;
guidata(hObject, handles);

set(handles.pbOutput,'Enable','On');

%% plot simulated data
axes(handles.axTS);
plot(1:length(fdata), fdata,'b-');
xlabel('Time, t');
ylabel('Amplitude, Y(t)');
xlim([1 length(fdata)]);


% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.figure1)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
                     ['Close ' get(handles.figure1,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.figure1)


% --- Executes on button press in pbOutput.
function pbOutput_Callback(hObject, eventdata, handles)
% hObject    handle to pbOutput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% File Saving prompts here:
[OutFileName,OutPathName] = uiputfile('*.txt','Specify an output file name');
if isequal(OutFileName,0) || isequal(OutPathName,0)
   %quit if no outfile is specified
    h=warndlg('No file created; Unable to save data', 'File Error!!'); uiwait(h);
    return;
end

%% Save data
fid=fopen([OutPathName OutFileName],'wt');
fprintf(fid,'%1.4f\r\n',handles.fdata);
fclose(fid);



function edParamH_Callback(hObject, eventdata, handles) %#ok<*INUSD>
% hObject    handle to edParamH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
temp = str2double(get(handles.edParamH,'String'));
if temp < .01
   temp = .01;
   set(handles.edParamH,'String', num2str(temp));
end

if temp > 0.99
    temp = 0.99;
    set(handles.edParamH,'String', num2str(temp));
end

% --- Executes during object creation, after setting all properties.
function edParamH_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edParamH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edParamLength_Callback(hObject, eventdata, handles)
% hObject    handle to edParamLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
temp = str2double(get(handles.edParamLength,'String'));
temp = round(temp);
set(handles.edParamLength,'String', num2str(temp));
if temp < 32
   temp = 32;
   set(handles.edParamLength,'String', num2str(temp));
end

if mod(log2(temp),1)~=0
     set(handles.edParamLength,'String', '1024');
end

% --- Executes during object creation, after setting all properties.
function edParamLength_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edParamLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes when selected object is changed in uipanel8.
function uipanel8_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel8 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
    case 'rbFGn'
        handles.FBnDo = 0;
        if (length(handles.fdata) > 2)
            fdata = diff(handles.fdata);
            handles.fdata = fdata;
            %% plot simulated data
            axes(handles.axTS);
            plot(1:length(fdata), fdata,'b-');
            xlabel('Time, t');
            ylabel('Amplitude, Y(t)');
            xlim([1 length(fdata)]);
        end
    case 'rbFBn'
        handles.FBnDo = 1;
        if (length(handles.fdata) > 2)
            fdata = cumsum(handles.fdata);
            handles.fdata = fdata;
            %% plot simulated data
            axes(handles.axTS);
            plot(1:length(fdata), fdata,'b-');
            xlabel('Time, t');
            ylabel('Amplitude, Y(t)');
            xlim([1 length(fdata)]);
        end
end
guidata(hObject, handles);


function fgn = simulate_fGn(H, lng)
% This function simulates fGn time series using Davies and Harte's (1997) algorithm. Steps
% in the algorighm were presented in Caccia et al. (1996). The
% implementation of steps 3 and 4 is borrowed from the ffgn function by 
% Yingchun Zhou (Jasmine), zhouyc@math.bu.edu and Stilian Stoev, sstoev@umich.edu)
    
N=lng;
if mod(log2(N),1)~=0 %Integer test, make sure no decimal remainder after dividing by 1
   %disp('Series must be an integer power of 2 in length.')
   N=2^nextpow2(N);
   %disp(['N set to ' num2str(N) ])
end

% N = 2^13;
M = 2*N;
n = 1;

sigma = 1;
tau = 0:M/2;
% H = .15;

%% Specify desired autocorrelation
gamma = (sigma^2/2).*(abs(tau+1).^(2*H)-2.*abs(tau).^(2*H)+abs(tau-1).^(2*H));
g = [gamma(1:end) fliplr(gamma(2:end-1))];

%% Step 1 (Calculate exact spectral power expected for the theoretical autocorrelation function)
S = real(fft(g));

%% Step 2 (Check if all spectral coefficients are greater than 0)
if min(S)<0,
    h = warndlg(' Some of the Sj are negative!', 'Analaysis Error'); uiwait(h);
    return;
end
S = abs(S);

%% Step 3 (Calculate randomized spectral amplitudes) 
z(:,1)=sqrt(2)*randn(n,1);
y(:,1)=z(:,1);

z(:,N+1)=sqrt(2)*randn(n,1);
y(:,N+1)=z(:,N+1);

a=randn(n,N-1);
b=randn(n,N-1);

z1=a+b*1i;
z(:,2:N)=z1; %#ok<NASGU>
y1=z1;

y(:,2:N)=y1;
y(:,N+2:2*N)=conj(y(:,N:-1:2));

y = y.*(ones(n,1)*sqrt(S));

%% Step 4 (use
fgn = real(fft(y')')/sqrt(4*N);
fgn = fgn(:,1:N);

% Random number generator for H = 0.5
if H == .5
    fgn = randn(N,1);
end

return

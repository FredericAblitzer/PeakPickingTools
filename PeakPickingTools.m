function varargout = PeakPickingTools(varargin)
% PeakPickingTools
% v 1.0 (April 14, 2022)
% 
% PeakPickingTools is a MATLAB toolbox to extract modal parameters from frequency response functions, using a robust peak-picking identification technique.
%
% It is an implementation of the method described in:
% Ablitzer F. 2021. Peak-picking identification technique for modal expansion of input impedance of brass instruments. Acta Acustica, 5, 53.
% https://doi.org/10.1051/aacus/2021046
%
%
% Copyright (c) 2022 Frédéric Ablitzer
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

% Last Modified by GUIDE v2.5 14-Apr-2022 13:52:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @PeakPickingTools_OpeningFcn, ...
    'gui_OutputFcn',  @PeakPickingTools_OutputFcn, ...
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


% --- Executes just before PeakPickingTools is made visible.
function PeakPickingTools_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PeakPickingTools (see VARARGIN)

% Choose default command line output for PeakPickingTools
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


% UIWAIT makes PeakPickingTools wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% add to path the folder containing user functions
addpath('dev')
% use f structure to make functions accessibles from outside
global f
f.update_plot = @update_plot;
f.update_list = @update_list;
f.FindAmplitudes = @FindAmplitudes;
f.Zmodal = @Zmodal;
f.Init = @Init;
f.ModifyFRF = @ModifyFRF;
f.AddMode = @AddMode;
f.RemoveMode = @RemoveMode;
f.RefineIdentification = @RefineIdentification;
f.ChangeFittingBandwidth = @ChangeFittingBandwidth;



% --- Outputs from this function are returned to the command line.
function varargout = PeakPickingTools_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in listboxModes.
function listboxModes_Callback(hObject, eventdata, handles)
% hObject    handle to listboxModes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listboxModes contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listboxModes

global s

update_list(handles)
update_plot(handles)


% --- Executes during object creation, after setting all properties.
function listboxModes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listboxModes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonAddMode.
function pushbuttonAddMode_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonAddMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global s

while true
    
    set(handles.textStatus,'String','Pick a peak... (right clic to stop)','Visible','on')
    [freq_0,~,button] = ginput(1);
    
    if button ~= 1
        set(handles.textStatus,'String','','Visible','off')
        break
    end
    
    AddMode(handles,freq_0,s.xi0);
    
end


% --- Executes on button press in pushbuttonRemove.
function pushbuttonRemove_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonRemove (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
n = get(handles.listboxModes,'Value');
RemoveMode(handles,n);

% --- Executes on button press in pushbuttonRefine.
function pushbuttonRefine_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonRefine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
RefineIdentification(handles)


% --- Executes on button press in radiobuttonModesR.
function radiobuttonModesR_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonModesR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobuttonModesR
if get(handles.radiobuttonModesC,'Value')
    SwitchRealComplex(handles);
end
set(handles.radiobuttonModesR,'Value',1)
set(handles.radiobuttonModesC,'Value',0)


% --- Executes on button press in radiobuttonModesC.
function radiobuttonModesC_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonModesC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobuttonModesC
if get(handles.radiobuttonModesR,'Value')
    SwitchRealComplex(handles);
end
set(handles.radiobuttonModesC,'Value',1)
set(handles.radiobuttonModesR,'Value',0)


function editBandwidth_Callback(hObject, eventdata, handles)
% hObject    handle to editBandwidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editBandwidth as text
%        str2double(get(hObject,'String')) returns contents of editBandwidth as a double

delta_freq = str2double(get(hObject,'String'));
ChangeFittingBandwidth(handles,delta_freq);




% --- Executes during object creation, after setting all properties.
function editBandwidth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editBandwidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --------------------------------------------------------------------
function uipushtoolOpen_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtoolOpen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global s



[FileName,PathName] = uigetfile('*.dat;*.txt;*.mat');

if strcmp(FileName(end-2:end),'dat') || strcmp(FileName(end-2:end),'txt')
    
    data = load([PathName FileName]);
    
    freq = data(:,1);
    Ze = data(:,2) + 1i*data(:,3);
    
    freq = freq(:).';
    
    s.freq_original = freq;
    s.Ze_original = Ze;
    
    Init(handles);
    
    ModifyFRF(handles);
    
    
else
    
    load([PathName FileName]);
    
    set(handles.editBandwidth,'String',num2str(s.delta_freq));
    
    set(handles.listboxModes,'Value',1)
    
    if s.modesC == 1
        set(handles.radiobuttonModesR,'Value',0)
        set(handles.radiobuttonModesC,'Value',1)
    else
        set(handles.radiobuttonModesR,'Value',1)
        set(handles.radiobuttonModesC,'Value',0)
    end
    
    ModifyFRF(handles);
    
end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% UPDATE PLOT
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function update_plot(handles)

global s

freq = s.freq;
Ze = s.Ze;
omega = 2*pi*freq;

df = s.freq(2)-s.freq(1);

fffreq = [df:df:s.fmax];
ooomega = 2*pi*fffreq;


COLOR_all = repmat(0.75*[1 1 1],length(s.OMEGA_n),1);
COLOR_s = repmat([0 0.75 0.75],length(s.OMEGA_n),1); % selected
color_Z = [0 0.4471 0.7412];







axes(handles.axesTop)


% autoscale or fixed scale ?
if isempty(s.magmin) || isempty(s.magmax)
    % prepare autoscale
    cla
    plot(freq,abs(Ze),'k','Linewidth',2)
    xlim([s.fmin s.fmax])
    set(gca,'YScale','log')
    ylim auto
    limY = ylim;
else
    limY = [s.magmin s.magmax];
end

cla
hold on

if strcmp(handles.menu_View_Hide_Fitting_Frequency_Domains.Checked,'off')
    
    % plot fitting bandwidth for all modes
    for n=1:length(s.OMEGA_n)
        rectangle('Parent',gca,'Position',[1/(2*pi)*s.OMEGA_n_0(n)-s.delta_freq/2 1e-100 s.delta_freq 1e100],'FaceColor',[0.9 1 0.9],'EdgeColor','none');
    end
    
    % plot fitting bandwidth for selected mode
    n = get(handles.listboxModes,'Value');
    if ~isempty(n)
        rectangle('Parent',gca,'Position',[1/(2*pi)*s.OMEGA_n_0(n)-s.delta_freq/2 1e-100 s.delta_freq 1e100],'FaceColor',[0.8 1 0.8],'EdgeColor','none');
    end
    
end

if strcmp(handles.menu_View_Hide_Modal_Contributions.Checked,'off')
    
    % plot the contribution of selected mode
    n = get(handles.listboxModes,'Value');
    if ~isempty(n)
        if ~s.locked(n)
            plot(fffreq,abs(Zmodal(ooomega,s.OMEGA_n(n),s.XI_n(n),s.A_n(n),s.B_n(n))),'-','Color',COLOR_s(n,:),'Linewidth',3)
        else
            plot(fffreq,abs(Zmodal(ooomega,s.OMEGA_n(n),s.XI_n(n),s.A_n(n),s.B_n(n))),'-','Color',[1 0.5 0],'Linewidth',3)
        end
    end
    
    % plot the contribution of all modes
    for n=1:length(s.OMEGA_n)
        if ~s.locked(n)
            plot(fffreq,abs(Zmodal(ooomega,s.OMEGA_n(n),s.XI_n(n),s.A_n(n),s.B_n(n))),'-','Color',COLOR_all(n,:),'Linewidth',1)
        else
            plot(fffreq,abs(Zmodal(ooomega,s.OMEGA_n(n),s.XI_n(n),s.A_n(n),s.B_n(n))),'-','Color',[1 0.5 0],'Linewidth',1)
        end
    end
    
end

% plot measured impedance
plot(freq,abs(Ze),'Color','k','Linewidth',2)

if strcmp(handles.menu_View_Hide_Reconstructed_FRF.Checked,'off')
    
    % plot reconstructed impedance
    if ~isempty(s.OMEGA_n)
        plot(fffreq,abs(Zmodal(ooomega,s.OMEGA_n,s.XI_n,s.A_n,s.B_n)),'-','Color',color_Z,'Linewidth',2)
    end
    
end

set(gca,'YScale','log')

xlim([s.fmin s.fmax])
ylim(limY);

box on
set(gca,'Layer','top')
xlabel('frequency (Hz)')
ylabel('magnitude')
%set(findall(gca,'-property','FontSize'),'Fontsize',12)

axes(handles.axesBottom)

cla
hold on

if strcmp(handles.menu_View_Hide_Fitting_Frequency_Domains.Checked,'off')
    
    % plot fitting bandwidth for all modes
    for n=1:length(s.OMEGA_n)
        rectangle('Parent',gca,'Position',[1/(2*pi)*s.OMEGA_n_0(n)-s.delta_freq/2 -pi s.delta_freq 2*pi],'FaceColor',[0.9 1 0.9],'EdgeColor','none');
    end
    
    % plot fitting bandwidth for selected mode
    n = get(handles.listboxModes,'Value');
    if ~isempty(n)
        rectangle('Parent',gca,'Position',[1/(2*pi)*s.OMEGA_n_0(n)-s.delta_freq/2 -pi s.delta_freq 2*pi],'FaceColor',[0.8 1 0.8],'EdgeColor','none');
    end
    
end


if strcmp(handles.menu_View_Hide_Modal_Contributions.Checked,'off')
    
    % plot the contribution of selected mode
    n = get(handles.listboxModes,'Value');
    if ~isempty(n)
        if ~s.locked(n)
            plot(fffreq,angle(Zmodal(ooomega,s.OMEGA_n(n),s.XI_n(n),s.A_n(n),s.B_n(n))),'-','Color',COLOR_s(n,:),'Linewidth',3)
        else
            plot(fffreq,angle(Zmodal(ooomega,s.OMEGA_n(n),s.XI_n(n),s.A_n(n),s.B_n(n))),'-','Color',[1 0.5 0],'Linewidth',3)
        end
    end
    
    % plot the contribution of all modes
    for n=1:length(s.OMEGA_n)
        if ~s.locked(n)
            plot(fffreq,angle(Zmodal(ooomega,s.OMEGA_n(n),s.XI_n(n),s.A_n(n),s.B_n(n))),'-','Color',COLOR_all(n,:),'Linewidth',1)
        else
            plot(fffreq,angle(Zmodal(ooomega,s.OMEGA_n(n),s.XI_n(n),s.A_n(n),s.B_n(n))),'-','Color',[1 0.5 0],'Linewidth',1)
        end
    end
    
end

% plot measured impedance
plot(freq,angle(Ze),'Color','k','Linewidth',2)

if strcmp(handles.menu_View_Hide_Reconstructed_FRF.Checked,'off')
    
    % plot reconstructed impedance
    if ~isempty(s.OMEGA_n)
        plot(fffreq,angle(Zmodal(ooomega,s.OMEGA_n,s.XI_n,s.A_n,s.B_n)),'-','Color',color_Z,'Linewidth',2)
    end
    
end

xlim([s.fmin s.fmax])
ylim([-1 1]*1.5*pi/2);

box on
set(gca,'Layer','top')
xlabel('frequency (Hz)')
ylabel('phase (rad)')

linkaxes([handles.axesTop,handles.axesBottom],'x')

%set(findall(gca,'-property','FontSize'),'Fontsize',12)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% UPDATE LIST
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function update_list(handles)

global s

if ~isempty(s.OMEGA_n)
    
    [~,i] = sort(s.OMEGA_n);
    
    s.OMEGA_n_0 = s.OMEGA_n_0(i);
    
    
    s.select = (sum(abs(ones(size(s.OMEGA_n_0))*s.freq-(s.OMEGA_n_0/(2*pi))*ones(size(s.freq))) < s.delta_freq/2,1))>0;
    
    s.locked = s.locked(i);
    
    s.OMEGA_n = s.OMEGA_n(i);
    s.XI_n = s.XI_n(i);
    s.A_n = s.A_n(i);
    s.B_n = s.B_n(i);
    
    n_prev = get(handles.listboxModes,'Value');
    N_prev = length(get(handles.listboxModes,'String'));
    
    liste = {};
    
    for n=1:length(s.OMEGA_n)
        if ~s.locked(n)
            if abs((s.OMEGA_n(n)-s.OMEGA_n_0(n))/(2*pi))<s.delta_freq/2
                liste{n} = sprintf('    mode %2.0i  %7.2f Hz  %5.2f %%',n,s.OMEGA_n(n)/(2*pi),s.XI_n(n)*100);
            else
                liste{n} = sprintf('(s) mode %2.0i  %7.2f Hz  %5.2f %%',n,s.OMEGA_n(n)/(2*pi),s.XI_n(n)*100);
            end
        else
            liste{n} = sprintf('(L) mode %2.0i  %7.2f Hz  %5.2f %%',n,s.OMEGA_n(n)/(2*pi),s.XI_n(n)*100);
        end
    end
    
    
    set(handles.listboxModes,'FontName','FixedWidth','FontWeight','bold','FontSize',12)
    set(handles.listboxModes,'String',liste)
    
    if length(s.OMEGA_n) > N_prev
        set(handles.listboxModes,'Value',find(i==length(s.OMEGA_n)));
    else
        set(handles.listboxModes,'Value',find(i==n_prev));
    end
    
    n = get(handles.listboxModes,'Value');
    
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INITIALIZATION
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Init(handles)

global s

% keep history of successive actions (up to 100, modify if needed)
s.iter = 0;
s.histACTION = {};
s.histRES = NaN*zeros(1,100);
s.histOMEGA_n = NaN*zeros(100,100);
s.histOMEGA_n_0 = NaN*zeros(100,100);
s.histXI_n = NaN*zeros(100,100);
s.histA_n = NaN*zeros(100,100);
s.histB_n = NaN*zeros(100,100);
s.histdelta_freq = NaN*zeros(1,100);


s.integrate = 0;
s.differentiate = 0;
s.changesign = 0;
s.invert = 0;
s.step = 1;


s.OMEGA_n_0 = [];
%s.select = freq==freq;

s.locked = [];

s.OMEGA_n = [];
s.XI_n = [];
s.A_n = [];
s.B_n = [];

s.fmin = 0;
s.fmax = min(s.freq_original(end),1000);

s.xi0 = 0.05;

s.magmin = [];
s.magmax = [];



s.delta_freq = 50;
set(handles.editBandwidth,'String',num2str(s.delta_freq));

s.modesC = 0;
set(handles.radiobuttonModesC,'Value',s.modesC)
set(handles.radiobuttonModesR,'Value',~s.modesC)


set(handles.listboxModes,'String',{})
set(handles.listboxModes,'Value',[])





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ADDING A MODE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function AddMode(handles,freq_0,xi_0)

global s

s.OMEGA_n_0 = [s.OMEGA_n_0 ; 2*pi*freq_0];

s.select = (sum(abs(ones(size(s.OMEGA_n_0))*s.freq-(s.OMEGA_n_0/(2*pi))*ones(size(s.freq))) < s.delta_freq/2,1))>0;

s.OMEGA_n = [s.OMEGA_n ; 2*pi*freq_0];
s.XI_n = [s.XI_n ; xi_0];

s.locked = [s.locked ; 0];

n = length(s.OMEGA_n);


X0 = [s.OMEGA_n(n) s.XI_n(n)];

Xsol = fminsearch(@(X) FindAmplitudes(...
    [ s.OMEGA_n([1:n-1]) ; X(1) ; s.OMEGA_n([n+1:end])],...
    [ s.XI_n([1:n-1]) ; X(2) ; s.XI_n([n+1:end])]),X0);%,optimset('Display','iter'))
omega_n = Xsol(1);
xi_n = Xsol(2);


s.OMEGA_n(n) = omega_n;
s.XI_n(n) = xi_n;


[RES,A_n,B_n] = FindAmplitudes(s.OMEGA_n,s.XI_n);
s.A_n = A_n;
s.B_n = B_n;

update_list(handles)
update_plot(handles)

s.iter = s.iter+1;
s.histRES(s.iter) = RES;
s.histACTION{s.iter} = 'A';
s.histOMEGA_n(1:length(s.OMEGA_n),s.iter) = s.OMEGA_n;
s.histXI_n(1:length(s.OMEGA_n),s.iter) = s.XI_n;
s.histA_n(1:length(s.OMEGA_n),s.iter) = s.A_n;
s.histB_n(1:length(s.OMEGA_n),s.iter) = s.B_n;
s.histOMEGA_n_0(1:length(s.OMEGA_n),s.iter) = s.OMEGA_n_0;
s.histdelta_freq(s.iter) = s.delta_freq;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% REMOVING A MODE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function RemoveMode(handles,n)

global s

if n == length(s.OMEGA_n)
    if n>1
        set(handles.listboxModes,'Value',n-1);
    else
        set(handles.listboxModes,'String',{})
        set(handles.listboxModes,'Value',[])
    end
end


s.OMEGA_n_0(n) = [];

if ~isempty(s.OMEGA_n_0)
    s.select = (sum(abs(ones(size(s.OMEGA_n_0))*s.freq-(s.OMEGA_n_0/(2*pi))*ones(size(s.freq))) < s.delta_freq/2,1))>0;
end

s.locked(n) = [];

s.OMEGA_n(n) = [];
s.XI_n(n) = [];
s.A_n(n) = [];
s.B_n(n) = [];



[RES,A_n,B_n] = FindAmplitudes(s.OMEGA_n,s.XI_n);
s.A_n = A_n;
s.B_n = B_n;

update_list(handles)
update_plot(handles)

s.iter = s.iter+1;
s.histRES(s.iter) = RES;
s.histACTION{s.iter} = 'R';
s.histOMEGA_n(1:length(s.OMEGA_n),s.iter) = s.OMEGA_n;
s.histXI_n(1:length(s.OMEGA_n),s.iter) = s.XI_n;
s.histA_n(1:length(s.OMEGA_n),s.iter) = s.A_n;
s.histB_n(1:length(s.OMEGA_n),s.iter) = s.B_n;
s.histOMEGA_n_0(1:length(s.OMEGA_n),s.iter) = s.OMEGA_n_0;
s.histdelta_freq(s.iter) = s.delta_freq;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CHANGING FITTING BANDWIDTH
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ChangeFittingBandwidth(handles,delta_freq)

global s

s.delta_freq = delta_freq;

if ~isempty(s.OMEGA_n)
    
    s.select = (sum(abs(ones(size(s.OMEGA_n_0))*s.freq-(s.OMEGA_n_0/(2*pi))*ones(size(s.freq))) < s.delta_freq/2,1))>0;
    
    [RES,A_n,B_n] = FindAmplitudes(s.OMEGA_n,s.XI_n);
    s.A_n = A_n;
    s.B_n = B_n;
    
    update_list(handles)
    update_plot(handles)
    
    s.iter = s.iter+1;
    s.histRES(s.iter) = RES;
    s.histACTION{s.iter} = 'B';
    s.histOMEGA_n(1:length(s.OMEGA_n),s.iter) = s.OMEGA_n;
    s.histXI_n(1:length(s.OMEGA_n),s.iter) = s.XI_n;
    s.histA_n(1:length(s.OMEGA_n),s.iter) = s.A_n;
    s.histB_n(1:length(s.OMEGA_n),s.iter) = s.B_n;
    s.histOMEGA_n_0(1:length(s.OMEGA_n),s.iter) = s.OMEGA_n_0;
    s.histdelta_freq(s.iter) = s.delta_freq;
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% REFINING THE IDENTIFICATION
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function RefineIdentification(handles)

global s

for n=1:length(s.OMEGA_n)
    
    set(handles.listboxModes,'Value',n)
    
    if ~s.locked(n)
        
        X0 = [1 1];
        
        Xsol = fminsearch(@(X) FindAmplitudes(...
            [ s.OMEGA_n([1:n-1]) ; X(1)*s.OMEGA_n(n) ; s.OMEGA_n([n+1:end])],...
            [ s.XI_n([1:n-1]) ; X(2)*s.XI_n(n) ; s.XI_n([n+1:end])])...
            ,X0,optimset('TolX',0.001,'Display','on'));
        omega_n = Xsol(1)*s.OMEGA_n(n);
        xi_n = Xsol(2)*s.XI_n(n);
        
        s.OMEGA_n(n) = omega_n;
        s.XI_n(n) = xi_n;
        
        
        [~,A_n,B_n] = FindAmplitudes(s.OMEGA_n,s.XI_n);
        s.A_n = A_n;
        s.B_n = B_n;
        
        update_list(handles)
        update_plot(handles)
        
    end
    
    
    
end


[RES,A_n,B_n] = FindAmplitudes(s.OMEGA_n,s.XI_n);
s.A_n = A_n;
s.B_n = B_n;

update_list(handles)
update_plot(handles)

s.iter = s.iter+1;
s.histRES(s.iter) = RES;
s.histACTION{s.iter} = 'F';
s.histOMEGA_n(1:length(s.OMEGA_n),s.iter) = s.OMEGA_n;
s.histXI_n(1:length(s.OMEGA_n),s.iter) = s.XI_n;
s.histA_n(1:length(s.OMEGA_n),s.iter) = s.A_n;
s.histB_n(1:length(s.OMEGA_n),s.iter) = s.B_n;
s.histOMEGA_n_0(1:length(s.OMEGA_n),s.iter) = s.OMEGA_n_0;
s.histdelta_freq(s.iter) = s.delta_freq;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% OTHER USEFUL FUNCTIONS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SwitchRealComplex(handles)

global s

s.modesC = ~s.modesC;

if length(s.OMEGA_n>0)
    [~,A_n,B_n] = FindAmplitudes(s.OMEGA_n,s.XI_n);
    s.A_n = A_n;
    s.B_n = B_n;
    
    update_list(handles)
    update_plot(handles)
end



function Z = Zmodal(omega,omega_n,xi_n,A_n,B_n)

global s

if ~s.modesC % modes reels
    
    Z = sum(...
        (A_n*(1i*omega))./...
        (omega_n.^2*ones(size(omega)) - ones(size(omega_n))*omega.^2 + 1i*2*(xi_n.*omega_n)*omega)...
        ,1);
    
else % modes complexes
    
    
    Z = sum(...
        (A_n*(1i*omega))./...
        (omega_n.^2*ones(size(omega)) - ones(size(omega_n))*omega.^2 + 1i*2*(xi_n.*omega_n)*omega)...
        +...
        (B_n*ones(size(omega)))./...
        (omega_n.^2*ones(size(omega)) - ones(size(omega_n))*omega.^2 + 1i*2*(xi_n.*omega_n)*omega)...
        ,1);
    
end



function [res,A_n,B_n] = FindAmplitudes(omega_n,xi_n)

global s

freq = s.freq(:).';
Ze = s.Ze(:).';

select = s.select;

omega = 2*pi*freq;

omega_n = omega_n(:);
xi_n = xi_n(:);


if ~s.modesC
    
    Ze_den = @(omega) (ones(size(omega_n))*(1i*omega))...
        ./...
        (omega_n.^2*ones(size(omega)) - ones(size(omega_n))*omega.^2 + 1i*2*(xi_n.*omega_n)*omega);
    
    A = [real(Ze_den(omega(select)).') ; imag(Ze_den(omega(select)).')];
    b = [real(Ze(select)).' ; imag(Ze(select)).'];
    
else
    
    Ze_den1 = @(omega) (ones(size(omega_n))*(1i*omega))...
        ./...
        (omega_n.^2*ones(size(omega)) - ones(size(omega_n))*omega.^2 + 1i*2*(xi_n.*omega_n)*omega);
    
    Ze_den2 = @(omega) (ones(size(omega_n))*ones(size(omega)))...
        ./...
        (omega_n.^2*ones(size(omega)) - ones(size(omega_n))*omega.^2 + 1i*2*(xi_n.*omega_n)*omega);
    
    A = [ [real(Ze_den1(omega(select)).') ; imag(Ze_den1(omega(select)).')] [real(Ze_den2(omega(select)).') ; imag(Ze_den2(omega(select)).')] ];
    b = [real(Ze(select)).' ; imag(Ze(select)).'];
    
end

x = A\b;
res = sum(abs(A*x-b).^2);
res = res/sum(select);

if ~s.modesC
    A_n = x;
    B_n = NaN*ones(size(x));
else
    N = length(x)/2;
    A_n = x(1:N);
    B_n = x(N+[1:N]);
end





% --------------------------------------------------------------------
function uipushtoolSave_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtoolSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global s

uisave('s');






function ModifyFRF(handles)

global s

s.Ze = s.Ze_original;
s.freq = s.freq_original;

if s.changesign
    s.Ze = -s.Ze;
    handles.menu_Tools_Change_the_sign.Checked = 'on';
else
    handles.menu_Tools_Change_the_sign.Checked = 'off';
end

if s.differentiate
    s.Ze = s.Ze.*(1j*2*pi*s.freq(:));
    handles.menu_Tools_Differentiate.Checked = 'on';
else
    handles.menu_Tools_Differentiate.Checked = 'off';
end

if s.integrate
    s.Ze = s.Ze./(1j*2*pi*s.freq(:));
    handles.menu_Tools_Integrate.Checked = 'on';
else
    handles.menu_Tools_Integrate.Checked = 'off';
end

if s.invert
    s.Ze = 1./s.Ze;
    handles.menu_Tools_Invert.Checked = 'on';
else
    handles.menu_Tools_Invert.Checked = 'off';
end


s.Ze = s.Ze(1:s.step:end);
s.freq = s.freq(1:s.step:end);

update_list(handles)
update_plot(handles)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Callbacks for menus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --------------------------------------------------------------------
function menu_Tools_Callback(hObject, eventdata, handles)
% hObject    handle to menu_Tools (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function menu_Tools_Differentiate_Callback(hObject, eventdata, handles)
% hObject    handle to menu_Tools_Differentiate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global s

s.differentiate = ~s.differentiate;
s.integrate = 0;

ModifyFRF(handles)


% --------------------------------------------------------------------
function menu_Tools_Integrate_Callback(hObject, eventdata, handles)
% hObject    handle to menu_Tools_Integrate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global s

s.integrate = ~s.integrate;
s.differentiate = 0;

ModifyFRF(handles)

% --------------------------------------------------------------------
function menu_Tools_Change_the_sign_Callback(hObject, eventdata, handles)
% hObject    handle to menu_Tools_Change_the_sign (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global s

s.changesign = ~s.changesign;

ModifyFRF(handles)







% --------------------------------------------------------------------
function menu_Tools_Modify_selected_pole_Callback(hObject, eventdata, handles)
% hObject    handle to menu_Tools_Modify_selected_pole (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global s

n = get(handles.listboxModes,'Value');

prompt = {'Frequency','Damping ratio'};
dlgtitle = 'Modify selected pole';
dims = [1 50; 1 50];
definput = {num2str(s.OMEGA_n(n)/(2*pi)),num2str(s.XI_n(n))};
answer = inputdlg(prompt,dlgtitle,dims,definput);

s.OMEGA_n(n) = 2*pi*str2num(answer{1});
s.XI_n(n) = str2num(answer{2});

[~,A_n,B_n] = FindAmplitudes(s.OMEGA_n,s.XI_n);
s.A_n = A_n;
s.B_n = B_n;

update_list(handles);
update_plot(handles);


% --------------------------------------------------------------------
function menu_Tools_Lock_Unlock_selected_pole_Callback(hObject, eventdata, handles)
% hObject    handle to menu_Tools_Lock_Unlock_selected_pole (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global s

n = get(handles.listboxModes,'Value');

if ~s.locked(n)
    s.locked(n) = 1;
    handles.menu_Tools_Change_the_sign.Checked = 'on';
else
    s.locked(n) = 0;
    handles.menu_Tools_Change_the_sign.Checked = 'off';
end

update_list(handles);
update_plot(handles);


% --------------------------------------------------------------------
function menu_Tools_Decimate_Callback(hObject, eventdata, handles)
% hObject    handle to menu_Tools_Decimate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global s

prompt = {'Enter a step to decimate the input data (1 = no decimation).'};
dlgtitle = 'Decimate';
dims = [1 50];
definput = {num2str(s.step)};
answer = inputdlg(prompt,dlgtitle,dims,definput);

s.step = str2num(answer{1});

ModifyFRF(handles)


% --------------------------------------------------------------------
function menu_Tools_Lock_all_poles_Callback(hObject, eventdata, handles)
% hObject    handle to menu_Tools_Lock_all_poles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global s

s.locked(:) = 1;

update_list(handles);
update_plot(handles);


% --------------------------------------------------------------------
function menu_Tools_Unlock_all_poles_Callback(hObject, eventdata, handles)
% hObject    handle to menu_Tools_Unlock_all_poles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global s

s.locked(:) = 0;

update_list(handles);
update_plot(handles);







% --------------------------------------------------------------------
function menu_View_Callback(hObject, eventdata, handles)
% hObject    handle to menu_View (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_View_set_axes_limits_Callback(hObject, eventdata, handles)
% hObject    handle to menu_View_set_axes_limits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global s

prompt = {'min Frequency','max Frequency','min Magnitude (leave empty for autoscale)','max Magnitude (leave empty for autoscale)'};
dlgtitle = 'Set axes limits';
dims = [1 50; 1 50; 1 50; 1 50];
definput = {num2str(s.fmin),num2str(s.fmax),num2str(s.magmin),num2str(s.magmax)};
answer = inputdlg(prompt,dlgtitle,dims,definput);

s.fmin = str2num(answer{1});
s.fmax = str2num(answer{2});

s.magmin = str2num(answer{3});
s.magmax = str2num(answer{4});

update_plot(handles);


% --------------------------------------------------------------------
function menu_Tools_Invert_Callback(hObject, eventdata, handles)
% hObject    handle to menu_Tools_Invert (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global s

s.invert = ~s.invert;

ModifyFRF(handles)


% --------------------------------------------------------------------
function menu_View_export_to_external_figure_Callback(hObject, eventdata, handles)
% hObject    handle to menu_View_export_to_external_figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

extfig = figure;

subplot(2,1,1)
pos1 = get(gca,'Position');
subplot(2,1,2)
pos2 = get(gca,'Position');

clf

copyobj(handles.axesTop,extfig);
copyobj(handles.axesBottom,extfig);

ax = get(gcf,'Children');

axes(ax(2)); % top
pos1(4) = 0.5;
pos1(2) = 0.45;
set(gca,'Position',pos1);

axes(ax(1)); % bottom
pos2(4) = 0.2;
set(gca,'Position',pos2);
xlabel('')
set(gca,'XTickLabel',{})

linkaxes(ax,'x');
set(gcf,'Units','pixels')
set(gcf,'Position',[0 0 700 500])

set(findall(get(gcf,'Children'),'-property','FontSize'),'Fontsize',14)



% --------------------------------------------------------------------
function meu_File_Callback(hObject, eventdata, handles)
% hObject    handle to meu_File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_Help_Callback(hObject, eventdata, handles)
% hObject    handle to menu_Help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_Help_About_Callback(hObject, eventdata, handles)
% hObject    handle to menu_Help_About (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

msg = cell(4,1);
msg{1} = sprintf('PeakPickingTools is a MATLAB toolbox to extract modal parameters from frequency response functions, using a robust peak-picking identification technique.\n');
msg{2} = sprintf('It is an implementation of the method described in:');
msg{3} = sprintf('Ablitzer F. 2021. Peak-picking identification technique for modal expansion of input impedance of brass instruments. Acta Acustica, 5, 53.\n');
msg{4} = sprintf('Copyright (c) 2022 Frédéric Ablitzer');
f = msgbox(msg, 'About','help');


% --------------------------------------------------------------------
function menu_File_Save_Callback(hObject, eventdata, handles)
% hObject    handle to menu_File_Save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uipushtoolSave_ClickedCallback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function menu_View_Hide_Fitting_Frequency_Domains_Callback(hObject, eventdata, handles)
% hObject    handle to menu_View_Hide_Fitting_Frequency_Domains (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if strcmp(handles.menu_View_Hide_Fitting_Frequency_Domains.Checked,'off')
    handles.menu_View_Hide_Fitting_Frequency_Domains.Checked = 'on';
else
    handles.menu_View_Hide_Fitting_Frequency_Domains.Checked = 'off';
end

update_plot(handles)



% --------------------------------------------------------------------
function menu_View_Hide_Modal_Contributions_Callback(hObject, eventdata, handles)
% hObject    handle to menu_View_Hide_Modal_Contributions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if strcmp(handles.menu_View_Hide_Modal_Contributions.Checked,'off')
    handles.menu_View_Hide_Modal_Contributions.Checked = 'on';
else
    handles.menu_View_Hide_Modal_Contributions.Checked = 'off';
end

update_plot(handles)


% --------------------------------------------------------------------
function menu_View_Hide_Reconstructed_FRF_Callback(hObject, eventdata, handles)
% hObject    handle to menu_View_Hide_Reconstructed_FRF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if strcmp(handles.menu_View_Hide_Reconstructed_FRF.Checked,'off')
    handles.menu_View_Hide_Reconstructed_FRF.Checked = 'on';
else
    handles.menu_View_Hide_Reconstructed_FRF.Checked = 'off';
end

update_plot(handles)


% --------------------------------------------------------------------
function menu_Dev_Callback(hObject, eventdata, handles)
% hObject    handle to menu_Dev (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_Dev_Execute_User_Function_1_Callback(hObject, eventdata, handles)
% hObject    handle to menu_Dev_Execute_User_Function_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

d = dir('dev');
list_functions = {};
for i=1:length(d)
    if ~d(i).isdir && strcmp(d(i).name(end-1:end),'.m')
        list_functions{length(list_functions)+1} = d(i).name;
    end
end

[indx,tf] = listdlg('PromptString',{'Select a user function to execute.'},...
    'SelectionMode','single',...
    'ListString',list_functions,...
    'Name','Execute user function',...
    'ListSize',[300,300]);

if tf
    feval(list_functions{indx}(1:end-2),handles)
end


% --------------------------------------------------------------------
function menu_Tools_Options_Callback(hObject, eventdata, handles)
% hObject    handle to menu_Tools_Options (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


global s

prompt = {'initial estimate for damping ratio'};
dlgtitle = 'Options';
dims = [1 50];
if isfield(s,'xi0')
    definput = {num2str(s.xi0)};
else
    definput = {num2str(0.05)};
end
answer = inputdlg(prompt,dlgtitle,dims,definput);

s.xi0 = str2num(answer{1});




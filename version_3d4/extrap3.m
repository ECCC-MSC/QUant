%
% Version Beta 3.4
% David S. Mueller 5/15/2012
% 

function varargout = extrap3(varargin)
% EXTRAP3 M-file for extrap3.fig
%      EXTRAP3, by itself, creates a new EXTRAP3 or raises the existing
%      singleton*.
%
%      H = EXTRAP3 returns the handle to a new EXTRAP3 or the handle to
%      the existing singleton*.
%
%      EXTRAP3('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EXTRAP3.M with the given input arguments.
%
%      EXTRAP3('Property','Value',...) creates a new EXTRAP3 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before extrap3_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to extrap3_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools File.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help extrap3

% Last Modified by GUIDE v2.5 30-Apr-2012 17:59:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @extrap3_OpeningFcn, ...
                   'gui_OutputFcn',  @extrap3_OutputFcn, ...
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

%==========================================================================
% --- Executes just before extrap3 is made visible.
function extrap3_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
%==========================================================================

% Choose default command line output for extrap3
handles.output = hObject;
%
% Initialize custom handles data structure variables
% --------------------------------------------------
handles.version='extrap - Beta 3.4';
set(handles.figure1,'Name',handles.version);
handles.draftUnits='M';
handles.threshold=20; %Threshold set to 20% of median # of pts.
handles.extents=[0 100];
handles.pathName=pwd;
handles.backup=handles;
%
% Update handles structure
% ------------------------
guidata(hObject, handles);



% *******************
% File Menu Callbacks
% *******************

%==========================================================================
function mmt_Callback(hObject, eventdata, handles)
%
% Loads and processes transect data from TRDI *.mmt file. The draft,
% pathname, and transect filenames are all that is used from the mmt file.
%==========================================================================

   set(handles.figure1,'Pointer','watch');
   drawnow;
% Clear gui and handles data structure
% ------------------------------------
handles=resetgui(handles);

% Prompt user for mmt file to process
% -----------------------------------
[filename,pathname]=uigetfile('*.mmt',  'Multiselect','on','SELECT FILE FOR EVALUATION',handles.pathName);

if ~filename==0
    set(handles.figure1,'Name',[handles.version '  ---  ' pathname]);

    % Load mmt file
    % -------------
    infile=strcat(pathname,filename);
    [~,~,MMT_Transects,~,MMT_Active_Config,~,MMT_Summary_BT,~,~,~,~,~,~]=mmt2mat(infile);  

    % Extract need information from mmt file
    % --------------------------------------
    idxChecked=find(MMT_Transects.Checked==1);
    fileName=MMT_Transects.Files(idxChecked);
    pathName=pathname;
    draft=MMT_Active_Config.Offsets_Transducer_Depth(idxChecked);
    beginDist=MMT_Active_Config.Edge_Begin_Shore_Distance(idxChecked);
    endDist=MMT_Active_Config.Edge_End_Shore_Distance(idxChecked);
    beginLeft=MMT_Summary_BT.Begin_Left(idxChecked);
    idxLeft= beginLeft==1;
    idxRight= beginLeft==0;
    startBank=cell(size(beginLeft));
    startBank(idxRight)={'Right'};
    startBank(idxLeft)={'Left'};
    nfiles=length(idxChecked);
    fileName{nfiles+1}='Measurement';

    % Populate the transect list box in the GUI
    % -----------------------------------------
    set(handles.listFiles,'String',fileName);
    drawnow;

    % Store data in handles data structure
    % ------------------------------------
    handles.fileName=fileName;
    handles.pathName=pathName;
    handles.draftUnits='M';
    handles.draft=draft;
    handles.fileType=1;
    handles.nfiles=nfiles;

    % Initialize transect data object
    % -------------------------------
    transData(1:nfiles)=OriginData();

    % Begin reading transect data
    % ----------------------------
    set(handles.puFit,'Value',1);
    set(handles.puTop,'Enable','off');
    set(handles.puBottom,'Enable','off');
    set(handles.edExp,'Enable','off');
    set(handles.cbOpt,'Enable','off'); 
    for ifile=1:nfiles

        % Read data file
        % --------------
        set(handles.listFiles,'Value',ifile);
        drawnow;
        transData(ifile)=OriginData(fileName{ifile},pathName,draft(ifile),beginDist(ifile),endDist(ifile),startBank(ifile));

        % Set Fit to Automatic
        % --------------------
        handles.selFit(ifile)=SelectFit();
        handles.selFit(ifile).fitMethod='Automatic';
    end
    handles.selFit(nfiles+1)=SelectFit();
    handles.selFit(nfiles+1).fitMethod='Automatic';

    % Store data in handles data structure
    % ------------------------------------
    handles.transData=transData;

    % Update handles structure
    % ------------------------
    guidata(hObject, handles);

    % Process transects
    % -----------------
    [normData,selFit,sensData]=processProfiles(handles);

    % Store processed data in handles data structure
    % ----------------------------------------------
    handles.normData=normData;
    handles.selFit=selFit;
    handles.sensData=sensData;

    % Update handles structure
    % ------------------------
    guidata(hObject, handles);

    % Plot data in graphs
    % -------------------
    plotTrans(handles)
    %
    % Update handles structure
    % ------------------------
    guidata(hObject, handles);
end
   set(handles.figure1,'Pointer','arrow');
   drawnow;

%==========================================================================
function raw_Callback(hObject, eventdata, handles)
%
% Loads and processes raw data files from TRDI (*.pd0 and *r.000). This
% function requires the user to enter the draft for each transect.
%==========================================================================
%
% Clear gui and handles data structure
% ------------------------------------

set(handles.figure1,'Pointer','watch');
drawnow;
% Clear gui and handles data structure
% ------------------------------------
handles=resetgui(handles);

% Open file select dialog and allow multiple files to be selected
% ---------------------------------------------------------------
[fileName,pathName]=uigetfile('*r.000; *.pd0',  'Multiselect','on','SELECT FILE(S) FOR EVALUATION',handles.pathName);
if ~fileName==0
    set(handles.figure1,'Name',[handles.version '  ---  ' pathName]);

    % Determine Number of Files to be Processed.
    % ------------------------------------------
    if  isa(fileName,'cell')
        nfiles=size(fileName,2); 
        fileName=sort(fileName);       
    else
        nfiles=1;
        fileName={fileName};
    end
    draft=zeros(nfiles,1);
    fileName{nfiles+1}='Measurement';

    % Populate the transect list box in the GUI
    % -----------------------------------------
    set(handles.listFiles,'String',fileName);
    drawnow;

    for ifile=1:nfiles

        % Request user to enter draft
        % ---------------------------
        if ifile==1
           units=questdlg('Select units for draft input','Draft Units', 'Meters','Feet',nan);
           draftUnits=units(1);
           draft(ifile)=0;
            if strcmp(draftUnits,'F') 
                draft(ifile)=str2double(inputdlg(['Enter draft in feet for ',char(fileName(ifile))],'Draft Input'));
                draft(ifile)=draft(ifile)./3.281;
            else
                draft(ifile)=str2double(inputdlg(['Enter draft in meters for ',char(fileName(ifile))],'Draft Input'));     
            end       
        elseif ifile>1
            draft(ifile)=draft(ifile-1);
            if strcmp(draftUnits,'F') 
                draft(ifile)=str2double(inputdlg(['Enter draft in feet for ',char(fileName(ifile))],'Draft Input',1,{num2str(draft(ifile-1).*3.281,4)}));
                draft(ifile)=draft(ifile)./3.281;
            else
                draft(ifile)=str2double(inputdlg(['Enter draft in meters for ',char(fileName(ifile))],'Draft Input',1,{num2str(draft(ifile-1),4)}));     
            end
        end
    end

    % Save data to handles structure
    % ------------------------------
    handles.fileName=fileName;
    handles.pathName=pathName;
    handles.draftUnits=draftUnits;
    handles.draft=draft;
    handles.fileType=1;
    handles.nfiles=nfiles;

    % Load data
    % ---------
    transData(1:nfiles)=OriginData();
    for ifile=1:nfiles

        % Read data file
        % --------------
        set(handles.listFiles,'Value',ifile);
        drawnow;    
        transData(ifile)=OriginData(fileName{ifile},pathName,draft(ifile),0,0,'Unknown');

        % Set Fit to Automatic
        % --------------------
        handles.selFit(ifile)=SelectFit();
        handles.selFit(ifile).fitMethod='Automatic';
    end
    handles.selFit(nfiles+1)=SelectFit();
    handles.selFit(nfiles+1).fitMethod='Automatic';

    % Store transData in handles structure
    % ------------------------------------
    handles.transData=transData;

    % Update handles structure
    % ------------------------
    guidata(hObject, handles);

    % Process transects
    % -----------------
    [normData,selFit,sensData]=processProfiles(handles);
    %
    % Store results in handles data structure
    % ---------------------------------------
    handles.normData=normData;
    handles.selFit=selFit;
    handles.sensData=sensData;

    % Update handles structure
    % ------------------------
    guidata(hObject, handles);

    % Plot data
    % ---------
    plotTrans(handles)


    % Update handles structure
    % ------------------------
    guidata(hObject, handles);
end

    % Update cursor
    % -------------
    set(handles.figure1,'Pointer','arrow');
    drawnow;

%==========================================================================
function mat_Callback(hObject, eventdata, handles)
% 
% Loads and processes Matlab output from RiverSurveyor Live
%==========================================================================
%
% Clear gui and handles data structure
% ------------------------------------

% Set cursor to busy
% ------------------
set(handles.figure1,'Pointer','watch');
drawnow;

% Clear gui and reset handles
% ---------------------------
handles=resetgui(handles);

% Get names and path of transects to be processed
% -----------------------------------------------
[fileName,pathName]=uigetfile('*.mat',  'Multiselect','on','SELECT FILE(S) FOR EVALUATION',handles.pathName);
if ~fileName==0
    set(handles.figure1,'Name',[handles.version '  ---  ' pathName]);

    % Determine Number of Files to be Processed.
    % ------------------------------------------
    if  isa(fileName,'cell')
        nfiles=size(fileName,2); 
        fileName=sort(fileName);       
    else
        nfiles=1;
        fileName={fileName};
    end
    fileName{nfiles+1}='Measurement';
    %
    % Update display
    % --------------
    set(handles.listFiles,'String',fileName);
    drawnow;
    %
    % Store data in handles data structure
    % ------------------------------------
    handles.fileName=fileName;
    handles.pathName=pathName;
    handles.draftUnits=nan;
    handles.draft=nan(nfiles,1);
    handles.fileType=3;
    handles.nfiles=nfiles;
    %
    % Load data
    % ---------
    transData(1:nfiles)=OriginData();
    for ifile=1:nfiles

        % Read data file
        % --------------
        set(handles.listFiles,'Value',ifile);
        drawnow;    
        transData(ifile)=OriginData(fileName{ifile},pathName);

        % Set Fit to Automatic
        % --------------------
        handles.selFit(ifile)=SelectFit();
        handles.selFit(ifile).fitMethod='Automatic';
    end
    handles.selFit(nfiles+1)=SelectFit();
    handles.selFit(nfiles+1).fitMethod='Automatic';

    % Store transData in handles data structure;
    % ------------------------------------------
    handles.transData=transData;

    % Update handles structure
    % ------------------------
    guidata(hObject, handles);
    %
    % Process transects
    % -----------------
    [normData,selFit,sensData]=processProfiles(handles);
    %
    % Store results in handles data structure
    % ---------------------------------------
    handles.normData=normData;
    handles.selFit=selFit;
    handles.sensData=sensData;

    % Update handles structure
    % ------------------------
    guidata(hObject, handles);

    % Plot data
    % ---------
     plotTrans(handles)

    % Update handles structure
    % ------------------------
    guidata(hObject, handles);
end
    set(handles.figure1,'Pointer','arrow');
    drawnow;
%==========================================================================
function SaveSum_Callback(hObject, eventdata, handles)
% 
% Saves a summary of the data to an ASCII file and appends data to an
% ongoing log stored as a Matlab file.
%==========================================================================
%
% Update display
% --------------

set(handles.figure1,'Pointer','watch');
drawnow;

%
% Prompt for user opinion and comments
% ------------------------------------
agree=questdlg('Do you agree with the automatic fit?','User Opinion','Yes','No','Yes');
comment=inputdlg('Enter comments:','User comments',10);
drawnow
%
% Prepare fit statistics output table
% -----------------------------------
nfiles=handles.nfiles;
for ifiles=1:nfiles+1
    fitOutTable{ifiles+1,1}=handles.selFit(ifiles).fileName;
    fitOutTable{ifiles+1,2}=handles.selFit(ifiles).topMethod;
    fitOutTable{ifiles+1,3}=handles.selFit(ifiles).botMethod;
    fitOutTable{ifiles+1,4}=handles.selFit(ifiles).exponent;
    fitOutTable{ifiles+1,5}=handles.selFit(ifiles).fitMethod;
    fitOutTable{ifiles+1,6}=handles.selFit(ifiles).fitrsqr;
    fitOutTable{ifiles+1,7}=handles.selFit(ifiles).topmaxdiff;
    fitOutTable{ifiles+1,8}=handles.selFit(ifiles).topr2;
    fitOutTable{ifiles+1,9}=handles.selFit(ifiles).botdiff;
    fitOutTable{ifiles+1,10}=handles.selFit(ifiles).botrsqr;
    fitOutTable{ifiles+1,11}=handles.selFit(ifiles).ppexponent;
    fitOutTable{ifiles+1,12}=handles.selFit(ifiles).nsexponent;
    fitOutTable{ifiles+1,13}=handles.threshold;
    fitOutTable{ifiles+1,14}=handles.normData(ifiles).dataExtent(1);
    fitOutTable{ifiles+1,15}=handles.normData(ifiles).dataExtent(2);
    fitOutTable{ifiles+1,16}=handles.selFit(ifiles).topr2;
    fitOutTable{ifiles+1,17}=handles.selFit(ifiles).topfitr2;
end
fitOutTable{1,1}='Filename';
fitOutTable{1,2}='Top';
fitOutTable{1,3}='Bottom';
fitOutTable{1,4}='Exponent';
fitOutTable{1,5}='Fit Type';
fitOutTable{1,6}='Fit r^2';
fitOutTable{1,7}='Top Diff';
fitOutTable{1,8}='Top Fit';
fitOutTable{1,9}='Bot Diff';
fitOutTable{1,10}='Bot r^2';
fitOutTable{1,11}='PP Opt';
fitOutTable{1,12}='NS Opt';
fitOutTable{1,13}='Threshold';
fitOutTable{1,14}='Lower Ext.';
fitOutTable{1,15}='Upper Ext.';
fitOutTable{1,16}='Top R^2';
fitOutTable{1,17}='Top Cust r';

%
% Create name and open ASCII output file
% --------------------------------------
outName=[handles.pathName datestr(now,'yyyymmddHHMMSS_extrap') '.txt'];
fid=fopen(outName, 'w');
%
% Write data to ASCII output file
% -------------------------------
fprintf(fid,datestr(now));
fprintf(fid,'\r\n');
fprintf(fid,handles.version);
fprintf(fid,'\r\n');
fprintf(fid,'\r\n');
% fprintf(fid,'************************************** \r\n');
% fprintf(fid,'Auto Recommendation for Extrapolations \r\n');
fprintf(fid,'************************************** \r\n');
fprintf(fid,'   Top Extrapolation: %s \r\n',fitOutTable{end,2});
fprintf(fid,'Bottom Extraploation: %s \r\n',fitOutTable{end,3});
fprintf(fid,'            Exponent: %6.4f \r\n',fitOutTable{end,4});
fprintf(fid,'            Fit Type: %s \r\n',fitOutTable{end,5});
fprintf(fid,'************************************** \r\n');
fprintf(fid,'\r\n');
fprintf(fid,'\r\n');
fprintf(fid,'Fit Analysis Statistics \r\n');
fprintf(fid,'----------------------- \r\n');
fprintf(fid,'%-30s %-10s %-10s %-10s %-10s  %-10s \r\n',fitOutTable{1,1:6});
for ifile=2:nfiles+2
    fprintf(fid,'%-30s %-10s %-10s %-10.4f %-10s %-10.4f \r\n',fitOutTable{ifile,1:6});
end
fprintf(fid,'\r\n');
fprintf(fid,'\r\n');
fprintf(fid,'%-30s %-10s %-10s %-10s %-10s \r\n',fitOutTable{1,1},fitOutTable{1,7:10});
for ifile=2:nfiles+2
    fprintf(fid,'%-30s %-10.4f %-10.4f %-10.4f %-10.4f \r\n',fitOutTable{ifile,1},fitOutTable{ifile,7:10});
end
fprintf(fid,'\r\n');
fprintf(fid,'\r\n');
fprintf(fid,'%-30s %-10s %-10s %-10s \r\n',fitOutTable{1,1},fitOutTable{1,13:15});
for ifile=2:nfiles+2
    fprintf(fid,'%-30s %-10.4f %-10.4f %-10.4f \r\n',fitOutTable{ifile,1},fitOutTable{ifile,13:15});
end
fprintf(fid,'\r\n');
fprintf(fid,'\r\n');
fprintf(fid,'Discharge Sensitivity Analysis \r\n');
fprintf(fid,'------------------------------ \r\n');
fprintf(fid,'Top        Bottom     Exponent   %% Difference \r\n');
 qtable=get(handles.q_table,'Data');
 rows=size(qtable,1);
for ifile=1:rows
fprintf(fid,'%-10s %-10s %-10s %12s \r\n',qtable{ifile,1:4});
end
fclose(fid);
%
% Prepare data for storing in Matlab file
% ---------------------------------------
for ifile=1:nfiles
    btVel=handles.transData(ifile).btVel;
    ensDeltaTime=handles.transData(ifile).ensDeltaTime;
    depthEns=handles.transData(ifile).depthEns;
    beamDepths=handles.transData(ifile).beamDepths;
    %
    % Compute width of cross section
    % ------------------------------
    numEns=size(btVel,2);
    ensDeltaTimeadj=ensDeltaTime;
    ensDeltaTimeadj(isnan(btVel(1,:)))=nan;
    k=1;
    for j=2:numEns
        if j+1<numEns && isnan(ensDeltaTimeadj(j))
            k=k+1;
        else
            ensDeltaTimeadj(j)=ensDeltaTimeadj(j).*k;
            k=1;
        end
    end
    xDist=btVel(1,:).*ensDeltaTimeadj;
    xDist=nancumsum(xDist);
    yDist=btVel(2,:).*ensDeltaTimeadj;
    yDist=nancumsum(yDist);
    dmg=sqrt(xDist.^2+yDist.^2);
    width(ifile)=handles.transData(ifile).beginDist+dmg(end)+handles.transData(ifile).endDist;
    %
    % Compute other channel characteristics
    % -------------------------------------
    meanDepth(ifile)=nanmean(depthEns); 
    maxDepth(ifile)=nanmax(depthEns);
    rough2(ifile)=nanmean(nanstd(beamDepths));
end
%
% Compute mean channel characteristics
% ------------------------------------
width(nfiles+1)=nanmean(width);
meanDepth(nfiles+1)=nanmean(meanDepth);
maxDepth(nfiles+1)=nanmean(maxDepth);
rough2(nfiles+1)=nanmean(rough2);
%
% Prepare stream characteristics table
% ------------------------------------
temp=fitOutTable(2:end,:);
temp2=fitOutTable{end,1};
temp3=fitOutTable{end-1,1};
temp4=[temp3(1:end-4) '_' temp2];
temp{end,1}=temp4;
qtemp{1,1}=temp4;
qtemp=[qtemp qtable(:,end)'];
for ifile=1:nfiles+1
    StreamCharNew(ifile,1)={temp(ifile,1)};
    StreamCharNew(ifile,2)={width(ifile)};
    StreamCharNew(ifile,3)={meanDepth(ifile)};
    StreamCharNew(ifile,4)={maxDepth(ifile)};
    StreamCharNew(ifile,5)={rough2(ifile)};
end
%


% Store data in new or existing Matlab file
% ------------------------------------------

    fitOutTable{end,1}=temp4;
    StreamChar=[{'Filename' 'Width' 'Mean Depth' 'Max Depth' 'Rough'} ;StreamCharNew];
    Agreement=[{'Filename' 'User Opinion' 'Comments'};[{temp4 agree} comment]];
    FitStats=fitOutTable;
    QStats=[{'Filename' 'PP' 'PPOpt' 'CNS' 'CNSOpt' '3PNS' '3PNSOpt' };[temp4 qtable(:,end)']];

    save ('extrap_summary.mat','FitStats','QStats','Agreement', 'StreamChar');

% 
% Update display
% --------------
   set(handles.figure1,'Pointer','arrow');
   drawnow;

%==========================================================================
function FileExit_Callback(hObject, eventdata, handles)
% 
% Close all figures and data then close the program
%==========================================================================
fclose all;
close all;

% ************************
% Configure Menu Callbacks
% ************************

%==========================================================================
function Threshold_Callback(hObject, eventdata, handles)
% 
% Displays input dialog for user to change the number of points that must
% be used in computing a median value, in order for that median to be used
% in the fitting process. The graphs and fits are updated.
%==========================================================================

   set(handles.figure1,'Pointer','watch');
   drawnow;
% Prompt for new threshold value, showing existing value as default
% -----------------------------------------------------------------
threshold=handles.threshold;
temp1=sprintf('Current threshold is %i \nEnter new threshold value:',threshold);
temp=inputdlg(temp1,'Threshold',1,{num2str(threshold)});
 
% Store new threshold in handles data structure
% ---------------------------------------------
handles.threshold=str2double(temp);

% Update handles structure
% ------------------------
guidata(hObject, handles);

% Reprocess data with new threshold value
% ---------------------------------------
[handles.normData,handles.selFit,handles.sensData]=processProfiles(handles);

% Update handles structure
% ------------------------
guidata(hObject, handles);

% Update graphs
% -------------
plotTrans(handles)
   set(handles.figure1,'Pointer','arrow');
   drawnow;

% Update handles structure
% ------------------------
guidata(hObject, handles);

%==========================================================================
function Discharge_Callback(hObject, eventdata, handles)
%
% Changes data processing to use discharge if not alread selected
%==========================================================================

   set(handles.figure1,'Pointer','watch');
   drawnow;
% Set menu display so that discharge is checked
% ---------------------------------------------
set(handles.Discharge,'Checked','on');
set(handles.Velocity,'Checked','off');

% Update handles structure
% ------------------------
guidata(hObject, handles);

% Reprocess data
% --------------
[normData,selFit,sensData]=processProfiles(handles);
%
% Store results in handles data structure
% ---------------------------------------
handles.normData=normData;
handles.selFit=selFit;
handles.sensData=sensData;
set(handles.figure1,'Pointer','arrow');
drawnow
%
% Update handles structure
% ------------------------
guidata(hObject, handles);
%
% Update display
% --------------
plotTrans(handles)
   

drawnow

%==========================================================================
function Velocity_Callback(hObject, eventdata, handles)
%
% Changes data processing to use velocity if not alread selected
%==========================================================================
%
% Set menu display so that velocity is checked
% --------------------------------------------

   set(handles.figure1,'Pointer','watch');
   drawnow;
set(handles.Discharge,'Checked','off');
set(handles.Velocity,'Checked','on');
%
% Update handles structure
% ------------------------
guidata(hObject, handles);
%
% Update display
% --------------
% set(handles.txtmsg,'String','PROCESSING');
%
% Reprocess data
% --------------
[normData,selFit,sensData]=processProfiles(handles);
%
% Store results in handles data structure
% ---------------------------------------
handles.normData=normData;
handles.selFit=selFit;
handles.sensData=sensData;
set(handles.figure1,'Pointer','arrow');
drawnow;
%
% Update handles structure
% ------------------------
guidata(hObject, handles);
%
% Update display
% --------------
plotTrans(handles)
   
drawnow

%==========================================================================
function subsection_Callback(hObject, eventdata, handles)
% 
% Allows the user to specify the upper and lower limit of cumulative
% discharge for which the user would like to analyze the data. I.E. 25 and
% 75 would allow the user to look at the center 50% for the flow.
%==========================================================================

% Set cursor to busy
% ------------------
set(handles.figure1,'Pointer','watch');
drawnow;
% Show existing extents in the user dialog for new extents
% ---------------------------------------------------------
extents=handles.extents;
temp=inputdlg({'Enter lower discharge limit in percent:';'Enter upper discharge limit in percent:'}...
    ,'Discharge Extents',1,{num2str(extents(1)) num2str(extents(2))});
 
% Store new extents in handles data structure
% ---------------------------------------------
handles.extents(1)=str2double(temp(1));
handles.extents(2)=str2double(temp(2));

% Update handles structure
% ------------------------
guidata(hObject, handles);

% Reprocess data with new threshold value
% ---------------------------------------
[handles.normData,handles.selFit,handles.sensData]=processProfiles(handles);

% Update handles structure
% ------------------------
guidata(hObject, handles);

% Update graphs
% -------------
plotTrans(handles)
   set(handles.figure1,'Pointer','arrow');
   drawnow;

% Update handles structure
% ------------------------
guidata(hObject, handles);

% *******************
% View Menu Callbacks
% *******************

%==========================================================================
function CellData_Callback(hObject, eventdata, handles)
% 
% Toggles display of the cell data on and off in the transect plot
%==========================================================================

   set(handles.figure1,'Pointer','watch');
   drawnow;
% Change displayed status of CellData menu item
% ---------------------------------------------
if strcmp(get(handles.CellData,'Checked'),'on');
    set(handles.CellData,'Checked','off');
else
    set(handles.CellData,'Checked','on');
end

% Update handles structure
% ------------------------
guidata(hObject, handles);

% Update display
% --------------
plotTrans(handles)
   set(handles.figure1,'Pointer','arrow');
drawnow

%==========================================================================
function SurfCells_Callback(hObject, eventdata, handles)
% 
% Toggles on and off the highlighting of the surface cells in the transect
% plot
%==========================================================================

   set(handles.figure1,'Pointer','watch');
   drawnow;
% Change displayed status of SurfCells menu item
% ----------------------------------------------
if strcmp(get(handles.SurfCells,'Checked'),'on');
    set(handles.SurfCells,'Checked','off');
else
    set(handles.SurfCells,'Checked','on');
end

% Update handles structure
% ------------------------
guidata(hObject, handles);

% Update display
% --------------
plotTrans(handles);
   set(handles.figure1,'Pointer','arrow');
drawnow

%==========================================================================
% --------------------------------------------------------------------
function TransFit_Callback(hObject, eventdata, handles)
% 
% Toggles the automatically selected fit to display in the transect plot
%==========================================================================

   set(handles.figure1,'Pointer','watch');
   drawnow;
% Change displayed status of TransFit menu item
% ---------------------------------------------
if strcmp(get(handles.TransFit,'Checked'),'on');
    set(handles.TransFit,'Checked','off');
else
    set(handles.TransFit,'Checked','on');
end

% Update handles structure
% ------------------------
guidata(hObject, handles);

% Update display
% --------------
plotTrans(handles)
   set(handles.figure1,'Pointer','arrow');
drawnow

%==========================================================================
function MeasFit_Callback(hObject, eventdata, handles)
% 
% Toggles display of the composite automatically selected fit in the
% transect plot
%==========================================================================

   set(handles.figure1,'Pointer','watch');
   drawnow;
% Update menu items to reflect selection
% --------------------------------------
if strcmp(get(handles.MeasFit,'Checked'),'on');
    set(handles.MeasFit,'Checked','off');
else
    set(handles.MeasFit,'Checked','on');
end

% Update handles structure
% ------------------------
guidata(hObject, handles);

% Update display
% --------------
plotTrans(handles)
   set(handles.figure1,'Pointer','arrow');
drawnow

%==========================================================================
function TransMed_Callback(hObject, eventdata, handles)
% 
% Turns on and off display of the median points for each transect in the
% composite measurement plot
%==========================================================================

   set(handles.figure1,'Pointer','watch');
   drawnow;
% Update menu items to reflect selection
% --------------------------------------
if strcmp(get(handles.TransMed,'Checked'),'on');
    set(handles.TransMed,'Checked','off');
else
    set(handles.TransMed,'Checked','on');
end

% Update handles structure
% ------------------------
guidata(hObject, handles);

% Update display
% --------------
plotTrans(handles)
   set(handles.figure1,'Pointer','arrow');
   
drawnow

%==========================================================================
function MeasMed_Callback(hObject, eventdata, handles)
% 
% Turns on and off display of the median points for composite measurement 
% data in the composite measurement plot
%==========================================================================


% Set cursor to busy
% ------------------
set(handles.figure1,'Pointer','watch');
drawnow;

% Update menu items to reflect selection
% --------------------------------------
if strcmp(get(handles.MeasMed,'Checked'),'on');
    set(handles.MeasMed,'Checked','off');
else
    set(handles.MeasMed,'Checked','on');
end

% Update handles structure
% ------------------------
guidata(hObject, handles);

% Update display
% --------------
plotTrans(handles)

% Return cursor to normal
% -----------------------
set(handles.figure1,'Pointer','arrow');
drawnow;


% *************
% GUI Callbacks
% *************
%==========================================================================
function cbOpt_Callback(hObject, eventdata, handles)
%==========================================================================

    % Retrieve data
    % -------------
    ifile=get(handles.listFiles,'Value');
    normData=handles.normData;
    selFit=handles.selFit;
    
    % Set new exponent and recompute
    % ------------------------------
    if get(handles.cbOpt,'Value')==1
        selFit(ifile).expMethod='Optimize';
    else
        selFit(ifile).expMethod='Manual';    
        expdlg=questdlg('Set exponent to 1/6th?','Exponent Setting','Yes','No','Yes');
        if strcmp(expdlg,'Yes')
            selFit(ifile).exponent=0.1667;
        end
    end
    
    set(handles.figure1,'Pointer','watch');
    drawnow;
    selFit(ifile)=SelectFit(normData(ifile),selFit(ifile).fitMethod,selFit(ifile));
    set(handles.edExp,'String',num2str(selFit(ifile).exponent,'%6.4f'));
    handles.selFit=selFit;
    
     % Compute and display discharge sensitivity
    % -----------------------------------------
    if strcmpi(handles.normData(ifile).dataType,'q')
        sensData=QSensitivity(handles.transData,selFit);
    else
        sensData=QSensitivity();
    end
    h=handles.q_table;
    DisplayTable(sensData,h);
    handles.sensData=sensData;

    % Update handles structure
    % ------------------------
    guidata(hObject, handles);

    % Update graphs
    % -------------
     plotTrans(handles)

    % Update handles structure
    % ------------------------
    guidata(hObject, handles);

    set(handles.figure1,'Pointer','arrow');
    drawnow;
    
%==========================================================================
function edExp_Callback(hObject, eventdata, handles)
%==========================================================================
    
    set(handles.figure1,'Pointer','watch');
    drawnow;
    % Retrieve data
    % -------------
    ifile=get(handles.listFiles,'Value');
    normData=handles.normData;
    selFit=handles.selFit;
    
    % Set new exponent and recompute
    % ------------------------------
    selFit(ifile).exponent=str2double(get(handles.edExp,'String'));
    selFit(ifile)=SelectFit(normData(ifile),selFit(ifile).fitMethod,selFit(ifile));
    handles.selFit=selFit;
    
     % Compute and display discharge sensitivity
    % -----------------------------------------
    if strcmpi(handles.normData(ifile).dataType,'q')
        sensData=QSensitivity(handles.transData,selFit);
    else
        sensData=QSensitivity();
    end
    h=handles.q_table;
    DisplayTable(sensData,h);
    handles.sensData=sensData;

    % Update handles structure
    % ------------------------
    guidata(hObject, handles);

    % Update graphs
    % -------------
     plotTrans(handles)

    % Update handles structure
    % ------------------------
    guidata(hObject, handles);

    set(handles.figure1,'Pointer','arrow');
    drawnow;
    
%==========================================================================
function puBottom_Callback(hObject, eventdata, handles)
%==========================================================================
    set(handles.figure1,'Pointer','watch');
    drawnow;
    % Retrieve data
    % -------------
    ifile=get(handles.listFiles,'Value');
    normData=handles.normData;
    selFit=handles.selFit;
    
    % Set new exponent and recompute
    % ------------------------------
    boti=get(handles.puBottom,'Value');
    botMethod=get(handles.puBottom,'String');
    selFit(ifile).botMethod=botMethod{boti};
    if get(handles.cbOpt,'Value')==1
        selFit(ifile).expMethod='Optimize';
    else
        selFit(ifile).exponent=0.1667;
    end
    if boti==2
        selFit(ifile).topMethod='Constant';
        set(handles.puTop,'Value',2);
    end
    
    selFit(ifile)=SelectFit(normData(ifile),selFit(ifile).fitMethod,selFit(ifile)); 
    set(handles.edExp,'String',num2str(selFit(ifile).exponent,'%6.4f'));
    handles.selFit=selFit;
    
    % Compute and display discharge sensitivity
    % -----------------------------------------
    if strcmpi(handles.normData(ifile).dataType,'q')
        sensData=QSensitivity(handles.transData,selFit);
    else
        sensData=QSensitivity();
    end
    h=handles.q_table;
    DisplayTable(sensData,h);
    handles.sensData=sensData;

    % Update handles structure
    % ------------------------
    guidata(hObject, handles);

    % Update graphs
    % -------------
     plotTrans(handles)

    % Update handles structure
    % ------------------------
    guidata(hObject, handles);
    
    set(handles.figure1,'Pointer','arrow');
    drawnow;
    
%==========================================================================
function puTop_Callback(hObject, eventdata, handles)
%==========================================================================
    
    set(handles.figure1,'Pointer','watch');
    drawnow;
    % Retrieve data
    % -------------
    ifile=get(handles.listFiles,'Value');
    normData=handles.normData;
    selFit=handles.selFit;
    
    % Set new exponent and recompute
    % ------------------------------
    topi=get(handles.puTop,'Value');
    topMethods=get(handles.puTop,'String');
    selFit(ifile).topMethod=topMethods{topi};
    if topi==1
        selFit(ifile).botMethod='Power';
        set(handles.puBottom,'Value',1);
    end
    selFit(ifile)=SelectFit(normData(ifile),selFit(ifile).fitMethod,selFit(ifile));
    handles.selFit=selFit;

    % Compute and display discharge sensitivity
    % -----------------------------------------
    if strcmpi(handles.normData(ifile).dataType,'q')
        sensData=QSensitivity(handles.transData,selFit);
    else
        sensData=QSensitivity();
    end
    h=handles.q_table;
    DisplayTable(sensData,h);
    handles.sensData=sensData;
    
    % Update handles structure
    % ------------------------
    guidata(hObject, handles);

    % Update graphs
    % -------------
     plotTrans(handles)

    % Update handles structure
    % ------------------------
    guidata(hObject, handles);
    
    set(handles.figure1,'Pointer','arrow');
    drawnow;
%==========================================================================
function puFit_Callback(hObject, eventdata, handles)
%==========================================================================
% ID data to be used
    % ------------------
   set(handles.figure1,'Pointer','watch');
   drawnow;
    j=get(handles.listFiles,'Value');
    selFit=handles.selFit;

    % Fit set to Automatic
    % --------------------
    if get(handles.puFit,'Value')==1
        set(handles.puTop,'Enable','off');
        set(handles.puBottom,'Enable','off');
        set(handles.edExp,'Enable','off');
        set(handles.cbOpt,'Enable','off'); 
        set(handles.cbOpt,'Value',0);
        handles.selFit(j).fitMethod='Automatic';
        handles.selFit(j).topMethod=handles.selFit(j).topMethodAuto;
        handles.selFit(j).botMethod=handles.selFit(j).botMethodAuto;
        handles.selFit(j).exponent=handles.selFit(j).exponentAuto;

        % Compute selected fit
        % --------------------
        selFit(j)=SelectFit(handles.normData(j),'Automatic',handles.selFit(j));
        
        % Update GUI
        % ----------
        if strcmp(handles.selFit(j).topMethodAuto,'Power')
            set(handles.puTop,'Value',1);
        else
            set(handles.puTop,'Value',2);
        end
        if strcmp(handles.selFit(j).botMethodAuto,'Power')
            set(handles.puBottom,'Value',1);
        else
            set(handles.puBottom,'Value',2);
        end
        set(handles.edExp,'String',num2str(handles.selFit(j).exponent,'%6.4f'));
        handles.selFit=selFit;
        plotTrans(handles); 
    else
        
        % Fit set to Manual
        % -----------------
        selFit(j).fitMethod='Manual';
        set(handles.puTop,'Enable','on');
        set(handles.puBottom,'Enable','on');
        set(handles.edExp,'Enable','on');
        set(handles.cbOpt,'Enable','on');   
        handles.selFit=selFit;
    end
    
     % Compute and display discharge sensitivity
    % -----------------------------------------
    if strcmpi(handles.normData(j).dataType,'q')
        sensData=QSensitivity(handles.transData,handles.selFit);
    else
        sensData=QSensitivity();
    end
    h=handles.q_table;
    DisplayTable(sensData,h);
    handles.sensData=sensData;
    
    % Update handles structure
    guidata(hObject, handles);  
    set(handles.figure1,'Pointer','arrow');
    drawnow;

%==========================================================================
function listFiles_Callback(hObject, eventdata, handles)
% 
% Allows the user to select which transect is plotted in the transect graph
%==========================================================================
   set(handles.figure1,'Pointer','watch');
   drawnow;
% Update graph
% ------------
plotTrans(handles)
%
% Update fit configuration for the transect
% -------------------------------------------------
selFit=handles.selFit;
itrans=get(handles.listFiles,'Value');
if strcmpi(selFit(itrans).fitMethod,'Automatic')
    set(handles.puFit,'Value',1);
    set(handles.puTop,'Enable','off');
    set(handles.puBottom,'Enable','off');
    set(handles.edExp,'Enable','off');
    set(handles.cbOpt,'Enable','off'); 
else
    set(handles.puFit,'Value',2);
    set(handles.puTop,'Enable','on');
    set(handles.puBottom,'Enable','on');
    set(handles.edExp,'Enable','on');
    set(handles.cbOpt,'Enable','on');
end
if strcmpi(selFit(itrans).topMethod,'Power')
    set(handles.puTop,'Value',1);
else
    set(handles.puTop,'Value',2);
end
if strcmpi(selFit(itrans).botMethod,'Power')
    set(handles.puBottom,'Value',1);
else
    set(handles.puBottom,'Value',2);
end
set(handles.edExp,'String',num2str(selFit(itrans).exponent,'%6.4f'));
if strcmpi(selFit(itrans).expMethod,'Optimize')
    set(handles.cbOpt,'Value',1)
else
    set(handles.cbOpt,'Value',0);
end
set(handles.figure1,'Pointer','arrow');
drawnow;

%**************************************************************************
%**************************************************************************
%                 No User Code Below Here
%**************************************************************************
%**************************************************************************

function File_Callback(hObject, eventdata, handles)
function Configure_Callback(hObject, eventdata, handles)
function View_Callback(hObject, eventdata, handles)
function DataType_Callback(hObject, eventdata, handles)
function varargout = extrap3_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Get default command line output from handles structure
varargout{1} = handles.output;
function trans_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listFiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function puFit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to puFit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function puTop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to puTop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function puBottom_CreateFcn(hObject, eventdata, handles)
% hObject    handle to puBottom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edExp_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function listFiles_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function [normData,selFit,sensData]=processProfiles(handles)
%
% This function process each transect and a composite data set to compute
% the normalized discharge or velocity and depth and then calls the object
% to automatically select an extrapolation method. The sensitivity to the
% various extrapolation methods to discharge is also computed.
%==========================================================================
%
% Retrieve data from handles data structure
% -----------------------------------------
nfiles=handles.nfiles;
fileName=handles.fileName;
threshold=handles.threshold;
selFit=handles.selFit;

% Initialize variables
% --------------------
% statsTable=cell(nfiles+1,8);
% statsTable{nfiles+1,1}='All';
transData=handles.transData;
%
% Determine data type to be processed
% ------------------------------------
if strcmp(get(handles.Discharge,'Checked'),'on')
    dataType='q';
else
    dataType='v';
end
normData(1:nfiles+1)=NormData();
%
% Begin processing each transect.
% ------------------------------
for ifile=1:nfiles
    %
    % Update stats table to show processing file
    % ------------------------------------------
%     statsTable(ifile,1)=fileName(ifile);
%     set(handles.fit_table,'Data',statsTable);
%     drawnow;
    %
    % Normalize data and create normData for class NormData
    % -------------------------------------------------------
    normData(ifile)=NormData(transData(ifile),dataType,threshold,handles.extents);
    %
    % Fit normalized data
    % --------------------------------------------------------
    if strcmpi(selFit(ifile).fitMethod,'Manual')
        selFit(ifile)=SelectFit(normData(ifile),selFit(ifile).fitMethod,selFit(ifile));
    else
        selFit(ifile)=SelectFit(normData(ifile),selFit(ifile).fitMethod,selFit(ifile));
    end


    set(handles.listFiles,'Value',ifile);
    drawnow;
end

% Create and empty object of class NormData for the composite profile
% -------------------------------------------------------------------
sumens(1)=0;
maxcells=0;
for ifile=1:nfiles
    bins(ifile)=size(normData(ifile).unitNormalized,1);
    if bins>maxcells
        maxcells=bins;
    end
    numEns(ifile)=size(normData(ifile).unitNormalized,2);
    sumens(ifile+1)=sumens(ifile)+numEns(ifile);
end

normData(nfiles+1).unitNormalized(1:maxcells,1:sumens)=nan;
normData(nfiles+1).cellDepthNormalized(1:maxcells,1:sumens)=nan;
normData(nfiles+1).fileName='Measurement';
normData(nfiles+1).dataType=upper(dataType);

% Build unitNormalized data property for composite profile
% ---------------------------------------------------------
for ifile=1:nfiles   
    normData(nfiles+1).unitNormalized(1:bins(ifile),1+sumens(ifile):sumens(ifile+1))=normData(ifile).unitNormalized;
    normData(nfiles+1).cellDepthNormalized(1:bins(ifile),1+sumens(ifile):sumens(ifile+1))=normData(ifile).cellDepthNormalized;
end
%idx=find(~isnan(normData(nfiles+1).unitNormalized));

%
% Compute remainder of properties for the composite profile
% ---------------------------------------------------------
normData(nfiles+1)=NormData(normData(nfiles+1),normData(nfiles+1).dataType,threshold,handles.extents);
%
% Fit composite profile
% ---------------------
ifile=nfiles+1;
if strcmpi(selFit(ifile).fitMethod,'Manual')
    selFit(ifile)=SelectFit(normData(ifile),selFit(ifile).fitMethod,selFit(ifile));
else
    selFit(ifile)=SelectFit(normData(ifile),selFit(ifile).fitMethod,selFit(ifile));
end

if (selFit(nfiles+1).topfitr2>0.9 || selFit(nfiles+1).topr2>0.9) && abs(selFit(nfiles+1).topmaxdiff)>0.2
    set(handles.msg1,'String','The measurement profile may warrant a 3-point fit at the top');
end

% Display automated fit selection
% -------------------------------
if strcmpi(selFit(nfiles+1).topMethod,'Power')
    set(handles.puTop,'Value',1);
else
    set(handles.puTop,'Value',2);
end
if strcmpi(selFit(nfiles+1).botMethod,'Power')
    set(handles.puBottom,'Value',1);
else
    set(handles.puBottom,'Value',2);
end
set(handles.edExp,'String',num2str(selFit(nfiles+1).exponent,'%6.4f'));
set(handles.listFiles,'Value',ifile);
drawnow;

% Compute and display discharge sensitivity
% -----------------------------------------
if strcmp(dataType,'q')
    sensData=QSensitivity(transData,selFit);
else
    sensData=QSensitivity();
end
h=handles.q_table;
DisplayTable(sensData,h);
% h=handles.fit_table;
% DisplayTable(selFit,h);


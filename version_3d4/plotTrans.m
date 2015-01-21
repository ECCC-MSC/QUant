function plotTrans(handles)
%
% This function updates the trans_axes of the GUI with the selected data and
% profile fits.
%
% Last modificaitons / validation 5/15/2012
%==========================================================================

% Get axes handle and clear axes
% -----------------------------
h=handles.axes_trans;
cla(h);

% Get data from handles data structure
% ------------------------------------
itrans=get(handles.listFiles,'Value');
threshold=handles.threshold;
nfiles=handles.nfiles;
selFit=handles.selFit;

% Set axes default title
% ----------------------
method='manual';

% Plot cell data
% --------------
if strcmp(get(handles.CellData,'Checked'),'on')
    plotRaw(handles.normData(itrans),h);
end

% Highlight surface cells
% -----------------------
if strcmp(get(handles.SurfCells,'Checked'),'on')
    plotSurface(handles.normData(itrans),h);
end
 
% Plot transect median points
% ---------------------------
if strcmp(get(handles.TransMed,'Checked'),'on') && itrans<nfiles+1
    plotMed(handles.normData(itrans),h);
elseif strcmp(get(handles.TransMed,'Checked'),'on') && itrans==nfiles+1
    for j=1:itrans
        hold on
        plotMed(handles.normData(j),h);
    end
end

% Plot measurement median points
% ------------------------------
if strcmp(get(handles.MeasMed,'Checked'),'on')
    plotMed(handles.normData(handles.nfiles+1),h);
end

% Plot transect profile fits
% --------------------------
if strcmp(get(handles.TransFit,'Checked'),'on')
    if itrans<nfiles+1
        plotfit(selFit(itrans),h,handles.transData(itrans).startBank);
    else
        for j=1:nfiles
            hold on
            plotfit(selFit(j),h,handles.transData(j).startBank);
        end
    end
            
end

% Plot composite profile fit
% --------------------------
if strcmp(get(handles.MeasFit,'Checked'),'on')
    plotfit(handles.selFit(handles.nfiles+1),h);
end

% Plot scale and labels
% ---------------------
plotScale(handles.normData(itrans),h);
plotLabel(handles.normData(itrans),h);

% Update number of points in each 5% increment
% --------------------------------------------
h=handles.zpts_table;
DisplayNos(handles.normData(itrans)',h);
drawnow;


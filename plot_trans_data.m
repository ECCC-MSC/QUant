function plot_trans_data(useData, navRef, show)

% navRef = 'gga', 'vtg' or 'bot'

% if you want to display file name, remove the semicolon

useData.filename;


% Assign variables from input data to local variables.
% ----------------------------------------------------
cellSize=useData.cellSize;
cellDepth=useData.cellDepth;

if navRef == 'bot'
    wVelx=useData.wVelx;
    wVely=useData.wVely;
    bVelx=useData.btVel(1,:);
    bVely=useData.btVel(2,:);
elseif navRef == 'gga'
    wVelx = useData.wVelx_reGGA;
    wVely = useData.wVely_reGGA;
    bVelx = useData.btVel_reGGA(1,:);
    bVely = useData.btVel_reGGA(2,:);
elseif navRef == 'vtg'
    wVelx = useData.wVelx_reVTG;
    wVely = useData.wVely_reVTG;
    bVelx = useData.btVel_reVTG(1,:);
    bVely = useData.btVel_reVTG(2,:);
end

depthEns=useData.depthEns;
ensDeltaTimeadj=useData.ensDeltaTime;
ensDeltaTime=useData.ensDeltaTime;

numEns=useData.numEns;
startEdge=useData.startBank;

leftDist = useData.lDist;
rightDist = useData.rDist;


leftCoef = useData.leftCoef;
leftNumEns2Avg = useData.leftNumEns2Avg;
rightCoef = useData.rightCoef; %
rightNumEns2Avg = useData.rightNumEns2Avg;



% Compute / initialize local variables
% ------------------------------------
k=0;
z=repmat(depthEns,size(cellDepth,1),1)-cellDepth;
z(isnan(wVelx))=nan;

% Loop to identify first and last valid cell in each ensemble,
% invalid ensembles, and invalid cells within valid ensembles.
% The range to the top of the top
% cell and bottom of the bottom cell is also computed.
% -----------------------------------------------------------
lgclIntQCells=true(size(cellSize));
idxTop=zeros(1,size(wVelx,2));
idxTop3=zeros(3,size(wVelx,2));
idxBot=zeros(1,size(wVelx,2));
idxValidEns=false(1,size(wVelx,2));
for jj=1:numEns
    test=find(~isnan(wVelx(:,jj)), 1);
    %
    % If entire ensemble has no valid data the ensemble is
    % marked invalid.
    % -----------------------------------------------------
    if ~isempty(test)
        idxValidEns(jj)=true;
        k=k+1;
        
        % Identify valid top and bottom cell and top 3 cells for
        % the 3-pt extrapolation method
        % ------------------------------------------------------
        idxTop(jj)=find(~isnan(wVelx(:,jj)),1,'first');
        temp=find(~isnan(wVelx(:,jj)),3,'first');
        if length(temp)>2
            idxTop3(:,jj)=temp;
        end
        idxBot(jj)=find(~isnan(wVelx(:,jj)),1,'last');
        
        % Identify invalid cells between top and bottom
        % ---------------------------------------------
        temp=find(isnan(wVelx(idxTop(jj):idxBot(jj),jj)))+idxTop(jj)-1;
        lgclIntQCells(temp,jj)=false;
        
        % Compute range for top and bottom extrapolation
        % ----------------------------------------------
        topRng(jj)=cellDepth(idxTop(jj),jj)-0.5.*cellSize(idxTop(jj),jj);
        
        botRng(jj)=depthEns(jj)-cellDepth(idxBot(jj),jj)-0.5*cellSize(idxBot(jj),jj);
        
    else
        
        topRng(jj)=nan;
        botRng(jj)=nan;
        
        % The delta time on the next valid ensemble is
        % increased by the delta time of the invalid ensemble
        % ---------------------------------------------------
        if jj+1<numEns && ~isnan(ensDeltaTimeadj(jj))
            ensDeltaTimeadj(jj+1)=ensDeltaTimeadj(jj)+ensDeltaTimeadj(jj+1);
        end
        ensDeltaTimeadj(jj)=nan;
    end
end

% Compute the cross product used in the other computations
% --------------------------------------------------------
xprod=((wVelx.*repmat(bVely,size(wVelx,1),1))-...
    (wVely.*repmat(bVelx,size(wVely,1),1)));




% Determine sign of cross product (see p 132 of Stephanie's
% first uOttawa lab book
% really, all you want is for the discharge to be positive in
% the end
%%%%%%%%%%%%


if strcmpi(startEdge,'Left')
    direc=1;
else
    direc=-1;
end
xprod=xprod.*direc;


length_refBT =  nancumsum(useData.ensDeltaTime.*(sqrt(useData.btVel(1,:).^2 + (useData.btVel(2,:).^2))));
vel2 = sqrt(wVelx.^2 + wVely.^2);

% Display the velocity data
% show the velocity magnitude with length of the track relative
% to the bottom as you would see it in Winriver II (note the
% colorbars are difference)

if show == 1
    figure
    hold on
    subplot 211
    hp = pcolor(length_refBT, cellDepth, vel2);
    hold on
    set(gca, 'ydir', 'rev')
    set(hp,'edgecolor','none')
    shading flat
    plot(length_refBT, useData.depthEns, 'k-.', 'linewidth', 1.5)
    xlabel(['Length (Ref: ', navRef, ') [m]'], 'fontsize', 16)
    ylabel('Depth [m]', 'fontsize', 16)
    set(gca, 'box', 'on')
    title('Velocity magnitude [m/s]')
    colorbar
    
end

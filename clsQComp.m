classdef clsQComp
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess=private)
        top
        middle
        bottom
        left
        right
        total
        intCells
        intEns
        
        navRef
    end
    
    methods
        function obj=clsQComp(transData,varargin)
        % Construct the clsQComp object
            if nargin>0
                % Determine number of transects
                nTransects=length(transData);
                % Preallocate class
                obj(nTransects)=clsQComp();
                
                % Process for each transect
                for n=1:nTransects
                    
                    % Assign variables from input data to local variables.
                    cellSize=transData(n).depths.depthCellSize_m;
                    cellDepth=transData(n).depths.depthCellDepth_m;
                    cellsAboveSL=transData(n).wVel.cellsAboveSL;
                    wVelx=transData(n).wVel.u_mps.*cellsAboveSL;
                    wVely=transData(n).wVel.v_mps.*cellsAboveSL;
                    bVelx=transData(n).btVel.u_mps;
                    bVely=transData(n).btVel.v_mps;
                    depthEns=transData(n).btDepths.Depth_m;
                    ensDeltaTimeadj=transData(n).sensors.ensDuration_sec;
                    ensDeltaTime=transData(n).sensors.ensDuration_sec;
                    %maxCells=size(wVelx,1);
                    numEns=size(wVelx,2);
                    startEdge=transData(n).startEdge;
                    leftDist=transData(n).edges.left.dist_m;
                    leftCoef=transData(n).edges.left.coef;
                    leftUserQ=transData(n).edges.left.userQ_cms;
                    leftNumEns2Avg=transData(n).edges.left.numEns2Avg;
                    rightDist=transData(n).edges.right.dist_m;
                    rightCoef=transData(n).edges.right.coef;
                    rightUserQ=transData(n).edges.right.userQ_cms;
                    rightNumEns2Avg=transData(n).edges.right.numEns2Avg;
                    obj(n).navRef=transData(n).wVel.navRef;
                    
                    
                    % Determine extrapolation methods and exponent
                    if ~isempty(varargin)
                        topMethod=varargin{1};
                        botMethod=varargin{2};
                        exponent=str2double(varargin{3});
                    else
                        topMethod=transData(n).extrap.topMethod;
                        botMethod=transData(n).extrap.botMethod;
                        exponent=transData(n).extrap.exponent;
                    end

                    % Compute / initialize local variables
                    k=0; 
                    z=repmat(depthEns,size(cellDepth,1),1)-cellDepth;
                    z(isnan(wVelx))=nan;

                    % Preallocate variables
                    lgclIntQCells=true(size(cellSize));
                    idxTop=zeros(1,size(wVelx,2));
                    idxTop3=zeros(3,size(wVelx,2));
                    idxBot=zeros(1,size(wVelx,2));
                    idxValidEns=false(1,size(wVelx,2));
                    topRng=nan(1,numEns);
                    botRng=nan(1,numEns);

                    % Loop to identify first and last valid cell in each ensemble,
                    % invalid ensembles, and invalid cells within valid ensembles.
                    % The range to the top of the top
                    % cell and bottom of the bottom cell is also computed.                    
                    for j=1:numEns
                        test=find(~isnan(wVelx(:,j)), 1);
                        
                        % If entire ensemble has no valid data the ensemble is
                        % marked invalid.
                        if ~isempty(test)
                           idxValidEns(j)=true;
                           k=k+1;
                           
                           % Identify valid top and bottom cell and top 3 cells for
                           % the 3-pt extrapolation method
                           idxTop(j)=find(~isnan(wVelx(:,j)),1,'first');
                           temp=find(~isnan(wVelx(:,j)),3,'first');
                           if length(temp)>2
                               idxTop3(:,j)=temp;
                           end
                           idxBot(j)=find(~isnan(wVelx(:,j)),1,'last');
                           
                           % Identify invalid cells between top and bottom
                           temp=find(isnan(wVelx(idxTop(j):idxBot(j),j)))+idxTop(j)-1;
                           lgclIntQCells(temp,j)=false;
                           
                           % Compute range for top and bottom extrapolation
                           topRng(j)=cellDepth(idxTop(j),j)-0.5.*cellSize(idxTop(j),j);
                           botRng(j)=depthEns(j)-cellDepth(idxBot(j),j)-0.5*cellSize(idxBot(j),j);

                        else

                            topRng(j)=nan;
                            botRng(j)=nan;

                            % The delta time on the next valid ensemble is
                            % increased by the delta time of the invalid ensemble
                            if j+1<numEns && ~isnan(ensDeltaTimeadj(j))
                                ensDeltaTimeadj(j+1)=ensDeltaTimeadj(j)+ensDeltaTimeadj(j+1);
                            end
                            ensDeltaTimeadj(j)=nan;
                        end
                    end

                    % Compute the cross product used in the other computations
                    xprod=((wVelx.*repmat(bVely,size(wVelx,1),1))-...
                        (wVely.*repmat(bVelx,size(wVely,1),1)));
                    
                    % Determine sign of cross product
                    if strcmp(startEdge,'Right')
                        dir=1;
                    else
                        dir=-1;
                    end
                    xprod=xprod.*dir;
                    
                    % Compute the measured or middle portion of the discharge using
                    % the delta time increase to account for invalid ensembles and
                    % the WinRiver II approach to computing discharge for invalid
                    % cells.
                    qMidCells=xprod.*cellSize.*repmat(ensDeltaTimeadj,size(wVely,1),1);
                    obj(n).middle=nansum(nansum(qMidCells));

                    % Compute the discharge for invalid cells
                     qIntCellsCells=obj.dischargeIntcells(botMethod,exponent,xprod,z,cellSize,depthEns,cellDepth,lgclIntQCells,ensDeltaTimeadj);
                     obj(n).intCells=nansum(nansum(qIntCellsCells));

                    % Update middle discharge
                    obj(n).middle=obj(n).middle+obj(n).intCells;

                    % Compute the top discharge
                    qTopEns=obj.dischargeTop(topRng,idxTop,idxTop3,topMethod,exponent,xprod,z,cellSize,cellDepth,depthEns,ensDeltaTimeadj);
                    obj(n).top=nansum(qTopEns);       

                    % Compute the bottom discharge
                    qBotEns=obj.dischargeBot(botRng,idxBot,botMethod,exponent,xprod,z,cellSize,depthEns,cellDepth,ensDeltaTimeadj);
                    obj(n).bottom=nansum(qBotEns);

                    % Compute right edge discharge
                    if isnan(rightUserQ)
                        obj(n).right=obj.dischargeEdge(wVelx,wVely,bVelx,bVely,depthEns,ensDeltaTimeadj,startEdge,'Right',rightCoef,rightDist,rightNumEns2Avg);
                    else
                        obj(n).right=rightUserQ;
                    end
                    
                    % Compute left edge discharge
                    if isnan(rightUserQ)
                        obj(n).left=obj.dischargeEdge(wVelx,wVely,bVelx,bVely,depthEns,ensDeltaTimeadj,startEdge,'Left',leftCoef,leftDist,leftNumEns2Avg);
                    else
                        obj(n).left=leftUserQ;
                    end
                    
                    % Compute total discharge
                    obj(n).total=obj(n).top+obj(n).middle+obj(n).bottom+obj(n).left+obj(n).right;

                    % Compute discharge for valid cells and ensembles only
                    qMidNIcells=xprod.*cellSize.*repmat(ensDeltaTime,size(wVely,1),1);
                    qMidNI=nansum(nansum(qMidNIcells));
                    qTopNIens=obj.dischargeTop(topRng,idxTop,idxTop3,topMethod,exponent,xprod,z,cellSize,cellDepth,depthEns,ensDeltaTime);
                    qTopNI=nansum(qTopNIens);
                    qBotNIens=obj.dischargeBot(botRng,idxBot,botMethod,exponent,xprod,z,cellSize,depthEns,cellDepth,ensDeltaTime);        
                    qBotNI=nansum(qBotNIens);

                    % Determine the discharge computed for invalid ensembles
                    obj(n).intEns=obj(n).middle+obj(n).top+obj(n).bottom-qMidNI-qTopNI-qBotNI;                   
                end % for n
            end % if nargin
        end % constructor      
    end % methods
    
    methods (Static)

        function qtop=dischargeTop(topRng,idxTop,idxTop3,topMethod,...
                exponent,xprod,z,cellSize,cellDepth,depthEns,deltat)
        % This function computes the extrapolated top discharge using
        % either constant or power or 3-pt methods.
        
            % Identify method
            switch topMethod
                
                case 'Power'                 
                    
                    coef=((exponent+1).*nansum(xprod.*cellSize))./...
                         nansum(((z+0.5.*cellSize).^(exponent+1))-...
                         ((z-0.5.*cellSize).^(exponent+1)));
                    qtop=deltat.*(coef./(exponent+1)).*...
                        (depthEns.^(exponent+1)-(depthEns-topRng).^(exponent+1));
                
                case 'Constant'
                    
                    nens=length(deltat);
                    qtop=nan(1,nens);
                    for j=1:nens
                        if idxTop(j)~=0
                            qtop(j)=deltat(j).*xprod(idxTop(j),j).*topRng(j);
                        end % if idxTop
                    end % for j
                    
                case '3-Point'
                    
                    % Determine number of bins available in each profile
                    nbins=ones(size(xprod));
                    nbins(isnan(xprod))=0;
                    nbins=nansum(nbins);
                    
                    % 3-pt is only applied to profiles with more than 6
                    % bins, otherwise constant is used
                    idxconstant=ones(size(nbins));
                    idxconstant(nbins>5)=0;
                    idx3pt=~idxconstant;
                    
                    % Compute constant fit where appropriate
                    deltatC=deltat;
                    deltatC(idx3pt)=nan;
                    nens=length(deltat);
                    qtopC=zeros(1,nens);
                    for j=1:nens
                        if idxTop(j)~=0
                            qtopC(j)=deltatC(j).*xprod(idxTop(j),j).*topRng(j);
                        end % if idxTop
                    end % for j
                    
                    % Compute 3-pt fit using linear least squares fit of 
                    % top 3 valid bins
                    qtop3=zeros(1,nens);                            
                    deltat3=deltat;
                    deltat3(logical(idxconstant))=0;
                    for j=1:nens
                        if idx3pt(j)==1
                            sumd=nansum(cellDepth(idxTop3(1:3,j),j));
                            sumd2=nansum(cellDepth(idxTop3(1:3,j),j).^2);
                            sumQ=nansum(xprod(idxTop3(1:3,j),j));
                            sumQd=nansum(xprod(idxTop3(1:3,j),j).*cellDepth(idxTop3(1:3,j),j));
                            delta=3*sumd2-sumd.^2;
                             A=(3.*sumQd-sumQ.*sumd)./delta;
                             B=(sumQ.*sumd2-sumQd.*sumd)./delta;
                    
                            % Compute discharge for 3-pt fit                      
                            Qo=(A.*topRng(j).^2)./2+B.*topRng(j);
                            qtop3(j)=deltat3(j).*Qo;
                        end % if idx3pt
                    end % for j
                    
                    % Combine constant and 3-pt discharges
                    qtop=nansum(qtop3)+nansum(qtopC);
            end % switch topMethod
        end % dischargeTop
        
        function qbot=dischargeBot(botRng,idxBot,botMethod,exponent,...
                xprod,z,cellSize,depthEns,cellDepth,deltat)
        % This function computes the extrapolated bottom discharge
        % using either the power or no slip method as defined in
        % WinRiver II.
        
            % Indentify bottom method
            switch botMethod
                
                case 'Power'
                    
                    coef=((exponent+1).*nansum(xprod.*cellSize))./...
                         nansum(((z+0.5.*cellSize).^(exponent+1))-((z-0.5.*cellSize).^(exponent+1)));
                    qbot=deltat.*(coef./(exponent+1)).*(botRng.^(exponent+1));    
                    
                case 'No Slip'
                    
                    % Valid data in the lower 20% of the water column or
                    % the last valid depth cell are used to compute the no
                    % slip power fit.
                    cutoffDepth=0.8.*depthEns;
                    depthOK=(cellDepth>repmat(cutoffDepth,size(cellDepth,1),1));
                    xprodOK=~isnan(xprod);
                    usens=depthOK.*xprodOK;
                    for j=1:length(deltat)
                        if idxBot(j)~=0
                            usens(idxBot(j),j)=1;
                        end % if idxBot
                    end % for j
                    usens(usens==0)=nan;
                    
                    % Create cross product and z arrays for the data to be
                    % used in no slip computations.
                    xprodns=usens.*xprod;
                    zns=usens.*z;
                    coef=((exponent+1).*nansum(xprodns.*cellSize))./...
                        nansum(((zns+0.5.*cellSize).^(exponent+1))-...
                        ((zns-0.5.*cellSize).^(exponent+1)));
                    
                    % Compute the bottom discharge of each profile
                    qbot=deltat.*(coef./(exponent+1)).*(botRng.^(exponent+1)); 
            end % switch botMethod
        end % function dischargeBot
        
        function edge=dischargeEdge(wVelx,wVely,bVelx,bVely,depthens,ensDeltaTimeadj,startEdge,edgeLoc,edgeCoef,edgeDist,edgeEns)
            
            % Find valid ensembles to average            
            temp=~isnan(wVelx);
            temp=nansum(temp,1);
            temp=temp~=0;
            idx=find(temp>0);
            temp=cumsum(temp);
            
            % Compute the average velocity and depth
            if strcmp(edgeLoc,startEdge)
                % Identify ensembles to average
                edgeIdx=find(temp==edgeEns);  
                % Average all cells for edge together
                xVel=wVelx(:,1:edgeIdx);
                yVel=wVely(:,1:edgeIdx);
                % Average all cells
                [~,edgeMag]=cart2pol(nanmean(xVel(1:end)),nanmean(yVel(1:end)));

                % Average velocity by bin and then ensemble
                wVelxAvg=nanmean(nanmean(wVelx(:,1:edgeIdx),2));
                wVelyAvg=nanmean(nanmean(wVely(:,1:edgeIdx),2));
                % Average depth on for those ensembles with valid
                % velocities
                depthAvg=nanmean(depthens(idx(1:edgeEns)));
            else
                % Identify ensembles to average
                edgeIdx=find((-1.*(temp-(temp(end)+1)))==edgeEns);
                % Average velocity by bin and then ensemble
                 wVelxAvg=nanmean(nanmean(wVelx(:,edgeIdx:end),2));
                 wVelyAvg=nanmean(nanmean(wVely(:,edgeIdx:end),2));

                % Average all cells together
                xVel=wVelx(:,edgeIdx:end);
                yVel=wVely(:,edgeIdx:end);
                % Average all cells
                [~,edgeMag]=cart2pol(nanmean(xVel(1:end)),nanmean(yVel(1:end)));
                % Average depth on for those ensembles with valid
                % velocities
                depthAvg=nanmean(depthens(idx(end-edgeEns+1:end)));
            end % if edgeLoc
            
            % Compute velocity direction and magnitude
            %[~, dsmedgeMag]=cart2pol(wVelxAvg,wVelyAvg);
            %edgeDir=rad2azdeg(edgeDir);
            
            % Compute the navigation direction and magnitude
            trackX=nancumsum(bVelx.*ensDeltaTimeadj);
            trackY=nancumsum(bVely.*ensDeltaTimeadj);
            if strcmp(startEdge,'Right')
                dir=1;
            else
                dir=-1;
            end
            xprod=(wVelxAvg.*trackY(end)-wVelyAvg.*trackX(end)).*dir;
            edge=edgeCoef.*depthAvg.*edgeMag.*edgeDist.*sign(xprod);
            
        end
        
        function qintcells=dischargeIntcells(botMethod,exponent,xprod,z,...
                           cellSize,depthEns,cellDepth,lgclIntQCells,deltat)
        % This function computes the discharge for the invalid cells using
        % the methods in WinRiver II. Power fit using the power fit
        % equation and no slip uses linear interpolation.
            
            % Compute cell range from streambed
            ncells=size(xprod,1);
            zadj=nan(size(z));
            zall=repmat(depthEns,size(cellDepth,1),1)-cellDepth;
            zadj(~lgclIntQCells)=zall(~lgclIntQCells);
            
            % Determine appropriate interpolation method
            switch botMethod
                
                case 'Power'
                    
                    coef=((exponent+1).*nansum(xprod.*cellSize))./...
                         nansum(((z+0.5.*cellSize).^(exponent+1))-((z-0.5.*cellSize).^(exponent+1)));
                     
                    qintcells=repmat(deltat,ncells,1).*...
                              repmat(coef,ncells,1).*...
                              cellSize.*zadj.^exponent;
                          
                case 'No Slip'
                    
                    % Initialize indices
                    xi=nan(size(xprod));
                    xi(~lgclIntQCells)=1;
                    qintcells=nan(size(xprod));
                    
                    % Use linear interpolation for invalid cells
                    for j=1:size(cellDepth,2)
                        idxY=~isnan(xprod(:,j));
                        idxX= ~isnan(xi(:,j));
                        if sum(idxX)>0 && sum(idxY)>1
                            qintcells(idxX,j)=interp1(cellDepth(idxY,j),xprod(idxY,j),cellDepth(idxX,j),'linear');
                        end % if
                    end % for j
                    
                    qintcells=qintcells.*cellSize.*repmat(deltat,ncells,1);
            end % switch botMethod
        end % function dischargeIntCells
        
        function avg=meanQ(obj,prop);
            nMax=length(obj);
            for n=1:nMax
                data(n)=obj(n).(prop);
            end
            avg=nanmean(data);
        end
        
        function cov=covQ(obj,prop)
            nMax=length(obj);
            for n=1:nMax
                data(n)=obj(n).(prop);
            end
            avg=nanmean(data);
            sd=nanstd(data);
            cov=abs(sd/avg);
        end
        
        function ranUncert95=ranUncertQ(obj,prop)
            nMax=length(obj);
            for n=1:nMax
                data(n)=obj(n).(prop);
            end
            avg=nanmean(data);
            sd=nanstd(data);
            cov=abs(sd/avg);
            ranUncert95=100.*abs(tinv(.025,nMax-1)).*cov./sqrt(nMax);
        end
            
            
    end % static methods 
    
end % class


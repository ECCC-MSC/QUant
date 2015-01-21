
classdef Discharge_sm
    %
    % Class definition for discharge computations. Constructor method
    % requires originData and fitData class variables. Currently methods
    % utilized in WinRiver II have been implemented, with the exception of
    % the 3-point top extrapolation method.
    % David S. Mueller 6/6/2011, Revised DSM 5/11/2012.
    %
    % Added 3-pt extrapolation, fixed several bugs, validated
    % Last modificaitons / validation 5/15/2012
    %
    properties(SetAccess='private')
        filename
        qTop         % discharge in top
        qMid         % discharge in middle
        qBot         % discharge in bottom
        qTot         % total discharge less edges
        qIntEns      % discharge interpolated for invalid ensembles
        qIntEns2     % if using ensembles before the bad ensembles as opposed to afterwards as is the standard
        qIntCells    % discharge interpolated for invalid cells
        qMeas        % discharge computed for valid cells
        left         % discharge in left edge
        right        % discharge in right edge
        navRef
        
    end
    
    methods
        % =================================================================
        function obj=Discharge_sm(useData, navRef, show, varargin)
            % show is for displyaing the measurements, you don't want it to do it each
            % time the monte carlo simulation is run
            
            % navRef = 'gga', 'vtg' or 'bot'
            
            
            % if you want to add uncertainty in the edge estimates
            if length(varargin) > 0
                varyLEV = varargin{1};
                varyREV = varargin{2};
                varyQEXTPTB = varargin{3};
                varyQINTENS = varargin{4};
                forQintUncert = varargin{5};
            else
                varyLEV = 0;
                varyREV = 0;
                varyQEXTPTB = 0;
                varyQINTENS = 0;
                forQintUncert = 0;
                
            end
            obj.navRef = navRef;
           
            
            if strcmp(useData.topMethod, 'Power');
                topMethod = 0;
            elseif strcmp(useData.topMethod, 'Constant');
                topMethod = 1;
            elseif strcmp(useData.topMethod, '3-point');
                topMethod = 2;
            end

            if strcmp(useData.botMethod, 'Power');
                botMethod = 0;
            elseif strcmp(useData.botMethod, 'No Slip');
                botMethod = 2;
            end

         
            % If NOT using extrap3.m
           % topMethod = useData.topMethod
           % botMethod = useData.botMethod
            
            exponent = useData.exponent;
            
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
            ensDeltaTimeadj = useData.ensDeltaTime;
            ensDeltaTimeadj2 = useData.ensDeltaTime;
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
                    if jj+1<numEns && ~isnan(ensDeltaTimeadj(jj)) % if before end of variable and
                        ensDeltaTimeadj(jj+1)=ensDeltaTimeadj(jj)+ensDeltaTimeadj(jj+1);
                    end
                    ensDeltaTimeadj(jj)=nan;
                end
            end
            
  
            
            % want to swap places of the larger ensemble time steps from
            % after the missing ensembles to before, a kind of fliplr for
            % the value surrounding NaNs.
            
            % Compute the cross product used in the other computations
            % --------------------------------------------------------
            xprod = ((wVelx.*repmat(bVely,size(wVelx,1),1))-...
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
            
            % in case there's something wrong with the mmt file and it
            % can't read the start and end bank correctly, this little
            % cheating works
          
            if nanmean(nanmean(xprod)) < 0 % the wrong start bank was used so the velocities have the wrong sign
                
                xprod = -xprod;
            end
            
           
            
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
                
                subplot 212
                hp = pcolor(length_refBT, cellDepth, xprod);
                hold on
                set(gca, 'ydir', 'rev')
                set(hp,'edgecolor','none')
                shading flat
                plot(length_refBT, useData.depthEns, 'k-.', 'linewidth', 1.5)
                xlabel(['Length (Ref: ', navRef, ') [m]'], 'fontsize', 16)
                ylabel('Depth [m]', 'fontsize', 16)
                set(gca, 'box', 'on')
                title('Cross product')
                colorbar
                
            end
            
            
            %%
            
            %===========================
            % Compute the measured or middle portion of the discharge using
            % the delta time increase to account for invalid ensembles and
            % the WinRiver II approach to computing discharge for invalid
            % cells.
            % -------------------------------------------------------------
            
           
            
            
            qMidCells=xprod.*cellSize.*repmat(ensDeltaTimeadj,size(wVely,1),1);
            obj.qMid=nansum(nansum(qMidCells));
            
            %========================
            % IF YOU HAVE MULTIPLE TRANSECTS
            % obj(n).qMid=nansum(nansum(qMidCells));
            
            % Compute the discharge for invalid cells
            % ---------------------------------------
            
            qIntCellsCells=obj.dischargeIntcells(botMethod,exponent,xprod,z,cellSize,depthEns,cellDepth,lgclIntQCells,ensDeltaTimeadj);
            obj.qIntCells=nansum(nansum(qIntCellsCells));
            %
            
            
            % Update middle discharge
            % -----------------------
            obj.qMid=obj.qMid+obj.qIntCells;
            
            
            % Compute the top discharge
            % -------------------------
            qTopEns=obj.dischargeTop(topRng,idxTop,idxTop3,topMethod,exponent,xprod,z,cellSize,cellDepth,depthEns,ensDeltaTimeadj);
            obj.qTop = (1 + varyQEXTPTB*random('norm', 0, useData.QextrapTBErrPct/100, 1,1))*nansum(qTopEns);
            
            % Compute the bottom discharge
            % ----------------------------
            qBotEns=obj.dischargeBot(botRng,idxBot,botMethod,exponent,xprod,z,cellSize,depthEns,cellDepth,ensDeltaTimeadj);
      
          
            obj.qBot = (1 + varyQEXTPTB*random('norm', 0, useData.QextrapTBErrPct/100, 1,1))*nansum(qBotEns);
            
            
            % Compute right edge discharge
            [obj.right] = obj.dischargeEdge(wVelx,wVely,bVelx,bVely,depthEns,ensDeltaTimeadj,startEdge,'Right',rightCoef,rightDist,rightNumEns2Avg, varyREV);
            
            
            % Compute left edge discharge
            % have added uncertainty to it
            [obj.left] = obj.dischargeEdge(wVelx,wVely,bVelx,bVely,depthEns,ensDeltaTimeadj,startEdge,'Left',leftCoef,leftDist,leftNumEns2Avg, varyLEV);
            
            
            % Compute total discharge
            
            obj.qTot = obj.qMid + obj.qTop + obj.qBot + obj.left + obj.right;
            
            
            
            % Compute discharge for valid cells and ensembles only
            % ----------------------------------------------------
            qMidNIcells=xprod.*cellSize.*repmat(ensDeltaTime,size(wVely,1),1);
            qMidNI=nansum(nansum(qMidNIcells));
            qTopNIens=obj.dischargeTop(topRng,idxTop,idxTop3,topMethod,exponent,xprod,z,cellSize,cellDepth,depthEns,ensDeltaTime);
            qTopNI=nansum(qTopNIens);
            qBotNIens=obj.dischargeBot(botRng,idxBot,botMethod,exponent,xprod,z,cellSize,depthEns,cellDepth,ensDeltaTime);
            qBotNI=nansum(qBotNIens);
            
            % Determine the discharge computed for invalid ensembles
            % ------------------------------------------------------
            obj.qIntEns=obj.qMid + obj.qTop + obj.qBot - qMidNI - qTopNI - qBotNI;
            
            
            if isempty(varargin) % if you're running it without distributions to get the reference discharge as winriver would calculate it
                
                % Compute the discharge for invalid cells if you use the valid
                % ensemble before the missing ensembles instead of the ones
                % after
                % - Simply flip things around in the adjusted time interval matric
                
                ensDeltaTimeadj2 = ensDeltaTimeadj;
                
                tstep = nanmean(ensDeltaTime);
                ffactor = round(ensDeltaTimeadj./tstep);
                [indUse] = find(ffactor>1);
                nEnsInt = ffactor(indUse)-1;
                
                for ii = 1:length(indUse)
                    ensDeltaTimeadj2(indUse(ii) - nEnsInt(ii)-1:indUse(ii)) = fliplr(ensDeltaTimeadj2(indUse(ii) - nEnsInt(ii)-1:indUse(ii)));
                end
                
                qIntCellsCells2=obj.dischargeIntcells(botMethod,exponent,xprod,z,cellSize,depthEns,cellDepth,lgclIntQCells,ensDeltaTimeadj2);
                qIntCells2=nansum(nansum(qIntCellsCells2));
                
                qMidCells2 = xprod.*cellSize.*repmat(ensDeltaTimeadj2,size(wVely,1),1);
                qMid2 = nansum(nansum(qMidCells2));
                
                qTop2 = nansum(obj.dischargeTop(topRng,idxTop,idxTop3,topMethod,exponent,xprod,z,cellSize,cellDepth,depthEns,ensDeltaTimeadj2));
                qBot2 = nansum(obj.dischargeBot(botRng,idxBot,botMethod,exponent,xprod,z,cellSize,depthEns,cellDepth,ensDeltaTimeadj2));
                
                obj.qIntEns2 = qMid2 + qIntCells2 + qTop2 + qBot2 - qMidNI - qTopNI - qBotNI;
            end
            
            
            
            obj.qTot = obj.qTot + varyQINTENS*obj.qIntEns*random('norm', 0, forQintUncert,1,1);
            
          
            % forQintUncert is the uncertainty as a fraction of qIntEns, 
            % it is calculated as  (Qnodist_reBot.qIntEns -  Qnodist_reBot.qIntEns2)/ Qnodist_reBot.qIntEns
            
            % Assign variables to object properties
            % -------------------------------------
            obj.filename=useData.filename;
            
            obj.qMeas=qMidNI;
            
        end
    end
    methods (Static)
        %==================================================================
        function qtop=dischargeTop(topRng,idxTop,idxTop3,topMethod,exponent,xprod,z,cellSize,cellDepth,depthEns,deltat)
            %
            % This function computes the extrapolated top discharge using
            % either constant or power or 3-pt methods.
            %==================================================================
            switch topMethod
                
                case 0 %'Power'
                    
                    coef=((exponent+1).*nansum(xprod.*cellSize))./...
                        nansum(((z+0.5.*cellSize).^(exponent+1))-((z-0.5.*cellSize).^(exponent+1)));
                    qtop=deltat.*(coef./(exponent+1)).*(depthEns.^(exponent+1)-(depthEns-topRng).^(exponent+1));
                    
                case 1 %'Constant'
                    
                    nens=length(deltat);
                    qtop=nan(1,nens);
                    for jj=1:nens
                        if idxTop(jj)~=0
                            qtop(jj)=deltat(jj).*xprod(idxTop(jj),jj).*topRng(jj);
                        end
                    end
                    
                case 2 %'3-Point'
                    
                    % Determine number of bins available in each profile
                    % --------------------------------------------------
                    nbins=ones(size(xprod));
                    nbins(isnan(xprod))=0;
                    nbins=nansum(nbins);
                    
                    % 3-pt is only applied to profiles with more than 6
                    % bins, otherwise constant is used
                    % -------------------------------------------------
                    idxconstant=ones(size(nbins));
                    idxconstant(nbins>5)=0;
                    idx3pt=~idxconstant;
                    
                    % Compute constant fit where appropriate
                    % --------------------------------------
                    deltatC=deltat;
                    deltatC(idx3pt)=nan;
                    nens=length(deltat);
                    qtopC=zeros(1,nens);
                    for jj=1:nens
                        if idxTop(jj)~=0
                            qtopC(jj)=deltatC(jj).*xprod(idxTop(jj),jj).*topRng(jj);
                        end
                    end
                    
                    % Compute 3-pt fit using linear least squares fit of
                    % top 3 valid bins
                    % ---------------------------------------------------
                    qtop3=zeros(1,nens);
                    deltat3=deltat;
                    deltat3(logical(idxconstant))=0;
                    for jj=1:nens
                        if idx3pt(jj)==1
                            sumd=nansum(cellDepth(idxTop3(1:3,jj),jj));
                            sumd2=nansum(cellDepth(idxTop3(1:3,jj),jj).^2);
                            sumQ=nansum(xprod(idxTop3(1:3,jj),jj));
                            sumQd=nansum(xprod(idxTop3(1:3,jj),jj).*cellDepth(idxTop3(1:3,jj),jj));
                            delta=3*sumd2-sumd.^2;
                            A=(3.*sumQd-sumQ.*sumd)./delta;
                            B=(sumQ.*sumd2-sumQd.*sumd)./delta;
                            
                            % Compute discharge for 3-pt fit
                            % ------------------------------
                            Qo=(A.*topRng(jj).^2)./2+B.*topRng(jj);
                            qtop3(jj)=deltat3(jj).*Qo;
                        end
                    end
                    
                    % Combine constant and 3-pt discharges
                    % ------------------------------------
                    qtop=nansum(qtop3)+nansum(qtopC);
            end
        end
        
        %==================================================================
        function qbot=dischargeBot(botRng,idxBot,botMethod,exponent,xprod,z,cellSize,depthEns,cellDepth,deltat)
            %
            % This function computes the extrapolated bottom discharge
            % using either the power or no slip method as defined in
            % WinRiver II.
            %==================================================================
            
            switch botMethod
                
                case 0 %'Power'
                    
                    coef=((exponent+1).*nansum(xprod.*cellSize))./...
                        nansum(((z+0.5.*cellSize).^(exponent+1))-((z-0.5.*cellSize).^(exponent+1)));
                    qbot=deltat.*(coef./(exponent+1)).*(botRng.^(exponent+1));
                    
                case 2 %'No Slip'
                    
                    % Valid data in the lower 20% of the water column or
                    % the last valid depth cell are used to compute the no
                    % slip power fit.
                    % ----------------------------------------------------
                    cutoffDepth=0.8.*depthEns;
                    depthOK=(cellDepth>repmat(cutoffDepth,size(cellDepth,1),1));
                    xprodOK=~isnan(xprod);
                    usens=depthOK.*xprodOK;
                    for jj=1:length(deltat)
                        if idxBot(jj)~=0
                            usens(idxBot(jj),jj)=1;
                        end
                    end
                    usens(usens==0)=nan;
                    
                    % Create cross product and z arrays for the data to be
                    % used in no slip computations.
                    % ----------------------------------------------------
                    xprodns=usens.*xprod;
                    zns=usens.*z;
                    coef=((exponent+1).*nansum(xprodns.*cellSize))./...
                        nansum(((zns+0.5.*cellSize).^(exponent+1))-((zns-0.5.*cellSize).^(exponent+1)));
                    
                    % Compute the bottom discharge of each profile
                    % --------------------------------------------
                    qbot=deltat.*(coef./(exponent+1)).*(botRng.^(exponent+1));
            end
        end
        
        %==================================================================
        
        %%
        function [edge]=dischargeEdge(wVelx,wVely,bVelx,bVely,depthEns,ensDeltaTimeadj,startEdge,edgeLoc,edgeCoef,edgeDist,edgeEns,varyEdgeVel)
            % have modified D. Mueller's code so that uncertainty is added to the
            % velocity used in the edge estimate. This uncertainty is calculated
            % from the standard deviation of the velocities use for the edge
            % estimate
            
            % Find valid ensembles to average
            temp=~isnan(wVelx);
            temp=nansum(temp,1);
            temp=temp~=0;
            idx=find(temp>0);
            temp=cumsum(temp);
            
            % Compute the average velocity and depth
            if strcmpi(edgeLoc,startEdge)
                % Identify ensembles to average
                edgeIdx=find(temp==edgeEns);
                % Average all cells for edge together
                xVel=wVelx(:,1:edgeIdx);
                yVel=wVely(:,1:edgeIdx);
                % Average all cells
                [~,edgeMag]=cart2pol(nanmean(xVel(1:end)),nanmean(yVel(1:end)));
                [~,edgeStd] = cart2pol(nanstd(xVel(1:end)),nanstd(yVel(1:end)));
                % SAM 2013-09-26: edgeStd is the standard deviation of
                % all measurements used for the edge veloctiy estimates
                
                % Average velocity by bin and then ensemble
                wVelxAvg=nanmean(nanmean(wVelx(:,1:edgeIdx),2));
                %wVelxStd=nanstd(nanstd(wVelx(:,1:edgeIdx),0,2));
                wVelyAvg=nanmean(nanmean(wVely(:,1:edgeIdx),2));

                % Average depth on for those ensembles with valid
                % velocities
              
                depthAvg=nanmean(depthEns(idx(1:edgeEns)));
                
                
            else
                % Identify ensembles to average
                edgeIdx=find((-1.*(temp-(temp(end)+1)))==edgeEns);
                % Average velocity by bin and then ensemble
                wVelxAvg=nanmean(nanmean(wVelx(:,edgeIdx:end),2));
                % wVelxStd=nanstd(nanstd(wVelx(:,edgeIdx:end),0,2));
                wVelyAvg=nanmean(nanmean(wVely(:,edgeIdx:end),2));
                
                
                % Average all cells together
                xVel=wVelx(:,edgeIdx:end);
                yVel=wVely(:,edgeIdx:end);
                
                % Average all cells
                [~,edgeMag] = cart2pol(nanmean(xVel(1:end)),nanmean(yVel(1:end)));
                [~,edgeStd] = cart2pol(nanstd(xVel(1:end)),nanstd(yVel(1:end)));
                
                % Average depth of those ensembles with valid velocities
               
                depthAvg=nanmean(depthEns(idx(end-edgeEns+1:end)));
               
                
            end % if edgeLoc
            
            % Compute velocity direction and magnitude
            %[~, dsmedgeMag]=cart2pol(wVelxAvg,wVelyAvg);
            %edgeDir=rad2azdeg(edgeDir);
            
            % Compute the navigation direction and magnitude
            trackX=nancumsum(bVelx.*ensDeltaTimeadj);
            trackY=nancumsum(bVely.*ensDeltaTimeadj);
            
            
            %see SMoore's lab book page 132
            %=====================
            if strcmpi(startEdge,'Left')
                direc=1;
            else
                direc=-1;
            end
            %====================
            
            xprod=(wVelxAvg.*trackY(end)-wVelyAvg.*trackX(end)).*direc;
            
            % cheat in case you got the wrong direction, this is assuming
            % that the average velocity vector for the edge goes in the same
            % direction as the flow direction
            %%tellStephAveXprodUsedForEdge = nanmean(nanmean(xprod))
            %keyboard
            if nanmean(nanmean(xprod)) < 0 % the wrong start bank was used so the velocities have the wrong sign
                xprod = -xprod;
            end
            
            % add uncertainty to the velocity used for the edge estimate
            % edgeStd is a percentage of the mean
            edge=edgeCoef.*depthAvg.*edgeMag.*edgeDist.*sign(xprod).*(1 + varyEdgeVel*random('norm', 0, edgeStd,1, 1));
            
        end
        
        %==================================================================
        %%
        
        function qintcells=dischargeIntcells(botMethod,exponent,xprod,z,...
                cellSize,depthEns,cellDepth,lgclIntQCells,deltat)
            %
            % This function computes the discharge for the invalid cells using
            % the methods in WinRiver II. Power fit using the power fit
            % equation and no slip uses linear interpolation.
            %==================================================================
            
            % Compute cell range from streambed
            % ---------------------------------
            
            ncells=size(xprod,1);
            zadj=nan(size(z));
            zall=repmat(depthEns,size(cellDepth,1),1)-cellDepth;
            zadj(~lgclIntQCells)=zall(~lgclIntQCells);
            
            switch botMethod
                
                case  0 %'Power'
                    
                    coef=((exponent+1).*nansum(xprod.*cellSize))./...
                        nansum(((z+0.5.*cellSize).^(exponent+1))-((z-0.5.*cellSize).^(exponent+1)));
                    
                    qintcells=repmat(deltat,ncells,1).*...
                        repmat(coef,ncells,1).*...
                        cellSize.*zadj.^exponent;
                    
                case 2 %'No Slip'
                    
                    % Initialize indices
                    % ------------------
                    xi=nan(size(xprod));
                    xi(~lgclIntQCells)=1;
                    qintcells=nan(size(xprod));
                    
                    % Use linear interpolation for invalid cells
                    % ------------------------------------------
                    for jj=1:size(cellDepth,2)
                        idxY=~isnan(xprod(:,jj));
                        idxX= ~isnan(xi(:,jj));
                        if sum(idxX)>0 && sum(idxY)>1
                            qintcells(idxX,jj)=interp1(cellDepth(idxY,jj),xprod(idxY,jj),cellDepth(idxX,jj),'linear');
                        end
                    end
                    
                    qintcells=qintcells.*cellSize.*repmat(deltat,ncells,1);
            end
        end
        
    end
    
end
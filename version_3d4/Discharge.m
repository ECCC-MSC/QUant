classdef Discharge
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
       qIntCells    % discharge interpolated for invalid cells
       qMeas        % discharge computed for valid cells
       qEdL         % discharge in left edge (not implemented)
       qEdR         % discharge in right edge (not implemented)
       topMethod    % top extrapolation method
       botMethod    % bottom extrapolation method
       exponent     % exponent for extrapolation method
       lShape       % shape of left edge (not implemented)
       rShape       % shape of right edge (not implemented)
   end
     
    methods
        % =================================================================
        function obj=Discharge(transData,topMethod,botMethod,exponent)
        %
        % Constructor method
        % =================================================================
            obj.qEdL=0;
            obj.qEdR=0;
            obj.lShape=nan;
            obj.rShape=nan;
            transData.filename
            
            % Assign variables from input data to local variables.
            % ----------------------------------------------------
            cellSize=transData.cellSize;
            cellDepth=transData.cellDepth;
            wVelx=transData.wVelx;
            wVely=transData.wVely;
            bVelx=transData.btVel(1,:);
            bVely=transData.btVel(2,:);
            depthEns=transData.depthEns;
            ensDeltaTimeadj=transData.ensDeltaTime;
            ensDeltaTime=transData.ensDeltaTime;
            maxCells=transData.maxCells;
            numEns=transData.numEns;
             
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
            for j=1:numEns
                test=find(~isnan(wVelx(:,j)), 1);
                %
                % If entire ensemble has no valid data the ensemble is
                % marked invalid.
                % -----------------------------------------------------
                if ~isempty(test)
                   idxValidEns(j)=true;
                   k=k+1;
                  
                   % Identify valid top and bottom cell and top 3 cells for
                   % the 3-pt extrapolation method
                   % ------------------------------------------------------
                   idxTop(j)=find(~isnan(wVelx(:,j)),1,'first');
                   temp=find(~isnan(wVelx(:,j)),3,'first');
                   if length(temp)>2
                       idxTop3(:,j)=temp;
                   end
                   idxBot(j)=find(~isnan(wVelx(:,j)),1,'last');
                   
                   % Identify invalid cells between top and bottom
                   % ---------------------------------------------
                   temp=find(isnan(wVelx(idxTop(j):idxBot(j),j)))+idxTop(j)-1;
                   lgclIntQCells(temp,j)=false;
                   
                   % Compute range for top and bottom extrapolation
                   % ----------------------------------------------
                   topRng(j)=cellDepth(idxTop(j),j)-0.5.*cellSize(idxTop(j),j);
                   botRng(j)=depthEns(j)-cellDepth(idxBot(j),j)-0.5*cellSize(idxBot(j),j);
                
                else
                    
                    topRng(j)=nan;
                    botRng(j)=nan;
                   
                    % The delta time on the next valid ensemble is
                    % increased by the delta time of the invalid ensemble
                    % ---------------------------------------------------
                    if j+1<numEns && ~isnan(ensDeltaTimeadj(j))
                        ensDeltaTimeadj(j+1)=ensDeltaTimeadj(j)+ensDeltaTimeadj(j+1);
                    end
                    ensDeltaTimeadj(j)=nan;
                end
            end
             
            % Compute the cross product used in the other computations
            % --------------------------------------------------------
            xprod=((wVelx.*repmat(bVely,size(wVelx,1),1))-...
                (wVely.*repmat(bVelx,size(wVely,1),1)));
            
            % Compute the measured or middle portion of the discharge using
            % the delta time increase to account for invalid ensembles and
            % the WinRiver II approach to computing discharge for invalid
            % cells.
            % -------------------------------------------------------------
            qMidCells=xprod.*cellSize.*repmat(ensDeltaTimeadj,size(wVely,1),1);
            obj.qMid=nansum(nansum(qMidCells));
            
            % Compute the discharge for invalid cells
            % ---------------------------------------
             qIntCellsCells=obj.dischargeIntcells(botMethod,exponent,xprod,z,cellSize,depthEns,cellDepth,lgclIntQCells,ensDeltaTimeadj);
             obj.qIntCells=nansum(nansum(qIntCellsCells));
            
            % Update middle discharge
            % -----------------------
            obj.qMid=obj.qMid+obj.qIntCells;
            
            % Compute the top discharge
            % -------------------------
            qTopEns=obj.dischargeTop(topRng,idxTop,idxTop3,topMethod,exponent,xprod,z,cellSize,cellDepth,depthEns,ensDeltaTimeadj);
            obj.qTop=nansum(qTopEns);       
            
            % Compute the bottom discharge
            % ----------------------------
            qBotEns=obj.dischargeBot(botRng,idxBot,botMethod,exponent,xprod,z,cellSize,depthEns,cellDepth,ensDeltaTimeadj);
            obj.qBot=nansum(qBotEns);
            
            % Compute total discharge
            % -----------------------
            obj.qTot=obj.qTop+obj.qMid+obj.qBot+obj.qEdL+obj.qEdR;
            
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
            obj.qIntEns=obj.qMid+obj.qTop+obj.qBot-qMidNI-qTopNI-qBotNI;
            
            % Assign variables to object properties
            % -------------------------------------
            obj.filename=transData.filename;
            obj.topMethod=topMethod;
            obj.botMethod=botMethod;
            obj.exponent=exponent;
            obj.lShape=nan;
            obj.rShape=nan;
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
                
                case 'Power'                 
                    
                    coef=((exponent+1).*nansum(xprod.*cellSize))./...
                         nansum(((z+0.5.*cellSize).^(exponent+1))-((z-0.5.*cellSize).^(exponent+1)));
                    qtop=deltat.*(coef./(exponent+1)).*(depthEns.^(exponent+1)-(depthEns-topRng).^(exponent+1));
                
                case 'Constant'
                    
                    nens=length(deltat);
                    qtop=nan(1,nens);
                    for j=1:nens
                        if idxTop(j)~=0
                            qtop(j)=deltat(j).*xprod(idxTop(j),j).*topRng(j);
                        end
                    end
                    
                case '3-Point'
                    
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
                    for j=1:nens
                        if idxTop(j)~=0
                            qtopC(j)=deltatC(j).*xprod(idxTop(j),j).*topRng(j);
                        end
                    end
                    
                    % Compute 3-pt fit using linear least squares fit of 
                    % top 3 valid bins
                    % ---------------------------------------------------
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
                            % ------------------------------
                            Qo=(A.*topRng(j).^2)./2+B.*topRng(j);
                            qtop3(j)=deltat3(j).*Qo;
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
                
                case 'Power'
                    
                    coef=((exponent+1).*nansum(xprod.*cellSize))./...
                         nansum(((z+0.5.*cellSize).^(exponent+1))-((z-0.5.*cellSize).^(exponent+1)));
                    qbot=deltat.*(coef./(exponent+1)).*(botRng.^(exponent+1));    
                    
                case 'No Slip'
                    
                    % Valid data in the lower 20% of the water column or
                    % the last valid depth cell are used to compute the no
                    % slip power fit.
                    % ----------------------------------------------------
                    cutoffDepth=0.8.*depthEns;
                    depthOK=(cellDepth>repmat(cutoffDepth,size(cellDepth,1),1));
                    xprodOK=~isnan(xprod);
                    usens=depthOK.*xprodOK;
                    for j=1:length(deltat)
                        if idxBot(j)~=0
                            usens(idxBot(j),j)=1;
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
                
                case 'Power'
                    
                    coef=((exponent+1).*nansum(xprod.*cellSize))./...
                         nansum(((z+0.5.*cellSize).^(exponent+1))-((z-0.5.*cellSize).^(exponent+1)));
                     
                    qintcells=repmat(deltat,ncells,1).*...
                              repmat(coef,ncells,1).*...
                              cellSize.*zadj.^exponent;
                          
                case 'No Slip'
                    
                    % Initialize indices
                    % ------------------
                    xi=nan(size(xprod));
                    xi(~lgclIntQCells)=1;
                    qintcells=nan(size(xprod));
                    
                    % Use linear interpolation for invalid cells
                    % ------------------------------------------
                    for j=1:size(cellDepth,2)
                        idxY=~isnan(xprod(:,j));
                        idxX= ~isnan(xi(:,j));
                        if sum(idxX)>0 && sum(idxY)>1
                            qintcells(idxX,j)=interp1(cellDepth(idxY,j),xprod(idxY,j),cellDepth(idxX,j),'linear');
                        end
                    end
                    
                    qintcells=qintcells.*cellSize.*repmat(deltat,ncells,1);
            end
        end
            
    end   
    
end
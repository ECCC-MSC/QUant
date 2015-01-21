classdef QSensitivity
    %
    % Class definition for discharge sensitivity computations. Constructor method
    % requires originData and fitData class variables. Currently methods
    % utilized in WinRiver II have been implemented, with the exception of
    % the 3-point top extrapolation method.
    % David S. Mueller 6/14/2011
    %
    %
    % Last modificaitons / validation 5/15/2012
    
    properties(SetAccess='private')
       qPPmean         % discharge power power 1/6
       qPPoptmean      % discharge power power optimized
       qCNSmean        % discharge constant no slip
       qCNSoptmean     % discharge constant optimized no slip
       q3pNSmean       % discharge 3-pt no slip
       q3pNSoptmean    % discharge 3-pt optimized no slip
       qPPoptperdiff   % power power optimized percent difference from power power 1/6     
       qCNSperdiff     % constant no slip percent difference from power power 1/6
       qCNSoptperdiff  % constant optimized no slip percent difference from power power 1/6
       q3pNSperdiff    % 3-point no slip percent difference from power power 1/6
       q3pNSoptperdiff % 3-point optimized no slip precent difference from power power 1/6
       ppExponent      % optimized power power exponent
       nsExponent      % optimized no slip exponent
       manTop
       manBot
       manExp
       qManmean
       qManperdiff
    end
     
    methods
        % =================================================================
        function obj=QSensitivity(transData,fitData)
        %
        % Constructor method
        % =================================================================        
        if nargin>0
            nfiles=size(transData,2);
            obj.ppExponent=fitData(1,nfiles+1).ppexponent;
            obj.nsExponent=fitData(1,nfiles+1).nsexponent;
            %
            % Compute discharge sensitivity
            % -----------------------------
            for j=1:nfiles
                qPP(j)=Discharge(transData(j),'Power','Power',0.1667);
                qPPopt(j)=Discharge(transData(j),'Power','Power',obj.ppExponent);
                qCNS(j)=Discharge(transData(j),'Constant','No Slip',0.1667);
                qCNSopt(j)=Discharge(transData(j),'Constant','No Slip',obj.nsExponent);
                q3PNS(j)=Discharge(transData(j),'3-Point','No Slip',0.1667);
                q3PNSopt(j)=Discharge(transData(j),'3-Point','No Slip',obj.nsExponent); 
                if strcmpi(fitData(1,nfiles+1).fitMethod,'Manual')
                    obj.manTop=fitData(1,nfiles+1).topMethod;
                    obj.manBot=fitData(1,nfiles+1).botMethod;
                    obj.manExp=fitData(1,nfiles+1).exponent;
                    qMan(j)=Discharge(transData(j),obj.manTop,obj.manBot,obj.manExp);
                end
            end
            %
            % Compute mean discharges
            % -----------------------
            obj.qPPmean=nanmean(abs([qPP.qTot]));
            obj.qPPoptmean=nanmean(abs([qPPopt.qTot]));
            obj.qCNSmean=nanmean(abs([qCNS.qTot]));
            obj.qCNSoptmean=nanmean(abs([qCNSopt.qTot]));
            obj.q3pNSmean=nanmean(abs([q3PNS.qTot]));
            obj.q3pNSoptmean=nanmean(abs([q3PNSopt.qTot]));  
            if strcmpi(fitData(1,nfiles+1).fitMethod,'Manual')
                obj.qManmean=nanmean(abs([qMan.qTot]));
            end
            %
            % Compute percent difference
            % --------------------------
            obj.qPPoptperdiff=((obj.qPPoptmean-obj.qPPmean)./obj.qPPmean).*100;
            obj.qCNSperdiff=((obj.qCNSmean-obj.qPPmean)./obj.qPPmean).*100; 
            obj.qCNSoptperdiff=((obj.qCNSoptmean-obj.qPPmean)./obj.qPPmean).*100;
            obj.q3pNSperdiff=((obj.q3pNSmean-obj.qPPmean)./obj.qPPmean).*100; 
            obj.q3pNSoptperdiff=((obj.q3pNSoptmean-obj.qPPmean)./obj.qPPmean).*100;
            if strcmpi(fitData(1,nfiles+1).fitMethod,'Manual')
                obj.qManperdiff=((obj.qManmean-obj.qPPmean)./obj.qPPmean).*100;
            end
        end
        end
        % =================================================================
        function Display(QSens)
        %
        % Display discharge sensitivity in separate window.
        % NOT CURRENTLY USED
        % =================================================================
            f=open ('DisplayQSens.fig');
            h=guihandles(f);
            set(h.txtExp2,'String',num2str(QSens.ppExponent,4));
            set(h.txtExp4,'String',num2str(QSens.nsExponent,4));
            set(h.txtExp6,'String',num2str(QSens.nsExponent,4));
            set(h.txtPer2,'String',num2str(QSens.qPPoptperdiff,2));
            set(h.txtPer3,'String',num2str(QSens.qCNSperdiff,2));
            set(h.txtPer4,'String',num2str(QSens.qCNSoptperdiff,2));
            set(h.txtPer5,'String',num2str(QSens.q3pNSperdiff,2));
            set(h.txtPer6,'String',num2str(QSens.q3pNSoptperdiff,2));

            drawnow;
        end
        % =================================================================
        function DisplayTable(QSens,h)
        %
        % Displays discharge sensitivity in predefined table with 
        % h as handle.
        % =================================================================
            qTable(1,1)={'Power'};
            qTable(1,2)={'Power'};
            qTable(1,3)={'0.1667'};
            qTable(1,4)={'Reference'};
            qTable(2,1)={'Power'};
            qTable(2,2)={'Power'};
            qTable(2,3)={num2str(QSens.ppExponent,'%6.4f')};
            qTable(2,4)={num2str(QSens.qPPoptperdiff,'%6.2f')};
            qTable(3,1)={'Constant'};
            qTable(3,2)={'No Slip'};
            qTable(3,3)={'0.1667'};
            qTable(3,4)={num2str(QSens.qCNSperdiff,'%6.2f')};                
            qTable(4,1)={'Constant'};
            qTable(4,2)={'No Slip'};
            qTable(4,3)={num2str(QSens.nsExponent,'%6.4f')};
            qTable(4,4)={num2str(QSens.qCNSoptperdiff,'%6.2f')};                
            qTable(5,1)={'3-Point'};
            qTable(5,2)={'No Slip'};
            qTable(5,3)={'0.1667'};
            qTable(5,4)={num2str(QSens.q3pNSperdiff,'%6.2f')}; 
            qTable(6,1)={'3-Point'};
            qTable(6,2)={'No Slip'};
            qTable(6,3)={num2str(QSens.nsExponent,'%6.4f')};
            qTable(6,4)={num2str(QSens.q3pNSoptperdiff,'%6.2f')};
            if ~isnan(QSens.manTop)
                qTable{7,1}=QSens.manTop;
                qTable{7,2}=QSens.manBot;
                qTable(7,3)={num2str(QSens.manExp,'%6.4f')};
                qTable(7,4)={num2str(QSens.qManperdiff,'%6.2f')};    
            end
            set(h,'Data',nan);
            set(h,'Data',qTable);
            
        end
            
    end
end
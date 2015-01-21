%possible options are
%            case 'fontsize'
%                fs = varargin{nn+1};
%            case 'xlabelfontsize'
%                xlabelfs = varargin{nn+1};
%            case 'ylabelfontsize'
%                ylabelfs = varargin{nn+1};
%            case 'zlabelfontsize'
%                zlabelfs = varargin{nn+1};
%            case 'titlefontsize'
 %               titlefs = varargin{nn+1};
%            case 'legendfontsize'
%                lfs = varargin{nn+1};
%            case 'markersize'
%                ms = varargin{nn+1};
%            case 'linewidth'
%                lw = varargin{nn+1};
%            case 'legendbox'
%                lb = varargin{nn+1};
%            case 'plotbox'
%                pb = varargin{nn+1};
%where the box takes 'on' or 'off'
%defaults to 
%all fs sizes 16, otherwise, fontsize is for plot boxes, all else specified
%individually
%boxes 'on'
%linewidth 3
%makersize 12


function makePresentationPerfect(figHandle,varargin)

%set defaults
global fs
global legendfs %legend
global xlabelfs
global ylabelfs
global zlabelfs
global titlefs
global lw
global ms
global lb
global pb
fs = 16;
legendfs = fs;
xlabelfs = fs;
ylabelfs = fs;
zlabelfs = fs;
titlefs = fs;
lw = 3;
ms = 20; % was 12 before
lb = 'on';
pb = 'on';

fH = figHandle;
if nargin > 1
    for nn = 1:2:length(varargin)
        lower(varargin{nn})
        switch lower(varargin{nn})
            case 'fontsize'
                fs = varargin{nn+1};
                legendfs = fs;
                xlabelfs = fs;
                ylabelfs = fs;
                zlabelfs = fs;
                titlefs = fs;
            case 'xlabelfontsize'
                xlabelfs = varargin{nn+1};
            case 'ylabelfontsize'
                ylabelfs = varargin{nn+1};
            case 'zlabelfontsize'
                zlabelfs = varargin{nn+1};
            case 'titlefontsize'
                titlefs = varargin{nn+1};
            case 'legendfontsize'
                legendfs = varargin{nn+1};
            case 'markersize'
                ms = varargin{nn+1};
            case 'linewidth'
                lw = varargin{nn+1};
            case 'legendbox'
                lb = varargin{nn+1};
            case 'plotbox'
                pb = varargin{nn+1};
            otherwise
                error(['unknown option ' varargin{nn}])
        end
    end
end
child = get(figHandle,'Children');

if ~isempty(child)
    for cc = 1:length(child)
        examineChild(child(cc))
    end
end



function examineChild(child)
global fs
global xlabelfs
global ylabelfs
global zlabelfs
global titlefs
global legendfs
global lw
global ms
global lb
global pb

    type = get(child,'Type')
    switch (type)
        case 'hggroup'
            set(child,'FontSize',fs);
        case 'axes'
            set(child,'FontSize',fs);
            if strcmp(lower(get(child,'Tag')),'legend')
                set(child,'Box',lb);
                set(child,'FontSize',legendfs);
            else
                set(child,'Box',pb);
            end
            xh = get(child,'XLabel');%labels are not children of axes, but have the handles in these fields
            yh = get(child,'YLabel');%labels are not children of axes, but have the handles in these fields
            zh = get(child,'ZLabel');%labels are not children of axes, but have the handles in these fields
            th = get(child,'Title');%labels are not children of axes, but have the handles in these fields
            set(xh,'FontSize',xlabelfs);
            set(yh,'FontSize',ylabelfs);
            set(zh,'FontSize',zlabelfs); 
            set(th,'FontSize',titlefs);
        case 'line'
            marker = get(child,'Marker');
            if strcmp(lower(marker),'none')
                set(child,'LineWidth',lw);
            else
                set(child,'MarkerSize',ms);
            end
        case 'text'
            parentType = get(get(child,'parent'),'Tag');
            if ~strcmp(lower(parentType),'legend')
                set(child,'FontSize',fs);;
            end
        otherwise
    end
    child2 = get(child,'Children');
    if ~isempty(child2)
        for cc = 1:length(child2)
            examineChild(child2(cc));
        end
    end

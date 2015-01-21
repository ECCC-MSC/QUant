%% SPD_Steph_dub

% Plots probability density across spectrum. Input array "s" has vector of
% frequency values as first column, and first row is ignored as another
% heading vector, i.e. data is s(2:r,2:c), where [r,c] = size(s).

[r,c] = size(s);        %size of input array
fvec = s(2:r,1);        %vector of frequencies

sv = s(2:r,2:c);        %data array

hind = 0.1;             %histogram (PDF) bin width
dbvec = 30:hind:160;    %dB values at which to calculate empirical PDF
ldbv = length(dbvec);   %number of bins in PDF

leq = 10*log10(mean(10.^(sv/10),2));            %calculate linear mean
p = prctile(sv,0:100,2);                        %array of percentiles
d = hist(sv.',dbvec)/(hind*(c-1));              %SPD array

d(d == 0) = NaN;                %suppress plotting of empty hist bins

[X,Y] = meshgrid(fvec,dbvec);   %axis arrays for SPD pcolor plot

%% Plot

figure(666)
set(figure(666),'color','w')
hold off
    
g = pcolor(X,Y,d);                          %SPD
set(g,'LineStyle','none')
colorbar

hold on
semilogx(fvec,leq,'m','linewidth',2)        %linear mean
semilogx(fvec,p(:,99),'k','linewidth',2)    %percentiles
semilogx(fvec,p(:,95),'k','linewidth',2)
semilogx(fvec,p(:,50),'k','linewidth',2)
semilogx(fvec,p(:,5),'k','linewidth',2)
semilogx(fvec,p(:,1),'k','linewidth',2)

caxis([0 0.05])
xlabel('Frequency [ Hz ]','fontsize',14,'fontname','Arial')
ylabel('PSD [ dB re 1 \muPa^2 Hz^-^1 ]','fontsize',14,'fontname','Arial')
set(gca,'XScale','log','TickDir','out','layer','top','fontsize',14,'fontname','Arial')
set(colorbar,'fontsize',14,'fontname','Arial')
ylabel(colorbar,'Empirical Probability Density','fontsize',14,'fontname','Arial')
ylim([40 140])
xlim([min(fvec) max(fvec)])

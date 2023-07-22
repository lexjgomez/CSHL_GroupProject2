function [] = movshonize(fs, wholefigure)
%Pascal Wallisch
%02/25/2011
%This function does the Movshon procedure on a figure. 
%Usage: movshonize(fontsize, applytowholefigure) where "apply" is 0 or 1
box off
set(gca,'TickDir','out')
set(gca,'XScale', 'linear')
set(gca,'fontsize',fs)
set(gca,'fontangle','oblique')
set(gca,'ticklength',[0.01 0.025])
set(gca,'layer','top')
axis square

if wholefigure == 1 
%Applies to whole figure
temp = findall(gcf,'Type','text');
set(temp,'FontSize',fs);
set(temp,'FontAngle','italic');
set(temp,'FontName','Helvetica');
%set(temp,'Interpreter','Latex');

end

end




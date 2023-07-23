function PSTH_GUI(spike_PSTH, trialsInfo)


global GUI PlottingData
PlottingData = {spike_PSTH, trialsInfo};
GUI = modalityPSTH_GUI();
set(GUI.neuronmenu, 'string', 1:size(spike_PSTH, 1));


end








function GUIHandles = modalityPSTH_GUI()

GUIHandles.Framework = figure('Position', [300 300 600 500], 'name', 'All Plots', 'numbertitle', 'off', 'MenuBar', 'none', 'Resize', 'off');
GUIHandles.neuronmenu = uicontrol('Style', 'popupmenu', 'Position', [60 430 150 25], 'Callback', {@update_neuronmenu_callback});

GUIHandles.Plotting1 = axes('Units', 'pixels', 'Position', [100 70 450 300],'tickdir','out');

end



%% Plotting function
function update_neuronmenu_callback(Obj, eventdata, handles)

global GUI PlottingData

windowsize = 0.1;
list = str2num(get(Obj, 'String'));
i = get(Obj, 'Value');
idx = list(i);

[spike_PSTH, trialsInfo] = PlottingData{1, 1:2};
thisNeuron = squeeze(spike_PSTH(i, :, :));

idx_L = find(trialsInfo.contrast < 0);
idx_None = find(trialsInfo.contrast == 0);
idx_R = find(trialsInfo.contrast > 0);

axes(GUI.Plotting1); 
cla
hold on

a = StandardErrorShade((thisNeuron(idx_L, 1:100)./windowsize),0.4,[0 0.4470 0.7410]);
b = StandardErrorShade((thisNeuron(idx_R, 1:100)./windowsize),0.4,[0.9290 0.6940 0.1250]);
c = StandardErrorShade((thisNeuron(idx_None, 1:100)./windowsize),0.4,[0.4660 0.6740 0.1880]);

yl = ylim;
line([50 50],yl,'Color','black','LineStyle','--');
% line([300 300],yl,'Color','black','LineStyle','--');

xlabel('Time from Stimulus On (s)');
xlim([0 100])
ylabel('Firing Rate (sp/s)');
xticks([0 50 100]);
xticklabels({'-0.5','0','0.5'});
legend([a, b, c], {'Left', 'Right', 'No Stimulus'});
set(gca,'box','off');
set(gca,'tickdir','out');



end
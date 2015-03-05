function changeTrialsPerScreen(newValue,currentLFs,start_trial)

global trialsPerScreen

trialsPerScreen=newValue;

makePlots(currentLFs,trialsPerScreen,start_trial,[],1)
function getDirectCaseDataList_display(out,num)
%Displays output of getDirectCaseDataList script
disp(['Main Cell ID: ' num2str(out{num}{1})])
disp(['Contacting Cell IDs: ' num2str(out{num}{2})])
for i = 1:length(out{num}{3});
	disp([num2str(out{num}{2}(i)) ' --> ' num2str(out{num}{3}{i})])
end

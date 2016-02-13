function getDirectCaseDataList_display(out,num)
%Displays output of getDirectCaseDataList script
disp(' ')
disp(['Main Cell ID: ' num2str(out{num}{1}) '. Main electrode: ' num2str(out{num}{2})])
disp(['Contacting Cell IDs: ' num2str(out{num}{3})])
for i = 1:length(out{num}{5});
	disp([num2str(out{num}{3}(i)) ' --> ' num2str(out{num}{5}{i}) '. Main electrode: ' num2str(out{num}{4}(i))])
end
disp(' ')

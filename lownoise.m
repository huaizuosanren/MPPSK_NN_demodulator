function datadef= lownoise (data)
% datadef is the predict initial data
% data is the initial data that have been interfeared by noise

datadef=data;
for di=3:length (data)
	if abs(data(di))>=2
		
		datadef(di)=(data(di-1)+data(di-2))/2;
	end
end


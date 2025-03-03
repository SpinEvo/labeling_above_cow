%***********************************************************************
%Method to choose a hadamard matrix to pick columns from to use as an ideal
%vessel encoding scheme
%***********************************************************************

function master_hadamard = fun_HW_ChooseHadamard(vesloc)

%choose which hadamard matrix to use as a base according to the number of
%vessels
if size(vesloc,2) <= 3 && size(vesloc,2) > 0
    master_hadamard = hadamard(4);
    
elseif size(vesloc,2) > 3 && size(vesloc,2) <= 7
    master_hadamard = hadamard(8);

elseif size(vesloc,2) > 7 && size(vesloc,2) <= 11
    master_hadamard = hadamard(12);
    
elseif size(vesloc,2) > 11 && size(vesloc,2) <= 15
    master_hadamard = hadamard(16);
    
else
    error('The number of vessels cannot exceed 16 or be less than 1')
end

random_order = randperm(size(master_hadamard,2));
master_hadamard = master_hadamard(:,random_order);

master_hadamard = master_hadamard(:,2:size(vesloc,2)+1);

return

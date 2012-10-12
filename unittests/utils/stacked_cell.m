% From the stack matrix Ja, build the cell hierarchized matrix J{}.
function Au = stacked_cell( A );

p=length(A);
Au = [];
for i=1:p
    Au = [Au;A{i}];
end

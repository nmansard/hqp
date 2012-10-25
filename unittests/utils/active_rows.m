function [A,b] = active_rows(h,data);

% active_rows returns the stack of the input active quantity.
%
%% Synopsis:
%    A     = active_rows(h)
%    [A,b] = active_rows(h)
%    A     = active_rows(h,data)
%    [A,b] = active_rows(h,data)
%
%% Input:
%    h      is a cell containing an "active" field that is an array of indices.
%    data   is a cell whose length is the same as h. Each cell k contains a
%              matrix (or colvec) of h(k).mmax rows. 
%
%% Output: 
%    A      the function selects the active rows of the cell "data".
%              (corresponding to h(k).active) and stack them to return.  If no 
%              data is input, then data is chosen from h.A. 
%    b      the stacked h.b
%
%

% --------------------------------------------------------------------
if nargin==1
    data=stacked_cell({h.A});
end
    
if strcmp(class(data),'cell')
    p=min(length(h),length(data));
else
    p=length(h);
end        
% --------------------------------------------------------------------

b=[];
A=[];
m=0;
for k=1:p
    ia = h(k).active;
    if strcmp(class(data),'cell')
        A=[A; data{k}(ia)];
    else
        A=[A;  data(m+ia,:)];
    end        
    b = [b ; h(k).b( h(k).activeb )];
    m=m+h(k).mmax;
end

function [ aset bounds ]= initset(btype,ainit, boundinit);

% initset: From the type of the type of the constraints, and maybe an initial
% set, activate all the TWIN (equality constraints).
%% Synopsis:
%     [ aset bounds ] = initset(btype,ainit, boundinit)
% 
%% Input:
%    btype       type of the constraint, in =,<=,>=,<=<=
%    ainit       initial guess for the active set.
%    boundinit   bounds initialy active.
%% Output:
%    aset        initial active set plus the equality constraint (if not initial
%                  guess, returns only the list of the equality constraints).
%    bounds      corresponding active bounds
%
%

% --- DEFAULT ARGUMENTS --------------------------------------------------------
if nargin==1
    for k=1:length(btype)
        ainit{k} = ones(0,1); boundinit{k}=ones(0,1);
    end
end
% ------------------------------------------------------------------------------

constants;

aset   = {};
bounds = {};
for k=1:length(btype)
    % Octave needs to now wether it is a row or a column vector.
    if size(ainit{k},2)==1      ainit{k}     = ainit{k}';     end
    if size(boundinit{k},2)==1  boundinit{k} = boundinit{k}'; end
        
    atwin     = setdiff(  find(btype{k}==Etwin),ainit{k}  );
    l         = length(atwin);

    aset{k}                = ainit{k};
    aset{k}(end+1:end+l)   = atwin;
    bounds{k}              = boundinit{k};
    bounds{k}(end+1:end+l) = Etwin;
end

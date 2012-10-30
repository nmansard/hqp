function rand_init(seed)

if exist('OCTAVE_VERSION')
    
    rand('state',seed);
    randn('state',rand*1000);
    rande('state',rand*1000);
    
else
   
    rng(seed,'twister');

end
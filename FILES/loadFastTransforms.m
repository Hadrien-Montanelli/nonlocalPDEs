function loadFastTransforms()

% First time through, load the Julia functions.
persistent loaded;

if isempty(loaded)
  jl_file = fullfile(fileparts(mfilename('fullpath')), 'fastTransforms.jl');
  jl.include(jl_file);
  jl_file = fullfile(fileparts(mfilename('fullpath')), 'evaluateLambda.jl');
  jl.include(jl_file);
  loaded = true;
end

end
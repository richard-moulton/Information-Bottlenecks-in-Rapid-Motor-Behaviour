%% A function replacing MATLAB's dir
%  Returns a listing of files in a directory, not including '.' and '..'
%  (these are returned by the library version of the function dir).
%
%  3 December, 2021

function listing = dir2(varargin)

if nargin == 0
    name = '.';
elseif nargin == 1
    name = varargin{1};
else
    error('Too many input arguments.')
end

listing = dir(name);

indices = [];
n    = 0;       % Number of matches against . and ..
k    = 1;       % Index of file being considered

while n < 2 && k <= length(listing)
    if any(strcmp(listing(k).name, {'.', '..'}))
        indices(end + 1) = k; %#ok<AGROW>
        n = n + 1;
    end
    k = k + 1;
end

listing(indices) = [];
end
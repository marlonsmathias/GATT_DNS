function y = diff(x,n,dim)
	if nargin < 3
		dim = find(size(x)>1,1,'first');
    if isempty(dim)
      y = [];
      return
    end
	end
	if nargin < 2 || isempty(n)
		n = 1;
	end
	if n >= size(x,dim)
		y = [];
		return
	end
	y = builtin('diff',x,n,dim);
end
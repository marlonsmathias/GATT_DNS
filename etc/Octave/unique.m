function [C,ia,ic] = unique(A,op1,op2)
  if nargin < 2
    op1 = [];
  end
  if nargin < 3
    op2 = [];
  end
  
  if strcmp(op1,'sorted')
    op1 = 'first';
  end
  if strcmp(op2,'sorted')
    op2 = 'first';
  end
  
  if ~strcmp(op1,'stable') && ~strcmp(op2,'stable')
    if isempty(op1)
      [C,ia,ic] = unique_builtin(A);
    elseif isempty(op2)
      [C,ia,ic] = unique_builtin(A,op1);
    else
      [C,ia,ic] = unique_builtin(A,op1,op2);
    end
    return
  end
  
  if strcmp(op1,'rows')
    [~,ia,ic] = unique_builtin(A,'rows','first');
    [ia,ind] = sort(ia);
    icnew = 0*ic;
    for i = 1:length(ia)
      icnew(ic==ind(i)) = i;
    end
    ic = icnew;
    C = A(ia,:);
  else
    [~,ia,ic] = unique_builtin(A,'first');
    [ia,ind] = sort(ia);
    icnew = 0*ic;
    for i = 1:length(ia)
      icnew(ic==ind(i)) = i;
    end
    ic = icnew;
    C = A(ia);
  end
  
  
end

function [y, i, j] = unique_builtin (x, varargin)

  if (nargin < 1)
    print_usage ();
  elseif (! (isnumeric (x) || islogical (x) || ischar (x) || iscellstr (x)))
    error ("unique: X must be an array or cell array of strings");
  endif

  if (nargin > 1)
    ## parse options
    if (! iscellstr (varargin))
      error ("unique: options must be strings");
    endif

    optrows  = any (strcmp ("rows", varargin));
    optfirst = any (strcmp ("first", varargin));
    optlast  = any (strcmp ("last", varargin));
    if (optfirst && optlast)
      error ('unique: cannot specify both "first" and "last"');
    elseif (optfirst + optlast + optrows != nargin-1)
      error ("unique: invalid option");
    endif

    if (optrows && iscellstr (x))
      warning ('unique: "rows" is ignored for cell arrays');
      optrows = false;
    endif
  else
    optrows = false;
    optfirst = false;
  endif

  if (issparse (x) && ! optrows && nargout <= 1)
    if (nnz (x) < numel (x))
      y = unique ([0; nonzeros(x)], varargin{:});
    else
      ## Corner case where sparse matrix is actually full
      y = unique (full (x), varargin{:});
    endif
    return;
  endif

  if (optrows)
    n = rows (x);
    dim = 1;
  else
    n = numel (x);
    dim = (rows (x) == 1) + 1;
  endif

  y = x;
  ## Special cases 0 and 1
  if (n == 0)
    if (! optrows && isempty (x) && any (size (x)))
      if (iscellstr (y))
        y = cell (0, 1);
      else
        y = zeros (0, 1, class (y));
      endif
    endif
    i = j = [];
    return;
  elseif (n == 1)
    i = j = 1;
    return;
  endif

  if (optrows)
    if (nargout > 1)
      [y, i] = sortrows (y);
    else
      y = sortrows (y);
    endif
    match = all (y(1:n-1,:) == y(2:n,:), 2);
    y(match,:) = [];
  else
    if (! isvector (y))
      y = y(:);
    endif
    if (nargout > 1)
      [y, i] = sort (y);
    else
      y = sort (y);
    endif
    if (iscellstr (y))
      match = strcmp (y(1:n-1), y(2:n));
    else
      match = (y(1:n-1) == y(2:n));
    endif
    y(match) = [];
  endif

  if (isargout (3))
    j = i;
    if (dim == 1)
      j(i) = cumsum ([1; ! match]);
    else
      j(i) = cumsum ([1, ! match]);
    endif
  endif

  if (isargout (2))
    idx = find (match);
    if (optfirst)
      idx += 1;   # in-place is faster than other forms of increment
    endif
    i(idx) = [];
  endif

endfunction
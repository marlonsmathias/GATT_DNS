function save(varargin)

args = cell(1,nargin);
args{1} = varargin{1};

for iii = 2:nargin
    if ~strcmp(varargin{iii}(1),'-')
        eval([varargin{iii} ' = evalin(''caller'',''' varargin{iii} ''');']);
    end
    
    if ~strcmp(varargin{iii},'-v7.3')
        args{iii} = varargin{iii};
    else
        args{iii} = '-v7';
    end
    
end

builtin('save',args{:});

end
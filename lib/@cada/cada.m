classdef cada
  % cada classdef file for use with the ADiGator package.
  %
  % The cada class is used internally by the ADiGator algorithm in order to
  % overload mathematical and organizational operations. It is not intended
  % to be used outside of the ADiGator algorithm.
  %
  % Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
  % Distributed under the GNU General Public License version 3.0
  properties
    % Unique integer identifier
    id
    % Function information (size, name of function variable, etc.)
    func
    % Derivative information (locations, name of derivative variable)
    deriv
  end
  
  % Constructor
  methods
    function x = cada(varID,func,deriv)
      % cada constructor function for the ADiGator algorithm. 
      %
      % Creates the overloaded object given the .func and .deriv
      % properties. If not all .func properties are defined, it sets them
      % to empty.
      if nargin == 1 && isstruct(varID)
        y = varID;
        x.id    = y.id;
        x.func  = y.func;
        x.deriv = y.deriv;
      else
        if ~isfield(func,'name')
          func.name = [];
        end
        if ~isfield(func,'size')
          func.size = [0 0];
        end
        if ~isfield(func,'zerolocs')
          func.zerolocs = [];
        end
        if ~isfield(func,'value')
          func.value = [];
        end
        x.id    = varID;
        x.func  = func;
        x.deriv = deriv;
      end
    end
  end
  
  % Utility ordinary methods called from ADiGator routines
  methods
    flag = cadaCheckForDerivs(x)
    NewVar = cadaPrintReMap(x,OverVar,VarCount)
    xOut = cadaOverMap(x,DAid)
    z    = cadaUnionVars(x,y)
    adigatorAnalyzeForData(FORCOUNT,dummyVar)
    adigatorPrintOutputIndices(x)
    x = adigatorStructAnalyzer(x,xStr,subsflag)
    varargout = adigatorVarAnalyzer(FunString,varargin)
  end
  
  
  methods 
    % Generalized derivative procedures
    [y,varargout] = cadaEmptyEval(varargin)
    y = cadaunarymath(x,zeroflag,callerstr)
    y = cadaunarylogical(x,callerstr,dim)
    y = cadacreatearray(callerstr,varargin)
    z = cadabinaryarraymath(x,y,xzeroflag,yzeroflag,callerstr)
  end 
  
  % Overloaded array creation operations
  methods
    function y = cell(varargin)
      % CADA overloaded version of function CELL.
      y = cadacreatearray('cell',varargin{:});
    end
    function y = eye(varargin)
      % CADA overloaded version of function EYE
      y = cadacreatearray('eye',varargin{:});
    end
    function y = false(varargin)
      % CADA overloaded version of function FALSE
      y = cadacreatearray('false',varargin{:});
    end
    function y = nan(varargin)
      % CADA overloaded version of function NAN.
      y = cadacreatearray('nan',varargin{:});
    end
    function y = ones(varargin)
      % CADA overloaded version of function ONES
      y = cadacreatearray('ones',varargin{:});
    end
    function y = true(varargin)
      % CADA overloaded version of function TRUE
      y = cadacreatearray('true',varargin{:});
    end
    function y = zeros(varargin)
      % CADA overloaded version of function ZEROS.
      y = cadacreatearray('zeros',varargin{:});
    end
  end
  
  % Overloaded unary math array operations
  methods
    
    
    function y = abs(x)
      % CADA overloaded ABS function
      global ADIGATOR
      if ADIGATOR.OPTIONS.COMPLEX
        y = sqrt(real(x).^2 + imag(x).^2);
      else
        y = cadaunarymath(x,1,'abs');
      end
    end
    function y = acos(x)
      % CADA overloaded ACOS function
      y = cadaunarymath(x,0,'acos');
    end
    function y = acosd(x)
      % CADA overloaded ACOSD function
      y = cadaunarymath(x,0,'acosd');
    end
    function y = acosh(x)
      % CADA overloaded ACOSH function
      y = cadaunarymath(x,0,'acosh');
    end
    function y = acot(x)
      % CADA overloaded ACOT function
      y = cadaunarymath(x,0,'acot');
    end
    function y = acotd(x)
      % CADA overloaded ACOTD function
      y = cadaunarymath(x,0,'acotd');
    end
    function y = acoth(x)
      % CADA overloaded ACOTH function
      y = cadaunarymath(x,0,'acoth');
    end
    function y = acsc(x)
      % CADA overloaded ACSC function
      y = cadaunarymath(x,0,'acsc');
    end
    function y = acscd(x)
      % CADA overloaded ACSCD function
      y = cadaunarymath(x,0,'acscd');
    end
    function y = angle(x)
      % CADA overloaded version of ANGLE
      y = atan2(imag(x),real(x));
    end
    function y = acsch(x)
      % CADA overloaded ACSCH function
      y = cadaunarymath(x,0,'acsch');
    end
    function y = asec(x)
      % CADA overloaded ASEC function
      y = cadaunarymath(x,0,'asec');
    end
    function y = asecd(x)
      % CADA overloaded ASECD function
      y = cadaunarymath(x,0,'asecd');
    end
    function y = asech(x)
      % CADA overloaded ASECH function
      y = cadaunarymath(x,0,'asech');
    end
    function y = asin(x)
      % CADA overloaded ASIN function
      y = cadaunarymath(x,1,'asin');
    end
    function y = asind(x)
      % CADA overloaded ASIND function
      y = cadaunarymath(x,1,'asind');
    end
    function y = asinh(x)
      % CADA overloaded ASINH function
      y = cadaunarymath(x,1,'asinh');
    end
    function y = atan(x)
      % CADA overloaded ATAN function
      y = cadaunarymath(x,1,'atan');
    end
    function y = atand(x)
      % CADA overloaded ATAND function
      y = cadaunarymath(x,1,'atand');
    end
    function y = atanh(x)
      % CADA overloaded ATANH function
      y = cadaunarymath(x,1,'atanh');
    end
    function y = ceil(x)
      % CADA overloaded CEIL function - derivative set to
      % zero
      global ADIGATOR
      NUMvod  = ADIGATOR.NVAROFDIFF;
      % Remove derivatives of x
      x.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
      y = cadaunarymath(x,1,'ceil');
      if ~ADIGATOR.RUNFLAG
        numvars = find(ADIGATOR.VARINFO.NAMELOCS(:,1),1,'last');
        cadaCancelDerivs(y.id,numvars);
      end
    end
    function y = conj(x)
    % CADA overloaded CONJ function
      y = cadaunarymath(x,1,'conj');
    end
    function y = cos(x)
      % CADA overloaded COS function
      y = cadaunarymath(x,0,'cos');
    end
    function y = cosd(x)
      % CADA overloaded COSD function
      y = cadaunarymath(x,0,'cosd');
    end
    function y = cosh(x)
      % CADA overloaded COSH function
      y = cadaunarymath(x,0,'cosh');
    end
    function y = cot(x)
      % CADA overloaded COT function
      y = cadaunarymath(x,0,'cot');
    end
    function y = cotd(x)
      % CADA overloaded COTD function
      y = cadaunarymath(x,0,'cotd');
    end
    function y = coth(x)
      % CADA overloaded COTH function
      y = cadaunarymath(x,0,'coth');
    end
    function y = csc(x)
      % CADA overloaded CSC function
      y = cadaunarymath(x,0,'csc');
    end
    function y = cscd(x)
      % CADA overloaded CSCD function
      y = cadaunarymath(x,0,'cscd');
    end
    function y = csch(x)
      % CADA overloaded CSCH function
      y = cadaunarymath(x,0,'csch');
    end
    function y = erf(x)
      % CADA overloaded ERF function
      y = cadaunarymath(x,1,'erf');
    end
    function y = exp(x)
      % CADA overloaded EXP function
      y = cadaunarymath(x,0,'exp');
    end
    function y = fix(x)
      % CADA overloaded FIX function - derivative set to
      % zero
      global ADIGATOR
      NUMvod  = ADIGATOR.NVAROFDIFF;
      % Remove derivatives of x
      x.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
      y = cadaunarymath(x,1,'fix');
      if ~ADIGATOR.RUNFLAG
        numvars = find(ADIGATOR.VARINFO.NAMELOCS(:,1),1,'last');
        cadaCancelDerivs(y.id,numvars);
      end
    end
    function y = floor(x)
      % CADA overloaded FLOOR function - derivative set to
      % zero
      global ADIGATOR
      NUMvod  = ADIGATOR.NVAROFDIFF;
      % Remove derivatives of x
      x.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
      y = cadaunarymath(x,1,'floor');
      if ~ADIGATOR.RUNFLAG
        numvars = find(ADIGATOR.VARINFO.NAMELOCS(:,1),1,'last');
        cadaCancelDerivs(y.id,numvars);
      end
    end
    function y = full(x)
      % CADA overloaded FULL function.
      y = cadaunarymath(x,1,'full');
    end
    function y = imag(x)
    % CADA overloaded IMAG function
      y = cadaunarymath(x,1,'imag');
    end
    function y = log(x)
      % CADA overloaded LOG function
      y = cadaunarymath(x,0,'log');
    end
    function y = log10(x)
      % CADA overloaded LOG10 function
      y = cadaunarymath(x,0,'log10');
    end
    function y = real(x)
      % CADA overloaded REAL function
      y = cadaunarymath(x,1,'real');
    end
    function y = round(x)
      % CADA overloaded ROUND function - derivative set to
      % zero
      global ADIGATOR
      NUMvod  = ADIGATOR.NVAROFDIFF;
      % Remove derivatives of x
      x.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
      y = cadaunarymath(x,1,'round');
      if ~ADIGATOR.RUNFLAG
        numvars = find(ADIGATOR.VARINFO.NAMELOCS(:,1),1,'last');
        cadaCancelDerivs(y.id,numvars);
      end
    end
    function y = sec(x)
      % CADA overloaded SEC function
      y = cadaunarymath(x,0,'sec');
    end
    function y = secd(x)
      % CADA overloaded SECD function
      y = cadaunarymath(x,0,'secd');
    end
    function y = sech(x)
      % CADA overloaded SECH function
      y = cadaunarymath(x,0,'sech');
    end
    function y = sign(x)
      % CADA overloaded SIGN function
      for Vcount = 1:length(x.deriv)
        if ~isempty(x.deriv(Vcount).nzlocs)
          warning('derivative of discontinuous sign function - making derivatives zero'); %#ok<WNTAG>
          x.deriv(Vcount).nzlocs = [];
          x.deriv(Vcount).name = [];
        end
      end
      y = cadaunarymath(x,1,'sign');
    end
    function y = sin(x)
      % CADA overloaded SIN function
      y = cadaunarymath(x,1,'sin');
    end
    function y = sind(x)
      % CADA overloaded SIND function
      y = cadaunarymath(x,1,'sind');
    end
    function y = sinh(x)
      % CADA overloaded SINH function
      y = cadaunarymath(x,1,'sinh');
    end
    function y = sqrt(x)
      % CADA overloaded SQRT function
      y = cadaunarymath(x,1,'sqrt');
    end
    function y = tan(x)
      % CADA overloaded TAN function
      y = cadaunarymath(x,1,'tan');
    end
    function y = tand(x)
      % CADA overloaded TAND function
      y = cadaunarymath(x,1,'tand');
    end
    function y = tanh(x)
      % CADA overloaded TANH function
      y = cadaunarymath(x,1,'tanh');
    end
    function y = uminus(x)
      % CADA overloaded UMINUS function
      y = cadaunarymath(x,1,'uminus');
    end
    function y = uplus(x)
      % CADA overloaded UPLUS function
      y = cadaunarymath(x,1,'uplus');
    end
  end
  
  % Overloaded unary logical operations
  methods
    function y = all(x,varargin)
      % CADA overloaded version of function ALL
      if nargin == 2
        dim = varargin{1};
      elseif x.func.size(1) == 1
        dim = 2;
      else
        dim = 1;
      end
      y = cadaunarylogical(x,'all',dim);
      
    end
    function y = any(x,varargin)
      % CADA overloaded version of function ANY
      if nargin == 2
        dim = varargin{1};
      elseif x.func.size(1) == 1
        dim = 2;
      else
        dim = 1;
      end
      y = cadaunarylogical(x,'any',dim);
    end
    function y = logical(x)
      % CADA overloaded version of function LOGICAL
      y = cadaunarylogical(x,'logical',0);
    end
    function y = not(x)
      % CADA overloaded version of function NOT
      y = cadaunarylogical(x,'not',0);
    end
  end
  
  % Overloaded binary math array operations
  methods
    function z = atan2(x,y)
      % CADA overloaded ATAN2 function
      z = cadabinaryarraymath(x,y,1,1,'atan2');
    end
    function z = ldivide(x,y)
      % CADA overloaded LDIVIDE function
      z = cadabinaryarraymath(y,x,0,1,'rdivide');
    end
    function z = mod(x,y)
      % CADA overloaded MOD function
      z = cadabinaryarraymath(x,y,1,0,'mod');
    end
    function z = plus(x,y)
      % CADA overloaded PLUS function
      z = cadabinaryarraymath(x,y,0,0,'plus');
    end
    function z = power(x,y)
      % CADA overloaded POWER function
      z = cadabinaryarraymath(x,y,1,1,'power');
    end
    function z = rdivide(x,y)
      % CADA overloaded RDIVIDE function
      z = cadabinaryarraymath(x,y,1,0,'rdivide');
    end
    function z = rem(x,y)
      % CADA overloaded REM function
      z = cadabinaryarraymath(x,y,1,0,'rem');
    end
    function z = times(x,y)
      % CADA overloaded TIMES function
      z = cadabinaryarraymath(x,y,1,1,'times');
    end
    function z = minus(x,y)
      % CADA overloaded MINUS function
      z = cadabinaryarraymath(x,y,0,0,'minus');
    end
    [y,varargout] = max(x,varargin)
    [y,varargout] = min(x,varargin)
  end
  
  % Overloaded binary logical array operations
  methods
    function z = and(x,y)
      % CADA overloaded version of function AND
      z = cadabinarylogical(x,y,'and');
    end
    function z = eq(x,y)
      % CADA overloaded version of function EQ
      z = cadabinarylogical(x,y,'eq');
    end
    function z = ge(x,y)
      % CADA overloaded version of function GE
      z = cadabinarylogical(x,y,'ge');
    end
    function z = gt(x,y)
      % CADA overloaded version of function GT
      z = cadabinarylogical(x,y,'gt');
    end
    function z = le(x,y)
      % CADA overloaded version of function LE
      z = cadabinarylogical(x,y,'le');
    end
    function z = lt(x,y)
      % CADA overloaded version of function LT
      z = cadabinarylogical(x,y,'lt');
    end
    function z = ne(x,y)
      % CADA overloaded version of function NE
      z = cadabinarylogical(x,y,'ne');
    end
    function z = or(x,y)
      % CADA overloaded version of function OR
      z = cadabinarylogical(x,y,'or');
    end
    function z = xor(x,y)
      % CADA overloaded version of function XOR
      z = cadabinarylogical(x,y,'xor');
    end
  end
  
  % Overloaded matrix operations
  methods
    z = mtimes(x,y)
    y = sum(x,varargin)
    z = mldivide(x,y)
    z = mrdivide(x,y)
    y = inv(x)
  end
  
  % Overloaded sizing opertions
  methods
    function out = end(x,dim,ndim)
      % CADA overloaded END
      if ndim == 1
        out = length(x);
      else
        out = size(x,dim);
      end
    end
    y = colon(varargin)
    y = nnz(x)
    y = numel(x)
    y = length(x)
    varargout = size(x,varargin)
  end
  
  % Overloaded organizational operations
  methods
    function y = fliplr(x)
      % CADA overloaded version of function FLIPLR
      global ADIGATOR
      % Size will cancel x's derivatives, do not let it..
      tmp = ADIGATOR.VARINFO.NAMELOCS;
      ind = size(x,2):-1:1;
      ADIGATOR.VARINFO.NAMELOCS = tmp;
      s.type = '()';
      s.subs = {':',ind};
      y = subsref(x,s);
    end
    function y = flipud(x)
      % CADA overloaded version of function FLIPUD
      global ADIGATOR
      % Size will cancel x's derivatives, do not let it..
      tmp = ADIGATOR.VARINFO.NAMELOCS;
      ind = size(x,1):-1:1;
      ADIGATOR.VARINFO.NAMELOCS = tmp;
      s.type = '()';
      s.subs = {ind,':'};
      y = subsref(x,s);
    end
    y = subsref(x,s)
    y = subsasgn(x,s,b)
    y = reshape(x,varargin)
    y = repmat(x,varargin)
    y = transpose(x)
    function y = ctranspose(x)
      % CADA overloaded CTRANSPOSE
      global ADIGATOR
      if ADIGATOR.OPTIONS.COMPLEX
        y = transpose(conj(x));
      else
        y = transpose(x);
      end 
    end
    y = horzcat(varargin)
    y = vertcat(varargin)
    y = diag(x,varargin)
    y = sparse(varargin)
    y = nonzeros(x)
  end
  
  % High level overloaded operations
  methods
    y = prod(x,varargin)
    yi = interp1(varargin)
    zi = interp2(x,y,z,xi,yi,varargin)
    yi = ppval(pp,xi)
    zi = adigatorEvalInterp2pp(pp,xi,yi)
    function z = dot(x,y,varargin)
      % CADA overloaded DOT - just calls z = sum(x.*y) - lazy.
      global ADIGATOR
      if ADIGATOR.OPTIONS.COMPLEX
        x = conj(x);
      end
      if nargin == 2
        z = sum(x.*y);
      else
        z = sum(x.*y,varargin{1});
      end
    end
    z = cross(x,y,varargin)
  end
  
  % Misc overloaded operations
  methods
    y = isempty(x)
    z = isequal(x,y,varargin)
    z = isequalwithequalnans(x,y,varargin)
    function s = num2str(x)
      % CADA overloaded version of function NUM2STR.
      s = num2str(x.func.value);
    end
    K = sub2ind(Asize,I,J)
  end
  
end
classdef adigatorInput
  % adigatorInput class - used only for specifying derivative and auxiliary
  % inputs to adigator(). Called by adigatorCreateDerivInput and
  % adigatorCreateAuxInput.
  %
  % Constructor syntax y = adigatorInput(func,deriv)
  %
  % See also: adigatorCreateDerivInput adigatorCreateAuxInput adigator
   
  properties
    % Function information - size and any fixed information
    func
    % Derivative information on size/name of variable of differentiation and non-zero locations
    deriv
  end
  
  methods
    function y = adigatorInput(func,deriv)
      % adigatorInput constructor function
      if ~isstruct(func)
        error('first func input to adigatorInput must be a structure with field .size');
      end
      if ~isfield(func,'size') || length(func.size) ~= 2
        error('adigatorInput must have defined function size of length 2');
      end
      if ~isfield(func,'value')
        func.value = [];
      end
      y.func = func;
      
      if isempty(deriv)
        y.deriv = [];
      elseif ~isstruct(deriv)
        error('deriv input must be structure');
      end
      
      if ~isempty(deriv)
        if ~isfield(deriv,'vodname') || ~isfield(deriv,'vodsize') || ~isfield(deriv,'nzlocs')
          error('deriv structure must be empty or contain fields: vodname, vodsize, nzlocs');
        end
      end
      
      ysize = func.size;
      for i = 1:length(deriv)
        xsize = deriv(i).vodsize;
        if length(xsize) ~= 2
          error('deriv.vodsize must be of length 2');
        end
        m = prod(ysize); n = prod(xsize);
        nzlocs = deriv(i).nzlocs;
        if size(nzlocs,2) ~= 2
          error(['nzlocs must be given in row/column format as an Nx2 ',...
            'integer matrix where first column corresponds to function ',...
            'values and second column corresponds to the variable ',...
            'of differentiation']);
        end
        if any(nzlocs(:,1) > m) || any(nzlocs(:,1) ~= abs(floor(nzlocs(:,1))))
          error('nzlocs either outside of bounds or not positive integers');
        elseif any(nzlocs(:,2) > n) || any(nzlocs(:,2) ~= abs(floor(nzlocs(:,2))))
          error('nzlocs either outside of bounds or not positive integers');
        end
        
        vodname = deriv(i).vodname;
        if ~ischar(vodname)
          error('vodname must be string with no spaces');
        end
      end
      
      y.deriv = deriv;
    end
    
  end
  
end
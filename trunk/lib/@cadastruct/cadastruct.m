classdef cadastruct
  % cadastruct classdef file for use with the ADiGator  package.
  %
  % Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
  % Distributed under the GNU General Public License version 3.0
  
  properties
    id 
    name 
    val
    arrayflag
  end
  
  % Constructor
  methods
    function y = cadastruct(x,name,id,arrayflag)
      if nargin == 1 && isstruct(x)
        y.id        = x.id;
        y.name      = x.name;
        y.val       = x.val;
        y.arrayflag = x.arrayflag;
      else
        y.id   = id;
        y.name = name;
        y.val  = x;
        y.arrayflag = arrayflag;
      end
    end
  end
  
  % Overloaded Methods
  methods
    y = ctranspose(x)
    y = horzcat(varargin)
    z = isequal(x,y,varargin)
    y = length(x)
    function y = numel(varargin)
      y = 1;
    end
    yi = ppval(pp,xi)
    y = repmat(x,varargin)
    y = reshape(x,varargin)
    varargout = size(x,varargin)
    y = struct(varargin)
    x = subsasgn(x,s,b)
    y = subsref(x,s)
    y = transpose(x)
    y = vertcat(varargin)
  end
  
  % Overloaded Utility Methods
  methods
    adigatorPrintOutputIndices(x);
    x = adigatorStructAnalyzer(x,xStr,subsflag)
    x = adigatorVarAnalyzer(FunString,x,xStr,subsflag)
    flag = cadaCheckForDerivs(x)
    xOut = cadaOverMap(x,DAid)
    y = cadaPrintReMap(x,y,varID)
    z = cadaUnionVars(x,y)
  end
  
  % Access Methods
  methods
    function val = cadaGetStruct(y)
      val = y.val;
    end
    function val = cadaGetStructID(y)
      val = y.id;
    end
    function aflag = cadaIsArray(x)
      aflag = x.arrayflag;
    end
    function [val, name, id, arrayflag] = cadastructDecomp(x)
      val       = x.val;
      name      = x.name;
      id        = x.id;
      arrayflag = x.arrayflag;
    end
  end
  
end
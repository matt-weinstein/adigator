classdef cadastruct
  % cadastruct classdef file for use with the ADiGator package.
  %
  % This class was introduced in order to adigator to handle cell and
  % structure array references/assignments. It is used internally by the
  % ADiGator algorithm and is not intended for general use.
  %
  % Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
  % Distributed under the GNU General Public License version 3.0
  
  properties
    % Unique integer ID
    id 
    % Name of the cell/structure in the printed derivative file
    name 
    % Structure/cell data containing CADA objects
    val
    % Flag stating whether data is cell/structure array or scalar structure
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
    yi = ppval(pp,xi)
    y = repmat(x,varargin)
    y = reshape(x,varargin)
    varargout = size(x,varargin)
    y = struct(varargin)
    x = subsasgn(x,s,b)
    y = subsref(x,s)
    y = transpose(x)
    y = vertcat(varargin)
    function y = isfield(x,f)
      % CADASTRUCT overloaded ISFIELD function
      y = isfield(x.val,f);
    end
    function y = numel(varargin)
      % CADASTRUCT overloaded NUMEL - always returns 1, cannot be used to
      % determine number of elements of cell/structure array
      y = 1;
    end
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
      % Value access for cadastruct
      val = y.val;
    end
    function val = cadaGetStructID(y)
      % ID access for cadastruct
      val = y.id;
    end
    function aflag = cadaIsArray(x)
      % Checks if object corresponds to a structure/cell array
      aflag = x.arrayflag;
    end
    function [val, name, id, arrayflag] = cadastructDecomp(x)
      % Accesses all cadastruct data
      val       = x.val;
      name      = x.name;
      id        = x.id;
      arrayflag = x.arrayflag;
    end
  end
  
end
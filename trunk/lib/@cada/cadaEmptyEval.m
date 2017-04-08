function [y,varargout] = cadaEmptyEval(varargin)
% This serves to simply increase ADIGATOR.VARINFO.COUNT and assign stuff
% to ADIGATOR.VARINFO.LASTOCC when in empty evaluations
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
global ADIGATOR
for Vcount = 1:nargin
  x = varargin{Vcount};
  if isa(x,'cada')
    ADIGATOR.VARINFO.LASTOCC(x.id,1) = ADIGATOR.VARINFO.COUNT+nargout-1;
  end
end
NUMvod = ADIGATOR.NVAROFDIFF;
y.id    = ADIGATOR.VARINFO.COUNT;
y.func  = struct('name',[],'size',[0 0],'zerolocs',[],'value',[]);
y.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
ADIGATOR.VARINFO.LASTOCC(ADIGATOR.VARINFO.COUNT,1) =...
  ADIGATOR.VARINFO.COUNT+nargout-1;
ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;
y = cada(y);

varargout = cell(nargout-1,1);
for Ocount = 1:nargout-1
  y1.id    = ADIGATOR.VARINFO.COUNT;
  y1.func  = struct('name',[],'size',[0 0],'zerolocs',[],'value',[]);
  y1.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
  ADIGATOR.VARINFO.LASTOCC(ADIGATOR.VARINFO.COUNT,1) =...
    ADIGATOR.VARINFO.COUNT+nargout-1-Ocount;
  ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;
  varargout{Ocount} = cada(y1);
end
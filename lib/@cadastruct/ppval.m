function yi = ppval(pp,xi)
% CADASTRUCT overloaded PPVAL
% This is just a wrapper for the CADA ppval.
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
global ADIGATOR
if ~isa(xi,'cada')
  NUMvod  = ADIGATOR.NVAROFDIFF;
  PFLAG   = ADIGATOR.PRINT.FLAG;
  func = struct('name',[],'size',size(xi),'zerolocs',[],...
    'value',xi);
  deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
  if PFLAG && ~ADIGATOR.EMPTYFLAG
    if numel(xi) < 11 && floor(xi(:))==xi(:)
      func.name = mat2str(xi);
    else
      func.name = cadamatprint(xi);
    end
  end
  xi = cada([],func,deriv);
end

if pp.arrayflag
  ADIGATOR.VARINFO.LASTOCC(pp.id,1) = ADIGATOR.VARINFO.COUNT;
  ADIGATOR.VARINFO.NAMELOCS(pp.id,3) = -Inf;
  fnames = fieldnames(pp.val);
  ppname = pp.name;
  for Fcount = 1:length(fnames)
    F = fnames{Fcount};
    v = pp.val.(F);
    if isa(v,'cada')
      v.func.name = [ppname,'.',F];
    end
  end
else
  fnames = fieldnames(pp.val);
  for Fcount = 1:length(fnames)
    F = fnames{Fcount};
    v = pp.val.(F);
    if isa(v,'cada')
      vid = v.id;
      ADIGATOR.VARINFO.LASTOCC(vid,1) = ADIGATOR.VARINFO.COUNT;
      ADIGATOR.VARINFO.NAMELOCS(vid,3) = -Inf;
    end
  end
end
yi = ppval(pp.val,xi);


end
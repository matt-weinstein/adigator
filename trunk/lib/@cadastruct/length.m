function y = length(x)
% CADASTRUCT overloaded LENGTH
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
global ADIGATOR
% just make a dummy CADA variable and use CADA length
NUMvod = ADIGATOR.NVAROFDIFF;
func  = struct('name',x.name,'size',size(x.val),'zerolocs',[],'value',[]);
deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));

xdummy = cada(x.id,func,deriv);
y = length(xdummy);
end
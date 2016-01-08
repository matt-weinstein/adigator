function display(x)
% Overloaded CADA display function
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0

fprintf('    id:')
disp(x.id)
fprintf('    func:\n');
disp(x.func)
fprintf('    deriv:\n');
disp(x.deriv)
function qd = getqd(t)
qd = zeros(2,length(t));
 qd(1,:) = 1.57*sin(6*t).*(1-exp(-.05*t.^4));
 qd(2,:) = 1.2*sin(1*t) .*(1-exp(-.05*t.^2));
% qd(1,:) = 0.7.*sin(2.*t).*(1-exp(-0.3.*t.^3));
% qd(2,:) = qd(1,:);
end
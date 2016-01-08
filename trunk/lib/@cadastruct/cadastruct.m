function y = cadastruct(x,name,id,arrayflag)
y.id   = id;
y.name = name;
y.val  = x;
y.arrayflag = arrayflag;
y = class(y,'cadastruct');
end
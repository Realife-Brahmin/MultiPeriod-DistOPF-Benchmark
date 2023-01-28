
function f = objfun(x, R, T, nb, Table)
% global Table;
% global nb;
% global T;
% global R;
%% onjective function I^2*R

f=0;
% f = x(1);
for i =2:nb
    Parent = find(i == T(:,2));
    
    if isempty(Parent)  
       f = f;
    else
        f =f+x((Table(Parent,5)))*R(T(Parent,1),T(Parent,2));
    end
    
end

 
 
 
    

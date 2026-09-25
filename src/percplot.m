function [pp] = percplot(x,Y,c,tr,lst,ec)

if size(x,2)==1
    x=x';
end
if size(Y,2)==3
    Y=Y(:,1:2:3);
elseif size(Y,1)==3
    Y=Y(1:2:3,:);
end
if size(Y,1)~=2
    Y=Y';
    if size(Y,1)~=2
        error('Expected Y to contain a percentile range: size(Y)=2(3)xN or Nx2(3)');
    end
end
 
xf=[x fliplr(x)];
Yf=[Y(2,:) fliplr(Y(1,:))];
 
pp=fill(xf,Yf,c,'LineStyle',lst,'EdgeColor',ec);
alpha(tr);
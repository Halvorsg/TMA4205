function img = mycolorwheel(n)

x = (-n:n)/n;
y = (-n:n)/n;

[XX,YY] = meshgrid(x,y);

circle = (XX.^2+YY.^2 <1);

UU = XX.*circle;
VV = YY.*circle;

img = mycomputeColor(UU,VV);
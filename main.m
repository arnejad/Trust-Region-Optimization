% This  code has been implemeted by Ali Nejad for educational purposes

f =     @(x1,x2) (1-x1).^2 + (10 * ((x2-x1.^2).^2));
g1 =    @(x1,x2) -40*x1*(x2-(x1^2)) - 2*(1 - x1);
g2 =    @(x1,x2) 20*(x2-x1^2);
H =     @(x1, x2) [-40*(x2)+(120*(x1^2))+2 -40*x1; -40*x1 20];

f_drawable = @(x1,x2) (ones(401)-x1).^2 + (10 * ((x2-x1.^2).^2));
% res = trustRegion(f, g1, g2, H, [2;2], "cauchyPoint", f_drawable, [ -1,3;-1,3]);
res = trustRegion(f, g1, g2, H, [2;2], "dogleg", f_drawable, [ -1,3;-1,3]);
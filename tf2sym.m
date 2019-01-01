function tfs = tf2sym(tf,n,m)
 syms s
 [num, den] = tfdata(tf(n,m));
 tfs = vpa(poly2sym(cell2mat(num),s)/poly2sym(cell2mat(den),s),4);
%  pretty(symb)
end
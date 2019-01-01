function[ans] = syms2tf(G)
 [symNum,symDen] = numden(G); %Get num and den of Symbolic TF
 TFnum = sym2poly(symNum); %Convert Symbolic num to polynomial
    TFden = sym2poly(symDen); %Convert Symbolic den to polynomial
 ans =tf(TFnum,TFden);
end
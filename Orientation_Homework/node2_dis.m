function s2 = node2_dis(r,F,E,L)
   [Q,sigma,R] = FEA(r,F,E,L);
   s2 = abs(Q(2));
end
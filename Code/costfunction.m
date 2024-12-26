function cost = costfunction(f1, f2, f3, f4, f5, u1, u2, u3)
syms u a1 b1 a2 b2 a3 b3 a4 b4;
f1approx=a1+b1*u;
f2approx=a2+b2*u;
f3approx=a3+b3*u;
f4approx=a4+b4*u;
cost=int((f1-f1approx)^2,0, 2)+int((f2-f1approx)^2, 2, 5)+int((f3-f2approx)^2, 5, 6.5)+int((f3-f3approx)^2, 6.5, 7)+int((f4-f3approx)^2, 7, 9)+int((f5-f3approx)^2, 9, 11)+int((f5-f4approx)^2, 11, 15);
end
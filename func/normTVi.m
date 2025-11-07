function val = normTVi(x)

grad = grads(x);
val = sum(sqrt(sum(grad.^2,4)),[1,2,3]);

end


function h = xs_response(dict, alpha, lam, j)
%XS_RESPONSE  Response h(lambda; alpha) of network j on its spectrum (DC mode 0).
h = dict.Phi(lam, j) * alpha;
h(lam(:) == 0) = 0;
end

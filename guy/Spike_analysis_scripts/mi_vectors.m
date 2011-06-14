function [mi_vector, mi_dec_vector, mi_prec_vector, mi_dec_prec_vector, separation_vector] = mi_vectors(scale, mixing)

mi_vector = [];
mi_dec_vector = [];
mi_prec_vector = [];
mi_dec_prec_vector = [];
separation_vector = [];


for bias = 20:-10:-20
    load(sprintf('/home/ucbtepi/code/network/data/f.5_20_-20/s%.2f/c%.2f/result_b%02d.mat', scale, mixing, bias));
    mi_vector = [mi_vector, mi(:,2)];
    mi_dec_vector = [mi_dec_vector, mi_dec(2:end)'];
    mi_prec_vector = [mi_prec_vector, mi_prec];
    mi_dec_prec_vector = [mi_dec_prec_vector, mi_dec_prec'];
    separation_vector = [separation_vector, separation'];
end
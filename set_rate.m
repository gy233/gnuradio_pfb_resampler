function [d_dec_rate,d_flt_rate]=set_rate(rate,d_int_rate)
d_dec_rate = floor(d_int_rate/rate);
d_flt_rate = (d_int_rate/rate)-d_dec_rate;
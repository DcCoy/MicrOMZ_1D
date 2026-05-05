function y = lin_profile(z, zt, zb, topval, botval)
y = topval + (botval - topval) * (z - zt)/(zb - zt);
end

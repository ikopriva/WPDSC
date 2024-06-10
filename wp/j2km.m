function [m, k] = j2km(j)
% Two dimensional offset [k, m] is calculated from linear offset j.
% Example: conversion table for decomposition level n = 2 
%         k=0  1  2  3
%
% m=0     j=0  1  4  5
%   1       2  3  6  7
%   2       8  9 12 13
%   3      10 11 14 15 
%
 j = uint16(j);
 k = uint16(0);
 k = bitset(k, 1, bitget(j,1));  k = bitset(k, 2, bitget(j,3));  k = bitset(k, 3, bitget(j,5));  k = bitset(k, 4, bitget(j,7));
 k = bitset(k, 5, bitget(j,9));  k = bitset(k, 6, bitget(j,11)); k = bitset(k, 7, bitget(j,13)); k = bitset(k, 8, bitget(j,15));
 m = uint16(0);
 m = bitset(m, 1, bitget(j,2));  m = bitset(m, 2, bitget(j,4));  m = bitset(m, 3, bitget(j,6));  m = bitset(m, 4, bitget(j,8));
 m = bitset(m, 5, bitget(j,10)); m = bitset(m, 6, bitget(j,12)); m = bitset(m, 7, bitget(j,14)); m = bitset(m, 8, bitget(j,16));
 k = double(k); m = double(m);
end

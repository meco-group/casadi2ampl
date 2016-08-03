#how to make auxiliary variables

reset;

param n:=5;
var x{1..n};
var cx{i in 1..n} = cos(x[i]);

minimize f: sum{i in 1..n}cx[i];

s.t.
	con{i in 1..n}: x[i]<=4;

option solver knitro;

solve;

display x;

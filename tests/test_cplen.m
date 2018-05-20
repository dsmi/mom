function test_cplen
%
% Test of cperlen function, it is called this way because of the name conflict.
% 
% Compare for a two-conductor case.
%

C2 = [ 1e-010  -8e-011 ;  -5e-011  2e-010 ];

C = cperlen(C2);

C_test = C2(1,1) - sum(C2(1,:))*sum(C2(:,1))/sum(sum(C2));

assertEquals(C_test, C, 1e-25);


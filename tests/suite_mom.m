function suite_mom
% suite_mom : testsuite which joins all the unittests of the
% mom solvers.
%
addTest('test_intg_helmsl2d');
addTest('test_intg_helmdl2d');
addTest('test_intg_laplace2d');
addTest('test_intg_lapdn2d');
addTest('test_mkmommat2');
addTest('test_find_edges2d');
addTest('test_find_eol2d');
addTest('test_ports2subs');
addTest('test_cperlen');
addTest('test_cperlen3');
addTest('test_cdiel');
addTest('test_roundz');
addTest('test_spiralr2d');
addTest('test_stripr2d');
addTest('test_cplen');
addTest('test_mkbases2d');
addTest('test_intg_lapsl2l');
addTest('test_intg_lapdn2l');
addTest('test_load_mesh');
addTest('test_intg_lapdn2s');
addTest('test_mstripc');

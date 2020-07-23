# Taken from https://travis-ci.com/github/fusion809/FunctionIntegrator.jl/jobs/364132770
# 10th entry is simppen
L_adaptive_simpsons = [0.486231; 0.064895; 0.000165; 0.037012; 0.032331; 0.030399; 0.023287; 0.023555; 0.099551; 9.451517; 0.039571; 0.000751; 0.000284];
L_chebyshev1 = [0.999828; 0.378005; 0.044058; 0.080038; 0.098033; 0.088635; 0.083370; 0.083933; 0.060513; 0.329102; 0.085079; 0.046761; 0.043942];
L_chebyshev2 = [0.122684; 0.060723; 0.043665; 0.116662; 0.095775; 0.125398; 0.100822; 0.099486; 0.062181; 8.654814; 0.085809; 0.049589; 0.044264];
L_chebyshev3 = [0.148326; 0.069086; 0.048887; 0.093295; 0.092409; 0.098670; 0.106251; 0.107016; 0.067388; 3.168286; 0.084671; 0.054494; 0.049837];
L_chebyshev4 = [0.197283; 0.067486; 0.048484; 0.109152; 0.105551; 0.098708; 0.093042; 0.091734; 0.084706; 4.084928; 0.084922; 0.052913; 0.049484];
L_jacobi = [20.378879; 0.535539; 0.078133; 0.481214; 0.225388; 2.794902; 0.530739; 0.424725; 0.459931; 21.860716; 0.155370; 1.836021; 0.575570];
L_legendre = [0.101542; 0.070673; 0.030698; 0.075401; 0.099265; 0.075489; 0.070451; 0.075890; 0.040959; 4.348325; 0.056232; 0.033010; 0.030431];
L_lobatto = [0.178613; 0.011121; 0.000107; 0.094005; 0.119870; 0.115679; 0.082939; 0.077295; 0.000196; 21.474618; 0.026952; 0.000390; 0.000152];
L_radau = [0.134821; 0.011121; 0.000051; 0.081257; 0.124391; 0.103346; 0.084129; 0.078259; 0.000089; 22.084138; 0.024211; 0.000374; 0.000090];
L_rectangle_midpoint = [0.074212; 0.000018; 0.000260; 0.018878; 0.020378; 0.019273; 0.018416; 0.018761; 0.000029; 1.586260; 0.028527; 0.000580; 0.000712];
L_rombergs = [0.458671; 0.224695; 0.042407; 0.070613; 0.068976; 0.073960; 0.085664; 0.072760; 0.066720; 0.558451; 0.137821; 0.038345; 0.043788];
L_simpsons = [0.024996; 0.000024; 0.000014; 0.020277; 0.030296; 0.022534; 0.020241; 0.023021; 0.000033; 13.500102; 0.001098; 0.000061; 0.000051];
L_simpsons38 = [0.053542; 0.040926; 0.000011; 0.023485; 0.033142; 0.025603; 0.023186; 0.038315; 0.018009; 13.818143; 0.041548; 0.000095; 0.000150];
L_trapezoidal = [0.150540; 0.035154; 0.000469; 0.016305; 0.022002; 0.017666; 0.015459; 0.017220; 0.012025; 12.852069; 0.082628; 0.001076; 0.001369];
N = length(L_simpsons);

# RMS of times
rms_adaptive_simpsons_wo_simppen = sqrt((L_adaptive_simpsons[1:9]'*L_adaptive_simpsons[1:9]+L_adaptive_simpsons[11:N]'*L_adaptive_simpsons[11:N])/(N-1));
rms_chebyshev1 = sqrt(L_chebyshev1'*L_chebyshev1/N);
rms_chebyshev1_wo_simppen = sqrt((L_chebyshev1[1:9]'*L_chebyshev1[1:9]+L_chebyshev1[11:N]'*L_chebyshev1[11:N])/(N-1));
rms_chebyshev2 = sqrt(L_chebyshev2'*L_chebyshev2/N);
rms_chebyshev2_wo_simppen = sqrt((L_chebyshev2[1:9]'*L_chebyshev2[1:9]+L_chebyshev2[11:N]'*L_chebyshev2[11:N])/(N-1));
rms_chebyshev3 = sqrt(L_chebyshev3'*L_chebyshev3/N);
rms_chebyshev3_wo_simppen = sqrt((L_chebyshev3[1:9]'*L_chebyshev3[1:9]+L_chebyshev3[11:N]'*L_chebyshev3[11:N])/(N-1));
rms_chebyshev4 = sqrt(L_chebyshev4'*L_chebyshev4/N);
rms_chebyshev4_wo_simppen = sqrt((L_chebyshev4[1:9]'*L_chebyshev4[1:9]+L_chebyshev4[11:N]'*L_chebyshev4[11:N])/(N-1));
rms_jacobi = sqrt(L_jacobi'*L_jacobi/N);
rms_jacobi_wo_simppen = sqrt((L_jacobi[1:9]'*L_jacobi[1:9]+L_jacobi[11:N]'*L_jacobi[11:N])/(N-1));
rms_legendre = sqrt(L_legendre'*L_legendre/N);
rms_legendre_wo_simppen = sqrt((L_legendre[1:9]'*L_legendre[1:9]+L_legendre[11:N]'*L_legendre[11:N])/(N-1));
rms_lobatto_wo_simppen = sqrt((L_lobatto[1:9]'*L_lobatto[1:9]+L_lobatto[11:N]'*L_lobatto[11:N])/(N-1));
rms_radau_wo_simppen = sqrt((L_radau[1:9]'*L_radau[1:9]+L_radau[11:N]'*L_radau[11:N])/(N-1));
rms_rectangle_midpoint = sqrt(L_rectangle_midpoint'*L_rectangle_midpoint/N);
rms_rectangle_midpoint_wo_simppen = sqrt((L_rectangle_midpoint[1:9]'*L_rectangle_midpoint[1:9]+L_rectangle_midpoint[11:N]'*L_rectangle_midpoint[11:N])/(N-1));
rms_rombergs_wo_simppen = sqrt((L_rombergs[1:9]'*L_rombergs[1:9]+L_rombergs[11:N]'*L_rombergs[11:N])/(N-1));
rms_simpsons_wo_simppen = sqrt((L_simpsons[1:9]'*L_simpsons[1:9]+L_simpsons[11:N]'*L_simpsons[11:N])/(N-1));
rms_simpsons38_wo_simppen = sqrt((L_simpsons38[1:9]'*L_simpsons38[1:9]+L_simpsons38[11:N]'*L_simpsons38[11:N])/(N-1));
rms_trapezoidal_wo_simppen = sqrt((L_trapezoidal[1:9]'*L_trapezoidal[1:9]+L_trapezoidal[11:N]'*L_trapezoidal[11:N])/(N-1));

rms_wo_simppen = [rms_chebyshev1_wo_simppen; rms_chebyshev2_wo_simppen; rms_chebyshev3_wo_simppen; rms_chebyshev4_wo_simppen; rms_jacobi_wo_simppen; rms_legendre_wo_simppen; rms_lobatto_wo_simppen; rms_radau_wo_simppen; rms_rectangle_midpoint_wo_simppen; rms_simpsons_wo_simppen; rms_trapezoidal_wo_simppen];

println("rms_chebyshev1 is $(rms_chebyshev1)")
println("rms_chebyshev2 is $(rms_chebyshev2)")
println("rms_chebyshev3 is $(rms_chebyshev3)")
println("rms_chebyshev4 is $(rms_chebyshev4)")
println("rms_jacobi     is $(rms_jacobi)")
println("rms_legendre   is $(rms_legendre)")

println("rms_adaptive_simpsons_wo_simppen,$(rms_adaptive_simpsons_wo_simppen)")
println("rms_chebyshev1_wo_simppen,$(rms_chebyshev1_wo_simppen)")
println("rms_chebyshev2_wo_simppen,$(rms_chebyshev2_wo_simppen)")
println("rms_chebyshev3_wo_simppen,$(rms_chebyshev3_wo_simppen)")
println("rms_chebyshev4_wo_simppen,$(rms_chebyshev4_wo_simppen)")
println("rms_jacobi_wo_simppen,$(rms_jacobi_wo_simppen)")
println("rms_legendre_wo_simppen,$(rms_legendre_wo_simppen)")
println("rms_lobatto_wo_simppen,$(rms_lobatto_wo_simppen)")
println("rms_radau_wo_simppen,$(rms_radau_wo_simppen)")
println("rms_rectangle_midpoint_wo_simppen,$(rms_rectangle_midpoint_wo_simppen)")
println("rms_rombergs_wo_simppen,$(rms_rombergs_wo_simppen)")
println("rms_simpsons_wo_simppen,$(rms_simpsons_wo_simppen)")
println("rms_simpsons38_wo_simppen,$(rms_simpsons38_wo_simppen)")
println("rms_trapezoidal_wo_simppen,$(rms_trapezoidal_wo_simppen)")

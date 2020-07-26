# Taken from https://travis-ci.com/github/fusion809/FunctionIntegrator.jl/jobs/365077894
# 10th entry is simppen
L_adaptive_simpsons = [0.455044; 0.051227; 0.000142; 0.037223; 0.031272; 0.027268; 0.032873; 0.024794; 0.092391; 11.442810; 0.040696; 0.001140; 0.000531];
L_chebyshev1 = [1.009558; 0.335515; 0.044585; 0.081733; 0.091252; 0.087307; 0.082136; 0.102419; 0.060185; 0.339226; 0.086334; 0.046033; 0.043327];
L_chebyshev2 = [0.132044; 0.066955; 0.045266; 0.095766; 0.112055; 0.103432; 0.098146; 0.110522; 0.063643; 8.680677; 0.092204; 0.048493; 0.042964];
L_chebyshev3 = [0.220306; 0.071286; 0.051752; 0.137612; 0.093434; 0.116630; 0.104549; 0.099629; 0.067163; 3.208961; 0.091863; 0.053932; 0.048825];
L_chebyshev4 = [0.136945; 0.070568; 0.051149; 0.092509; 0.092965; 0.098516; 0.093687; 0.112979; 0.065209; 4.127343; 0.091688; 0.052680; 0.048865];
L_jacobi = [20.866783; 0.489243; 0.078353; 0.543790; 0.219776; 2.871823; 0.514692; 0.425720; 0.474459; 21.683630; 0.164953; 1.855387; 0.561469];
L_legendre = [0.106171; 0.078795; 0.031839; 0.074370; 0.074632; 0.077306; 0.074217; 0.069197; 0.041277; 4.304367; 0.059233; 0.032885; 0.030395];
L_lobatto = [0.190303; 0.012930; 0.000149; 0.085664; 0.140379; 0.115814; 0.095905; 0.090178; 0.000142; 21.309380; 0.026758; 0.000402; 0.000159];
L_radau = [0.141603; 0.022506; 0.000051; 0.086428; 0.110123; 0.099484; 0.090253; 0.079221; 0.000128; 21.947621; 0.024431; 0.000389; 0.000146];
L_rectangle_midpoint = [0.073207; 0.000017; 0.000273; 0.018856; 0.020056; 0.019106; 0.018988; 0.018468; 0.000033; 1.635390; 0.029941; 0.000537; 0.000739];
L_rombergs = [0.502167; 0.249403; 0.043215; 0.076397; 0.084500; 0.065737; 0.083264; 0.070981; 0.062676; 0.559453; 0.139889; 0.037316; 0.041226];
L_simpsons = [0.025735; 0.000020; 0.000016; 0.021620; 0.029775; 0.023050; 0.024037; 0.022381; 0.000023; 13.338188; 0.001048; 0.000050; 0.000041];
L_simpsons38 = [0.055831; 0.042373; 0.000011; 0.025192; 0.031424; 0.025814; 0.025907; 0.024228; 0.016936; 13.800786; 0.039795; 0.000060; 0.000034];
L_trapezoidal = [0.152162; 0.036755; 0.000453; 0.017747; 0.021522; 0.018655; 0.018124; 0.015779; 0.012151; 12.730575; 0.082469; 0.000966; 0.001362];
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

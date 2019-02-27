import numpy as np
from legendre_to_interpolating import convert_to_interpolating

in_file = open("I_b-spline-deriv1.txt", "r")
out_file = open('L_b-spline-deriv1.txt', 'w')

# Input files with the following style:
#
# 2
# -1.131377801839763e-01 -6.219856949725088e-01
# 2.811997871134302e-01 -7.276481276751023e-01
# -3.000187374447446e-01 5.886938720395576e-01
# -5.886938720395574e-01 3.000187374447447e-01
# 7.276481276751025e-01 -2.811997871134304e-01
# 6.219856949725088e-01 1.131377801839761e-01
# 3
# 1.601834450047005e-01 -7.046051955774445e-01 -3.431687211510057e-01
# -9.630562745712626e-02 4.343434343434349e-01 -9.777508000495582e-01
# -6.488830489667494e-02 2.893566139391730e-01 -2.950557118863133e-01
# -5.351722243945190e-02 4.979592136259431e-01 5.039576531510103e-01
# -1.159796279726884e+00 -3.384671404501071e-16 1.159796279726884e+00
# -5.039576531510097e-01 -4.979592136259433e-01 5.351722243945141e-02
# 2.950557118863135e-01 -2.893566139391732e-01 6.488830489667509e-02
# 9.777508000495581e-01 -4.343434343434345e-01 9.630562745712651e-02
# 3.431687211510054e-01 7.046051955774449e-01 -1.601834450047005e-01
# ....


conv_to_interpolating = False

for i in range(2, 21):
    check = in_file.readline()
    if (int(check) != i):
        print("Wrong order ", "check: ", check, "i: ", i)
    k0 = np.empty([i, i])
    k1 = np.empty([i, i])
    k2 = np.empty([i, i])
    for j in range(i):
        a = in_file.readline()
        k0[j, :] = np.fromstring(a, dtype=np.float64, sep=' ')
    for j in range(i):
        a = in_file.readline()
        k1[j, :] = np.fromstring(a, dtype=np.float64, sep=' ')
    for j in range(i):
        a = in_file.readline()
        k2[j, :] = np.fromstring(a, dtype=np.float64, sep=' ')
    k0I = convert_to_interpolating(i, k0, conv_to_interpolating)  # Add kwarg True to go
    k1I = convert_to_interpolating(i, k1, conv_to_interpolating)  # from interpolating to the
    k2I = convert_to_interpolating(i, k2, conv_to_interpolating)  # Legendre
    out_file.write(str(i) + '\n')
    for j in range(i):
        st = np.array2string(k0I[j,:], formatter={'float_kind':lambda x: "%.15e" % x}, max_line_width=10e10)
        st1 = st.replace('[', '')
        st2 = st1.replace(']', '')
        out_file.write(st2+'\n')
    for j in range(i):
        st = np.array2string(k1I[j,:], formatter={'float_kind':lambda x: "%.15e" % x}, max_line_width=10e10)
        st1 = st.replace('[', '')
        st2 = st1.replace(']', '')
        out_file.write(st2 + '\n')
    for j in range(i):
        st = np.array2string(k2I[j,:], formatter={'float_kind':lambda x: "%.15e" % x}, max_line_width=10e10)
        st1 = st.replace('[', '')
        st2 = st1.replace(']', '')
        out_file.write(st2 + '\n')
out_file.close()

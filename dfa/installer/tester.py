import numpy as np
import dfa
import matplotlib.pyplot as plt
import scipy.signal as ss
np.random.seed(42)

print("1. Test")
np.random.seed(42)
X = np.random.randn(1000)
X = np.abs(ss.hilbert(X))
#X = np.ones(1000)
#np.array([[1.,1.], [2., 1.], [3., 1.2], [4., 1.], [5., 0.8], [6., 1.]])
scales, fluct_to_array, coeff_to_array= dfa.play(X,5,9,0.25)
print(fluct_to_array)
print(coeff_to_array)
#we need to swap the order here, due to the way coeff are computed!
coeff_tmp = coeff_to_array[1]
coeff_to_array[1] = coeff_to_array[0]
coeff_to_array[0] = coeff_tmp
fluctfit = 2**np.polyval(coeff_to_array,np.log2(scales))
plt.loglog(scales, fluct_to_array, 'bo')
plt.loglog(scales, fluctfit, 'r', label=r'$\alpha$ = %0.2f'%coeff_to_array[0])
plt.title('DFA')
plt.xlabel(r'$\log_{10}$(time window)')
plt.ylabel(r'$\log_{10}$<F(t)>')#
plt.legend()
plt.show()

"""
nmf.play(X,2,100, 42,1)
#X, n_components, max_iter, random_state, verbose
print("2. test")
X = np.random.rand(10,10)
nmf.play(X,2,100,42,1)
#then try with 6 components and the same with 95 columns
print("3. test")
X = np.random.rand(1025,2)
nmf.play(X,2, 100, 42,1)
print("4.test")
X=np.random.rand(1025,90)
nmf.play(X, 2, 100, 42,1)
print("5. test")
X=np.random.rand(1025,90)
nmf.play(X, 6, 100, 42,1)
"""

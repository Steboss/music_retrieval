import numpy as np
import nmf
print("1. Test")
#X = np.array([[1.,1.], [2., 1.], [3., 1.2], [4., 1.], [5., 0.8], [6., 1.]])
#nmf.play(X,2,1000, 42,1)
#X, n_components, max_iter, random_state, verbose
print("2. test")
#make a 1025 x 95 array
X = np.random.rand(1025,95)
#now pass this to th enmf  with 2 components
nmf.play(X,2,1000,42,0)
#then try with 6 components

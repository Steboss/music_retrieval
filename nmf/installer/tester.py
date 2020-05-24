import numpy as np
import nmf
np.random.seed(42)

print("1. Test")
X = np.array([[1.,1.], [2., 1.], [3., 1.2], [4., 1.], [5., 0.8], [6., 1.]])
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

#modify files in cd /home/steboss/.local/lib/python3.6/site-packages/sklearn/decomposition/_nmf.py

import numpy as np
X = np.array([[1,1], [2, 1], [3, 1.2], [4, 1], [5, 0.8], [6, 1]])
#X = np.random.rand(1025,10)
from sklearn.decomposition import non_negative_factorization
W, H, n_iter = non_negative_factorization(X, n_components=2, init='random', random_state=42, max_iter=100)

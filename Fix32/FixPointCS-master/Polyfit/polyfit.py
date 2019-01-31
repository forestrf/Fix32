#
# FixPointCS
#
# Copyright(c) 2018 Jere Sanisalo, Petri Kero
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#

# OBSOLETE -- remez.py is used instead

# Utility script for generating fitted polynomials to approximate various
# trascendental functions for use in Fixed64.cs.
#
# The current algorithm optimize least squares error (LSE), which produces
# good results. Better polynomials could be generated by optimizing for the
# worst case error, for example with the Remez algorithm.
#
# The generated code may require some modifications like renaming the input
# variable or dropping some zero constants. Also, not all the generated
# polynomials can be calculated with s2.30 due to overflows.

import numpy as np
import numpy.polynomial.polynomial as poly
import matplotlib.pyplot as plt

FIT_SAMPLES = 500     # number of samples to use for fitting polynomials
ERR_SAMPLES = 50000   # number of samples to use for computing error in fitted polynomials

# see: https://stackoverflow.com/questions/15191088/how-to-do-a-polynomial-fit-with-fixed-points
def polyfitWithFixedPoints(n, x, y, xf, yf):
  mat = np.empty((n + 1 + len(xf),) * 2)
  vec = np.empty((n + 1 + len(xf),))
  x_n = x**np.arange(2 * n + 1)[:, None]
  yx_n = np.sum(x_n[:n + 1] * y, axis=1)
  x_n = np.sum(x_n, axis=1)
  idx = np.arange(n + 1) + np.arange(n + 1)[:, None]
  mat[:n + 1, :n + 1] = np.take(x_n, idx)
  xf_n = xf**np.arange(n + 1)[:, None]
  mat[:n + 1, n + 1:] = xf_n / 2
  mat[n + 1:, :n + 1] = xf_n.T
  mat[n + 1:, n + 1:] = 0
  vec[:n + 1] = yx_n
  vec[n + 1:] = yf
  params = np.linalg.solve(mat, vec)
  return params[:n + 1]

def genQmul(inputName, ndx, order):
  if ndx == order:
    return f"Qmul30(C{ndx}, {inputName})"
  else:
    return f"Qmul30({genQmul(inputName, ndx+1, order)} + C{ndx}, {inputName})"

def fitFunc(name, order, func, domain, fixed, inputName):
  (mn, mx) = domain
  x = mn + (mx - mn) * (np.arange(FIT_SAMPLES * 1.0) / (FIT_SAMPLES - 1))
  y = func(x)

  # fit polynomial (with fixed points)
  coefs = polyfitWithFixedPoints(order, x, y, fixed, func(np.array(fixed)))
  fitted = poly.Polynomial(coefs)

  # compute max error
  ex = mn + (mx - mn) * (np.arange(ERR_SAMPLES * 1.0) / (ERR_SAMPLES - 1))
  ey = func(ex)
  err = ey - fitted(ex)
  errmax = np.max(err)

  # print polynomial error and code
  print(f'  order {order}: {-np.log2(errmax):0.4} bits')

  for ndx in range(len(coefs)):
    coef = coefs[ndx]
    if np.abs(coef) < 1e-12:
      coef = 0
    print(f"    const int C{ndx} = {int(coef * (1 << 30))}; // {coef}")
  
  print(f"    int y = {genQmul(inputName, 1, order)} + C0;")
  print()

  return fitted

FUNCS = {
  'sin': ((-1.0, 1.0), 'z', lambda x: np.sin(x * 0.5*np.pi), [-1.0, 0.0, 1.0]),
#  'cos': ((0.0, 1.0), 'z', lambda x: np.cos(x * 0.5*np.pi), [0.0, 1.0]),
  'atan': ((0.0, 1.0), 'k', lambda x: np.arctan(x), [0.0]),
  'rcpm1': ((0.0, 1.0), 'k', lambda x: 1.0 / (x + 1), [0.0, 1.0]),
  'sqrt': ((1.0, 2.0), 'n', lambda x: np.sqrt(x), [1.0, 2.0]),
  'rsqrtm1': ((0.0, 1.0), 'k', lambda x: 1.0 / np.sqrt(x + 1), [0.0, 1.0]),
  'logm1': ((0.0, 1.0), 'k', lambda x: np.log(x+1), [0.0, 1.0]),
  'exp2': ((0.0, 1.0), 'k', lambda x: np.exp2(x), [0.0, 1.0]),
}

for name in FUNCS.keys():
  if name != 'logm1':
    continue

  # extract func domain and function
  (domain, inputName, func, fixed) = FUNCS[name]
  (mn, mx) = domain

  print('')
  print(f'{name}({mn} .. {mx}):')

  for order in [3, 4, 5, 6, 7, 8, 9]:
    fitted = fitFunc(name, order, func, domain, fixed, inputName)
    if order >= 9:
      mne = mn + 0.0001
      x = mne + (mx - mne) * (np.arange(ERR_SAMPLES * 1.0) / (ERR_SAMPLES - 1))
      plt.plot(x, (fitted(x) - func(x)) / func(x))

  plt.show()

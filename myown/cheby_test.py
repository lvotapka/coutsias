'''
Created on May 26, 2017

@author: lvotapka

This module runs some preliminary tests to solve a 1D differential equation
using a spectral method according to the following paper:

"Integration preconditioners for differential operators in spectral tau-methods"
by Coutsias, Hagstrom, Hesthaven, and Torres 1995

I'm going to attempt to recreate the equation solvers outlined in the above
paper.
'''

import numpy as np
import numpy.polynomial.chebyshev as cheby
from numpy.polynomial import Polynomial as P
from numpy.polynomial import Chebyshev as T
import matplotlib.pyplot as plt

N = 64 # dimensionality of the vectors
N_bar = 128
a = -1.0
b = 1.0
epsilon = 1e-14 # allowable error



def make_u_true_soln_exp1(N, a, b, delta=5e-4, x0=0.0):
    '''Creates a vector u=u(x) that approximates the following function:
    
    u(x) = (1 / sqrt(delta)) * exp(-(x - x0)**2 / delta)
    
    Input:
      N: dimensionality of u vector
      a: left boundary of domain
      b: right boundary of domain
      delta: an effective variance of Gaussian-like function u
      x0: location of the maximum of function u
      
    Output:
      u: a vector that approximates function u(x)
    '''
    
    u = np.zeros(N)
    for i in range(N):
        x = a + i * (b-a)/N # assign location of x based on index i
        u[i] = (1 / np.sqrt(delta)) * np.exp(-(x - x0)**2 / delta) # compute value of function
    return u

def broydens_method(A_N, Y_N, a, b, epsilon=epsilon):
    '''Use Broyden's method for solving a differential equation.
    
    Input:
      A_N: a matrix representing some operator
      Y_N: a vector representing the rhs of the equation
      a: left boundary of the domain
      b: right boundary of the domain
      
    Output:
      v: the solution to the equation described
    '''
    max_iter = 100
    #v_0 = s_0 = r_0 = Y_N # initialize vectors
    v_k = s_k = r_k = Y_N # initialize vectors
    gamma_k = np.dot(s_k, s_k)
    k = 1
    while np.sqrt(np.dot(r_k, r_k)) >= epsilon:
        r_k = Y_N - np.dot(A_N * v_k)
        z_k_current = r_k
        for j in range(2,k):
            pass # stopped here...
    


def plot_1D_function(f, a, b, title="function"):
    '''Plots a function f(x) approximated by vector f.
    
    Input:
      f: a vector of values
      a: left boundary of domain
      b: right boundary of domain
      title: optional title of plot
      
    Output:
      None
    '''
    
    x = np.arange(a, b, (b-a)/N)
    plt.plot(x, f, 'g')
    plt.title(title)
    plt.xlim(a,b)
    plt.show()
    return

if __name__ == "__main__":
    x = np.linspace(a, b, N_bar)
    u_true = make_u_true_soln_exp1(N_bar, a, b)
    u_true_cheb = T.fit(x, u_true, N)
    print u_true_cheb
    
    
    '''
    p = P([1,2,3]) # a conventional power series
    print p.coef
    print p.domain
    print p.window
    print p+p
    print p-p
    print p*p
    print p//P([-1,1])
    print p%P([-1,1])
    x = np.arange(5)
    print p(x)
    x = np.arange(6).reshape(3,2)
    print p(x)
    print p(p)
    print p.roots()
    print p.integ()
    print p.integ(2)
    print p.integ(lbnd=-1)
    print p.integ(lbnd=-1, k=1)
    print p.deriv()
    print p.deriv(2)
    p = P.fromroots([1,2,3])
    print p
    print p.convert(kind=T)
    
    for i in range(5):
        t = T.basis(i)
        print P.cast(t)
    ''
    x = np.linspace(-2,2,100)
    for i in range(6): ax = plt.plot(x, T.basis(i)(x), lw=2, label="$T_%d$"%i)
    plt.legend(loc="upper left")
    plt.show()
    ''
    np.random.seed(12)
    x = np.linspace(0, 2*np.pi, 20)
    y = np.sin(x) + np.random.normal(scale=0.1, size=x.shape)
    p = T.fit(x, y, 5)
    plt.plot(x, y, 'o')
    xx, yy = p.linspace()
    plt.plot(xx, yy, lw=2)
    print p.domain
    print p.window
    plt.show()
    '''
    #plot_1D_function(u, a, b, title='Plot of goal function')
    
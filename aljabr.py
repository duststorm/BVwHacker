#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
:Authors:
    Manuel Bastioni,
    Marc Flerackers

:Version: 1.0
:Copyright: MakeHuman Team 2001-2011
:License: GPL3 

This module contains the most common 3D algebraic operations used in MakeHuman (http://www.makehuman.org/).
These are mostly the vector and matrix operations core to any 3D application. For efficiency and speed, all matrix
operation will be done thru flat arrays. Function with matrices as flat arrays are written with underscore "_", whilst
functions with matrices as list of lists will have the same name without the underscore.

The name is a tribute to *Al-jabr wa'l muqabalah* the most important paper of Mohammed ibn-Musa al-Khuwarizmi (VII - VIII sec d.C.)
The paper was so important that Al-jabr is the root of modern word *algebra* and al-Khuwarizmi is the root of word *algorithm*.

Categories:
  Vector Operations
  Matrix Operations
  Quaternions
  Geometric Operations
  Various Functions
  
"""

from math import sqrt, cos, sin, tan, atan2, fabs, acos, pi, exp
from random import random

machine_epsilon = 1.0e-16
degree2rad = pi/180.0

"""

.. note::

    A triple of Euler angles can be applied/interpreted in 24 ways, which can
    be specified using a 4 character string or encoded 4-tuple:

      *Axes 4-string*: e.g. 'sxyz' or 'ryxy'

      - first character : rotations are applied to 's'tatic or 'r'otating frame
      - remaining characters : successive rotation axis 'x', 'y', or 'z'

      *Axes 4-tuple*: e.g. (0, 0, 0, 0) or (1, 1, 1, 1)

      - inner axis: code of axis ('x':0, 'y':1, 'z':2) of rightmost matrix.
      - parity : even (0) if inner axis 'x' is followed by 'y', 'y' is followed
        by 'z', or 'z' is followed by 'x'. Otherwise odd (1).
      - repetition : first and last axis are same (1) or different (0).
      - frame : rotations are applied to static (0) or rotating (1) frame.
"""
_NEXT_AXIS = [1, 2, 0, 1]
_AXES2TUPLE = {
    'sxyz': (0, 0, 0, 0), 'sxyx': (0, 0, 1, 0), 'sxzy': (0, 1, 0, 0),
    'sxzx': (0, 1, 1, 0), 'syzx': (1, 0, 0, 0), 'syzy': (1, 0, 1, 0),
    'syxz': (1, 1, 0, 0), 'syxy': (1, 1, 1, 0), 'szxy': (2, 0, 0, 0),
    'szxz': (2, 0, 1, 0), 'szyx': (2, 1, 0, 0), 'szyz': (2, 1, 1, 0),
    'rzyx': (0, 0, 0, 1), 'rxyx': (0, 0, 1, 1), 'ryzx': (0, 1, 0, 1),
    'rxzx': (0, 1, 1, 1), 'rxzy': (1, 0, 0, 1), 'ryzy': (1, 0, 1, 1),
    'rzxy': (1, 1, 0, 1), 'ryxy': (1, 1, 1, 1), 'ryxz': (2, 0, 0, 1),
    'rzxz': (2, 0, 1, 1), 'rxyz': (2, 1, 0, 1), 'rzyz': (2, 1, 1, 1)}

#Vector Operations

def vsub(u, v):
    """
    This function returns the difference between two vectors of the same dimension. Works also for flat matrices
    
    :param u: the subrahend
    :param v: the minuend
    :type u: float iterable
    :type v: float iterable
    :return: The resulting vector, vect1-vect2
    :rtype: double array
    
    """
    ret = []
    for i in xrange(len(u)):
        ret.append(u[i]-v[i])
    return ret

def vadd(*vlist):
    """
    This function sums several vectors of the same dimension. If for instance one has vectors v1,v2,v3,v4 all four having dimension n, then one can use
    vadd(v1,v2,v3,v4). This works for arbitrary number of vectors (with the same dimension), vadd(v1) is also valid. Works also for flat matrices
         
    :param vlist: the sequence without paranthesis, that determines all the vectors to be added together.
    :type vlist: a sequence of list of integers of doubles
    :return: the sum of vectors to be added
    :rtype: double or integer array

    """
    returnValue=[]
    for i in xrange(len(vlist[0])):
        a=0
        for j in xrange(len(vlist)):
            a=a+vlist[j][i]
        returnValue.append(a)
    return returnValue

def vmul(vect, s):
    """
    This function returns the vector result of multiplying each entries of a vector by a scalar. Works also for flat matrices
    
    :param vect: the vector to be multiplied with the scalar value
    :param s: the scalar value
    :type vect: double or integer iterable
    :type s: double or integer
    :return: The resulting vector s(vect1)
    :rtype: double iterable

    """
    ret=[]
    for x in vect:
        ret.append(x*s)
    return ret

def vdot(u, v):
    """
    This function returns the dot product between two vectors of the same dimension
    
    :param u: The first vector
    :param v: The second vector
    :type u: float or integer iterable
    :type v: float or integer iterable
    :return: dot-Product of u and v
    :rtype: double or integer
    """
    a=0
    for i in xrange(len(u)):
        a=a+u[i]*v[i]
    return a

def vlen(v):
    """
    This function returns the norm (length) of a vector (as a float).

    :rtype: double
    :return: euclidean norm of v
    :type  vect: float or integer iterable
    :param vect: The vector
    """
    return sqrt(vdot(v,v))


def vnorm(vect):
    """
    This function returns a normalized vector ie a unit length
    vector pointing in the same direction as the input vector.  This performs
    essentially the same function as vunit(vect) except that this function
    handles potential zero length vectors.

    :rtype: double array
    :return: normalized form of  vect 
    :type  vect: double iterable
    :param vect: The vector - in the format [x,y,z]
        (or [x,y,z,0] for affine transformations in an homogeneous space).
    """

    length = vlen(vect)

    # Keep the program from blowing up by providing an acceptable
    # value for vectors whose length may be calculated too close to zero.

    if length == 0.0:
        return len(vect)*[0.0]

    # Dividing each element by the length will result in a
    # unit normal vector.
    #ret = array('d')
    ret = []
    for x in vect:
        ret.append(x/length)
    return ret


def vdist(vect1, vect2):
    """
    This function returns the euclidean distance (the straight-line distance)
    between two vector coordinates.
    The distance between two points is the length of the vector joining them.

    :rtype: double
    :return: euclidean distance between  vect1  and  vect2  in 3D space
    :type  vect1: double iterable
    :param vect1: The vector - in the format [x,y,z]
        (or [x,y,z,0] for affine transformations in an homogeneous space).
    :type  vect2: double iterable
    :param vect2: The vector - in the format [x,y,z]
        (or [x,y,z,0] for affine transformations in an homogeneous space).
    """
    return vlen(vsub(vect1,vect2))


def vcross(vect1, vect2):
    """
    This function returns the cross product of two vectors.

    :rtype: double list
    :return: cross product M{vect1 S{times} vect2}
    :type  vect1: double list
    :param vect1: The vector - in the format [x,y,z]
        (or [x,y,z,0] for affine transformations in an homogeneous space).
    :type  vect2: double list
    :param vect2: The vector - in the format [x,y,z]
        (or [x,y,z,0] for affine transformations in an homogeneous space).
    """

    return [vect1[1] * vect2[2] - vect1[2] * vect2[1], vect1[2] * vect2[0] - vect1[0] * vect2[2], vect1[0] * vect2[1] - vect1[1] * vect2[0]]

def pseudoGrammSchmidt(v, w):
  """
  Given two linearly indeopendent vectors in 3D, this method perform the gramm-schmidt orthogonormalization of the set of vectors.
  The output is a vector normal to the first vector and belonging to the plain defined by the two vectors.
  See http://en.wikipedia.org/wiki/Gram%E2%80%93Schmidt_process.
  
  :rtype:    array of 3 doubles
  :return:   normal vector to the first input vector and belonging to the plain in which the two input vector generate
  :type  v:  array of 3 doubles
  :param v:  first input vector
  :type  w:  array of 3 doubles
  :param w:  first input vector
  """
  return vsub(w, vmul(v, vdot(w,v)/vdot(w,w)))

def isPositive(vec, point):
  """
  Given a 3D vector and a point in space we want to determine whether the point lies in the direction of the vector or not. In other words, if the vector represent an exis of 
  another coordinate system (with center still at 0,0,0) and if we project the point onto that axis, we want to determine if the projection is positive with respect to that 
  axis or not.
  
  :rtype:    a bool
  :return:   true if point lies in the direction of the vector
  :type  vec: array of 3 doubles
  :param vec: a vector in 3D
  :type  point: array of 3 doubles
  :param point: a point in 3d space
  """
  return (fabs(acos(vdot(vec,point)/(vlen(vec)*vlen(point)))) <= pi/2)

  
#Matrix Operations

def mmul(m2, m1):
    """
    .. todo::
    
        Still to comment.
    
    """
    
    return [m1[0] * m2[0]  + m1[4] * m2[1]  + m1[8]  * m2[2] ,
            m1[1] * m2[0]  + m1[5] * m2[1]  + m1[9]  * m2[2] ,
            m1[2] * m2[0]  + m1[6] * m2[1]  + m1[10] * m2[2] ,
            m1[3] * m2[0]  + m1[7] * m2[1]  + m1[11] * m2[2]  + m2[3],
            m1[0] * m2[4]  + m1[4] * m2[5]  + m1[8]  * m2[6],
            m1[1] * m2[4]  + m1[5] * m2[5]  + m1[9]  * m2[6],
            m1[2] * m2[4]  + m1[6] * m2[5]  + m1[10] * m2[6],
            m1[3] * m2[4]  + m1[7] * m2[5]  + m1[11] * m2[6]  + m2[7],
            m1[0] * m2[8]  + m1[4] * m2[9]  + m1[8]  * m2[10],
            m1[1] * m2[8]  + m1[5] * m2[9]  + m1[9]  * m2[10],
            m1[2] * m2[8]  + m1[6] * m2[9]  + m1[10] * m2[10],
            m1[3] * m2[8]  + m1[7] * m2[9]  + m1[11] * m2[10] + m2[11],
            0.0, 0.0, 0.0, 1.0]

def flatten(M):
    """
    For readability it is easier to write matrices as list of list of doubles. In most cases we do this. But for speed and efficiency,
    we it is best to have these matrices as an (flattened matrix) array. This function converts a list of list into an array.
    
    :rtype: array
    :return: an array object
    :type  M: double iterable
    :param M: Matrix to convert
    """
    #N=array('d')
    N=[]
    for i in xrange(len(M)):
        for j in xrange(len(M[0])):
            N.append(M[i][j])
    return N

def _unFlatten(M,rows,cols):
    N = []
    #N=array('d')
    for i in xrange(rows):
        row = []
        n=i*cols
        for j in xrange(cols):
            row.append(M[n+j])
        N.append(row)
    return N

def zeros(*shape):
    """
    This function returns an multidimensional zero-matrix (row-major, list of lists) or zero-vector (list of doubles). For instance: If you want to have a zero-vector of 3-dimensions you type
    zeros(3). If you want a 2x3 zero-matrix, we write zeros(2,3).

    :rtype:    list of double lists
    :return:   a matrix represented as list of lists. Each entry of the list represents a row of the matrix (if this is a nxm matrix). The representation is a row-major order.
    :type  shape:  sequence of integers (e.g. 2,3 or 2)
    :param shape:  this represent the dimensions (in integer tuples) of the output matrix (e.g. for 2x3 matrix shape is 2,2)
    """
    if len(shape) == 0:
        return 0.0
    car = shape[0]
    cdr = shape[1:]
    return [zeros(*cdr) for i in xrange(car)]

def _unitMatrix(n):
    """
    This function returns an nxn unit matrix of doubles.

    :rtype:    array of doubles
    :return:   an nxn flat unit-matrix, row-major order.
    :type  n:  integer
    :param n:  the size of the row of the unit-matrix
    """
    M=array('d')
    for i in xrange(n):
        for j in xrange(n):
            if (i==j): M.append(1.0)
            else: M.append(0.0)
    return M

def _transpose(M,rows=0,cols=0):
    """
    This function returns the transpose of a flat matrix

    :rtype:    double array
    :return:   a matrix that is the transpose of the input matrix (row-major)
    :type  M:  iterable of doubles or integers
    :param M:  the input flat matrix (row-major) that we want to transpose
    :type  rows:  integer
    :param rows:  number of rows that M has (as a row-major matrix)
    :type  cols:  integer
    :param colss: number of columns that M has (as a row-major matrix)

    """
    #ret = array('d')
    ret = []
    for i in xrange(cols):
      for j in xrange(rows):
        ret.append(M[i+j*cols])
    return ret

def _vmulv(u,v):
    """
    This function returns the matrix uv^T (where T here means transpose).

    :rtype:    array of doubles
    :return:   flat matrix uv^T (row-major)
    :type  u:  double iterable
    :param u:  the vector multiplied from left
    :type  v:  double iterable
    :param v:  the vector multiplied whose adjoint is multiplied from right
    """
    #M=array('d')
    M=[]
    for i in xrange(len(u)):
        for j in xrange(len(v)):
            M.append(u[i]*v[j])
    return M

# Warning: Unfinished! 
def _QR(M,n):
    """
    
    .. warning::
    
        **Unfinished!**
        
    QR-Decomposition of a flat singular square matrix using Householder transformations. 

    :rtype:    tuple of array of doubles
    :return:   a tuple of flat matrices first matrix is an array representing Q, second matrix represents R for the QR-decomposition of M
    :type  M:  array of doubles
    :param M:  flat square matrix (row-major) that we want to take the QR-decomposition
    :type  n:  integer
    :param n:  dimension of the square matrix M
    """
    A=M[:] #deep copy for a flat iterable. warning [:] does shallow copy for multidimensional iterables
    R=n*n*array('d',[0]) #zero matrix
    for j in xrange(n):
        m=n-j
        x=array('d')
        e=m*array('d',[0])
        e[0]=1.0
        for i in xrange(m):
            x.append(A[i])
        v = vadd(x,vmul(vlen(x),e))
        d=vlen(v) #nonzero because A is singular
        d=2/(d*d)
        P=vsub(unitMatrix(m),vmul(vmulv(v,v),d))
        #A=_mmul(P,A,n,n,n)

        B=array('d')
        #smart matrix matrix multiplication extracting the lower submatrix into B
        #i.e. : removing the first row and column of the multiplication of P and A and assigning B to it
        # see how Householder Transformation are created in QR-Decomposition!
        for i in xrange(m):
            m=i*(n-j)
            for j2 in xrange(m):
                a=0
                for k in xrange(m):
                    a=a+P[m+k]*A[k*m+j2]
                if j2==0:
                    R[j2+j+j*n]=a
                else: B.append(a)
        A=B
        #A= A
    #v is not zero because the matrix is singular
    #jocapsco: Todo .. finish this...

def _mmul(M,N,rowsM,colsM,colsN):
    """
    This is the naive matrix multiplication. There are faster matrix multiplication algorithms (like those by
    Strassen (http://en.wikipedia.org/wiki/Strassen_algorithm) or
    Coppersmith-Winograd (http://en.wikipedia.org/wiki/Coppersmith-Winograd_algorithm. But fast algorithms will make our
    code uneccessarily long and complicated and for small sized matrix (in 3D programming most matrix
    operation are limited to 3x3 matrices) the performance improvement is insignifcant.

    :rtype:    array of doubles
    :return:   a flat mxp matrix reprenting the product of M and N
    :type  M:  array of doubles
    :param M:  flat mxn matrix (row-major), that is supposed to be the left-multiplier
    :type  rowsM:  integer
    :param rowsM:  number of rows of M
    :type  colsM:  integer
    :param colsM:  number of columns of M = number of rows of N
    :type  colsN:  integer
    :param colsN:  number of columns of N
    """
    #P=array('d')
    P=[]
    for i in xrange(rowsM):
        n=i*colsM
        for j in xrange(colsN):
            a=0
            for k in xrange(colsM):
                a=a+M[n+k]*N[k*colsM+j]
            P.append(a)
    return P

    

#Quaternions

# Quaternions are of the form (x,y,z,w)
def qmul(q1, q2):
  return [q1[1]*q2[2] - q1[2]*q2[1] + q1[0]*q2[3] + q2[0]*q1[3], 
          q1[2]*q2[0] - q2[2]*q1[0] + q1[1]*q2[3] + q2[1]*q1[3], 
          q1[0]*q2[1] - q2[0]*q1[1] + q1[2]*q2[3] + q2[2]*q1[3], 
          q1[3]*q2[3] - q1[0]*q2[0] - q1[1]*q2[1] - q1[2]*q2[2]] 

def quaternionVectorTransform(q, v):
    return [q[3]*q[3]*v[0] + 2*q[1]*q[3]*v[2] - 2*q[2]*q[3]*v[1] + q[0]*q[0]*v[0] + 2*q[1]*q[0]*v[1] + 2*q[2]*q[0]*v[2] - q[2]*q[2]*v[0] - q[1]*q[1]*v[0],
            2*q[0]*q[1]*v[0] + q[1]*q[1]*v[1] + 2*q[2]*q[1]*v[2] + 2*q[3]*q[2]*v[0] - q[2]*q[2]*v[1] + q[3]*q[3]*v[1] - 2*q[0]*q[3]*v[2] - q[0]*q[0]*v[1],
            2*q[0]*q[2]*v[0] + 2*q[1]*q[2]*v[1] + q[2]*q[2]*v[2] - 2*q[3]*q[1]*v[0] - q[1]*q[1]*v[2] + 2*q[3]*q[0]*v[1] - q[0]*q[0]*v[2] + q[3]*q[3]*v[2]]

def axisAngleToQuaternion(axis, angle):
    s = sin(angle/2.0)
    qx = axis[0] * s
    qy = axis[1] * s
    qz = axis[2] * s
    qw = cos(angle/2.0)
    return (qx, qy, qz, qw)
    
def quaternionTranslationToDual(q, t):
    return [q,
            [0.5 * ( t[0] * q[3] + t[1] * q[2] - t[2] * q[1]),
             0.5 * (-t[0] * q[2] + t[1] * q[3] + t[2] * q[0]),
             0.5 * ( t[0] * q[1] - t[1] * q[0] + t[2] * q[3]),
            -0.5 * ( t[0] * q[0] + t[1] * q[1] + t[2] * q[2])]]

# todo: correct to row-major and test validity
def dualToMatrix(d):
    # Since the rotation part is a unit quaternion, we don't need to divide I think
    #length = vdot(d[0], d[0])
    x, y, z, w = d[0]
    t1, t2, t3, t0 = d[1]
    m = [[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]
        
    m[0][0] = w*w + x*x - y*y - z*z
    m[1][0] = 2.0*x*y - 2.0*w*z
    m[2][0] = 2.0*x*z + 2.0*w*y
    m[0][1] = 2.0*x*y + 2.0*w*z
    m[1][1] = w*w + y*y - x*x - z*z
    m[2][1] = 2.0*y*z - 2.0*w*x
    m[0][2] = 2.0*x*z - 2.0*w*y
    m[1][2] = 2.0*y*z + 2.0*w*x
    m[2][2] = w*w + z*z - x*x - y*y
    
    m[3][0] = -2.0*t0*x + 2.0*t1*w - 2.0*t2*z + 2.0*t3*y
    m[3][1] = -2.0*t0*y + 2.0*t1*z + 2.0*t2*w - 2.0*t3*x
    m[3][2] = -2.0*t0*z - 2.0*t1*y + 2.0*t2*x + 2.0*t3*w
    
    m[0][3] = 0.0
    m[1][3] = 0.0
    m[2][3] = 0.0
    m[3][3] = 1.0
    
    #mdiv(m, length)
    
    return m

#Note: Quaternions have to of normalized form
# Quaternions are of the form (x,y,z,w)
def quaternion2Matrix(q):
    m = [ [1.0, 0.0, 0.0],
          [0.0, 1.0, 0.0],
          [0.0, 0.0, 1.0]]  # will be a 3x3 euler rotation matrix
    m[0][0] = float(q[3]*q[3] + q[0]*q[0] - q[1]*q[1] - q[2]*q[2])
    m[0][1] = 2.0*(q[0]*q[1]-q[3]*q[2])
    m[0][2] = 2.0*(q[0]*q[2]+q[3]*q[1])

    m[1][0] = 2.0*(q[1]*q[0]+q[3]*q[2])
    m[1][1] = float(q[3]*q[3]-q[0]*q[0]+q[1]*q[1]-q[2]*q[2])
    m[1][2] = 2.0*(q[1]*q[2]-q[3]*q[0])

    m[2][0] = 2.0*(q[2]*q[0]-q[3]*q[1])
    m[2][1] = 2.0*(q[2]*q[1]+q[3]*q[0])
    m[2][2] = float(q[3]*q[3]-q[0]*q[0]-q[1]*q[1]+q[2]*q[2])
    return m

def matrix2Quaternion(m):
    q = [0.0, 0.0, 0.0, 1.0];
    t = 1.0 + m[0][0] + m[1][1] + m[2][2]
    r = 0.0;
    i = 0;

    if (t == 0.0):
      return q
    elif (t > 0):
      r = sqrt(1.0 + m[0][0] + m[1][1] + m[2][2]);
    else:
      if ((m[0][0] > m[1][1]) and (m[0][0] > m[2][2])):
        i = 1
      elif (m[1][1] > m[2][2]):
        i = 2
      else:
        i = 3
        
      r = 2 * m[i - 1][i - 1] - m[0][0] - m[1][1] - m[2][2];
  

    sgn = 1 - 2 * (i % 2)

    q[(i + 3) % 4] = r / 2.0
    q[(sgn + i + 3) % 4] = (m[2][1] - m[1][2]) / (2 * r)
    q[(2 * sgn + i + 3) % 4] = (m[0][2] - m[2][0]) / (2 * r)
    q[(3 * sgn + i + 3) % 4] = (m[1][0] - m[0][1]) / (2 * r)
    
    return q
 
def euler2Quaternion(e, axes='sxyz'):
  return matrix2Quaternion(_unFlatten(euler2matrix(e, axes),4,4))
    
def quaternionLerp(q1, q2, alpha):
    
    return vnorm([q1[0] + alpha * (q2[0] - q1[0]),
                  q1[1] + alpha * (q2[1] - q1[1]),
                  q1[2] + alpha * (q2[2] - q1[2]),
                  q1[3] + alpha * (q2[3] - q1[3])])

'''    
def quaternionSlerp2(q1, q2, alpha):
    
    dot = vdot(q1, q2)
    
    if dot > 0.1:
        return vnorm([q1[0] + alpha * (q2[0] - q1[0]),
                      q1[1] + alpha * (q2[1] - q1[1]),
                      q1[2] + alpha * (q2[2] - q1[2]),
                      q1[3] + alpha * (q2[3] - q1[3])])
                      
    dot = max(-1.0, min(dot, 1.0))
    theta0 = acos(dot)
    theta = theta0 * alpha

    q = vnorm([q2[0] - alpha * q1[0],
               q2[1] - alpha * q1[1],
               q2[2] - alpha * q1[2],
               q2[3] - alpha * q1[3]])

    return vadd(vmul(q1, cos(theta)), vmul(q, sin(theta)))
'''
            
def quaternionSlerp(q1, q2, alpha):
        
    cosHalfTheta = q1[3] * q2[3] + q1[0] * q2[0] + q1[1] * q2[1] + q1[2] * q2[2]
    
    if abs(cosHalfTheta) >= 1.0:
        return q1

    halfTheta = acos(cosHalfTheta)
    sinHalfTheta = sqrt(1.0 - cosHalfTheta * cosHalfTheta)

    if abs(sinHalfTheta) < 0.001:
        return [q1[0] * 0.5 + q2[0] * 0.5,
                q1[1] * 0.5 + q2[1] * 0.5,
                q1[2] * 0.5 + q2[2] * 0.5,
                q1[3] * 0.5 + q2[3] * 0.5]

    ratioA = sin((1 - t) * halfTheta) / sinHalfTheta;
    ratioB = sin(t * halfTheta) / sinHalfTheta; 

    return [q1[0] * ratioA + q2[0] * ratioB,
            q1[1] * ratioA + q2[1] * ratioB,
            q1[2] * ratioA + q2[2] * ratioB,
            q1[3] * ratioA + q2[3] * ratioB]

# Axis is normalized, angle is in radians
def axisAngleToEuler(x, y, z, angle):

    s = sin(angle)
    c = cos(angle)
    t = 1-c

    if (x*y*t + z*s) > 0.998:
        
        heading = 2.0 * atan2(x*sin(angle/2.0), cos(angle/2.0))
        attitude = pi/2.0
        bank = 0.0
        return heading, attitude, bank

    if (x*y*t + z*s) < -0.998:
        
        heading = -2.0*atan2(x*sin(angle/2.0),cos(angle/2.0))
        attitude = -pi/2.0
        bank = 0.0
        return heading, attitude, bank
        
    heading = atan2(y * s- x * z * t , 1.0 - (y*y+ z*z ) * t)
    attitude = asin(x * y * t + z * s)
    bank = atan2(x * s - y * z * t , 1.0 - (x*x + z*z) * t)
    return heading, attitude, bank

"""
Geometric Operations
"""

def mulmatvec3x3(m, vect):
    """
    This function returns a 3D vector which consists of the 3D input
    vector multiplied by a 3x3 matrix.
    
    :rtype:    double iterable
    :return:   3D vector
    :type  m:  double iterable
    :param m:  The matrix to multiply
    :type  vect:  double iterable
    :param vect:  The vector - in the format[x,y,z]
        (or [x,y,z,0] for affine transformations in an homogeneous space)
  
    """


    r = [0.0, 0.0, 0.0]
    r[0] = vect[0] * m[0][0] + vect[1] * m[1][0] + vect[2] * m[2][0]
    r[1] = vect[0] * m[0][1] + vect[1] * m[1][1] + vect[2] * m[2][1]
    r[2] = vect[0] * m[0][2] + vect[1] * m[1][2] + vect[2] * m[2][2]
    return r


def makeRotEulerMtx3D(rx, ry, rz):
    """
    This function returns a 3x3 euler rotation matrix based on the 3 angles
    rx, ry and rz.
    
    :rtype:    double iterable
    :return:   3x3 euler rotation matrix
    :type  rx:  float
    :param rx:  The angle of rotation (in radians) around the x-axis
    :type  ry:  float
    :param ry:  The angle of rotation (in radians) around the x-axis
    :type  rz:  float
    :param rz:  The angle of rotation (in radians) around the x-axis

    """

    SRX = sin(rx)
    SRY = sin(ry)
    SRZ = sin(rz)
    CRX = cos(rx)
    CRY = cos(ry)
    CRZ = cos(rz)

    return [[CRY * CRZ, CRY * SRZ, -SRY], [(CRZ * SRX) * SRY - CRX * SRZ, CRX * CRZ + (SRX * SRY) * SRZ, CRY * SRX], [SRX * SRZ + (CRX * CRZ) * SRY, (CRX * SRY) * SRZ
             - CRZ * SRX, CRX * CRY]]

def makeTransform(Rxyz, Txyz):
    """
    Makes a flat 4x4 homogenous transform matrix (row-major). Using xyz rotations and position
    
    :rtype:    double iterable
    :return:   4x4 roto-translation matrix
    :type  Rxyz:  double iterable
    :param Rxyz:  A list or rotations around X,Y and Z axis, in radians.
    :type  Txyz:  double iterable
    :param Txyz:  The translation vector along X,Y and Z axis.
    
    """
   
    sx = sin(Rxyz[0])
    sy = sin(Rxyz[1])
    sz = sin(Rxyz[2])
    cx = cos(Rxyz[0])
    cy = cos(Rxyz[1])
    cz = cos(Rxyz[2])

    return [cy*cz           , cy*sz           , -sy  , Txyz[0], 
            cz*sx*sy - cx*sz, cx*cz + sx*sy*sz, cy*sx, Txyz[1],
            sx*sz + cx*cz*sy, cx*sy*sz - cz*sx, cx*cy, Txyz[2],
            0.0             , 0.0             , 0.0  , 1.0]


def makeRotEulerMtx2D(theta, rotAxe):
    """
    This function returns a 3x3 euler matrix that rotates a point on
    a plane perpendicular to a specified rotational axis.

    :rtype:    double iterable
    :return:   3x3 euler rotation matrix
    :type  theta:  float
    :param theta:  The angle of rotation (in radians)
    :type  rotAxe:  string
    :param rotAxe:  The axis of rotation, which can be "X", "Y" or "Z".
 
    """

    if rotAxe == 'X':
        Rmtx = makeRotEulerMtx3D(theta, 0, 0)
    elif rotAxe == 'Y':
        Rmtx = makeRotEulerMtx3D(0, theta, 0)
    elif rotAxe == 'Z':
        Rmtx = makeRotEulerMtx3D(0, 0, theta)
    return Rmtx


def makeRotMatrix(angle, axis):
    """
    This function returns a 3x3 transformation matrix that represents a
    rotation through the specified angle around the specified axis.
    This matrix is presented in Graphics Gems (Glassner, Academic Press, 1990),
    and discussed here: http://www.gamedev.net/reference/programming/features/whyquats/

    :rtype:    double iterable
    :return:   3x3 euler matrix
    :type  angle: float
    :param angle: The angle of rotation (rad) around the specified axis
    :type  axis: double iterable
    :param axis: A 3d vector [x,y,z] defining the axis of rotation
        (this should already be normalized to avoid strange results)

    """

    a = angle
    x = axis[0]
    y = axis[1]
    z = axis[2]
    t = 1 - cos(a)
    c = cos(a)
    s = sin(a)
    M11 = (t * x) * x + c
    M12 = (t * x) * y + s * z
    M13 = (t * x) * z - s * y
    M21 = (t * x) * y - s * z
    M22 = (t * y) * y + c
    M23 = (t * y) * z + s * x
    M31 = (t * x) * z + s * y
    M32 = (t * y) * z - s * x
    M33 = (t * z) * z + c
    return [[M11, M12, M13], [M21, M22, M23], [M31, M32, M33]]

def rotMatrix2Matrix4(m):
  return [ m[0][0], m[0][1], m[0][2], 0.0,
           m[1][0], m[1][1], m[1][2], 0.0,
           m[2][0], m[2][1], m[2][2], 0.0,
           0.0, 0.0, 0.0, 1.0]
    
def makeUnit():
    
    return [1.0, 0.0, 0.0, 0.0,
            0.0, 1.0, 0.0, 0.0,
            0.0, 0.0, 1.0, 0.0,
            0.0, 0.0, 0.0, 1.0]
    
def makeTranslation(x, y, z):
    
    return [1.0, 0.0, 0.0, x,
            0.0, 1.0, 0.0, y,
            0.0, 0.0, 1.0, z,
            0.0, 0.0, 0.0, 1.0]
            
def getTranslation(m):
  """
   get the translation vector of a homogeneous 4x4 row-major transformation matrix
  """
  return [m[3], m[7], m[11]]
    
def makeRotation(axis, angle):

    c = cos(angle)
    s = sin(angle)
    t = 1 - c
    x, y, z = axis
    return [x*x * t + c,   y*x * t - z*s, x*z * t + y*s, 0.0,
            x*y * t + z*s, y*y * t + c,   y*z * t - x*s, 0.0,
            x*z * t - y*s, y*z * t + x*s, z*z * t + c,   0.0,
            0.0,           0.0,           0.0,           1.0]

def makeScale(scale):
        
    if type(scale) is float:
      scale = [scale,scale,scale]
    
    return [scale[0],   0.0, 0.0, 0.0,
            0.0, scale[1],   0.0, 0.0,
            0.0, 0.0, scale[2],   0.0,
            0.0, 0.0, 0.0, 1.0]
            
def mtransform(m, v):  
       
    return[m[0]*v[0] + m[1]*v[1] + m[2]*v[2]  + m[3],
           m[4]*v[0] + m[5]*v[1] + m[6]*v[2]  + m[7],
           m[8]*v[0] + m[9]*v[1] + m[10]*v[2] + m[11]]



def invTransform(m):
    """
    A fast way to inverse a homogenous 4x4 flat row-major transformation matix
    we use the fact that rotations are orthogonal 
    
    .. note::
    
        there shouldnt be scaling in the matrix)
    
    """
    rinv =   [m[0], m[4], m[8], 0,
              m[1], m[5], m[9], 0,
              m[2], m[6], m[10],0,
              0.0, 0.0, 0.0, 1.0]
    t =  mtransform(rinv, [-m[3],-m[7],-m[11]])
    return [m[0], m[4], m[8], t[0],
            m[1], m[5], m[9], t[1],
            m[2], m[6], m[10],t[2],
            0.0, 0.0, 0.0, 1.0]

# uses flat row-major 4x4 transformation matrices. Returned angles are in radians            
def matrix2euler(m, ):
    """
    See: http://www.lfd.uci.edu/~gohlke/code/transformations.py.html
    """
    _NEXT_AXIS = [1, 2, 0, 1]
    firstaxis, parity, repetition, frame = _AXES2TUPLE[axes.lower()]
    i = firstaxis
    j = _NEXT_AXIS[i+parity]
    k = _NEXT_AXIS[i-parity+1]

    if repetition:
        sy = sqrt(m[i+4*j]*m[i + 4*j] + m[i+4*k]*m[i+4*k])
        if sy > machine_epsilon:
            ax = atan2( m[i+4*j],  m[i+4*k])
            ay = atan2( sy,        m[i+4*i])
            az = atan2( m[j+4*i], -m[k+4*i])
        else:
            ax = math.atan2(-m[j+4*k],  m[j+4*j])
            ay = math.atan2( sy,        m[i+4*i])
            az = 0.0
    else:
        cy = math.sqrt(m[i+4*i]*m[i+4*i] + m[j+4*i]*m[j+4*i])
        if cy > machine_epsilon:
            ax = math.atan2( m[k+4*j],  m[k+4*k])
            ay = math.atan2(-m[k+4*i],  cy)
            az = math.atan2( m[j+4*i],  m[i+4*i])
        else:
            ax = math.atan2(-m[j+4*k],  m[j+4*j])
            ay = math.atan2(-m[k+4*i],  cy)
            az = 0.0
            
    if parity:
        ax, ay, az = -ax, -ay, -az
    if frame:
        ax, az = az, ax
    return [ax, ay, az]

#angles are radians!, returns flat matrix!
def euler2matrix(rotation, axes='sxyz'):
    """
    Return homogeneous rotation matrix from Euler angles and axis sequence.
    see: http://www.lfd.uci.edu/~gohlke/code/transformations.py.html
    
    """
    ai, aj, ak = rotation[0], rotation[1], rotation[2]
    firstaxis, parity, repetition, frame = _AXES2TUPLE[axes]
    i = firstaxis
    j = _NEXT_AXIS[i+parity]
    k = _NEXT_AXIS[i-parity+1]

    if frame:
        ai, ak = ak, ai
    if parity:
        ai, aj, ak = -ai, -aj, -ak

    si, sj, sk = sin(ai), sin(aj), sin(ak)
    ci, cj, ck = cos(ai), cos(aj), cos(ak)
    cc, cs = ci*ck, ci*sk
    sc, ss = si*ck, si*sk

    m = makeUnit()
    if repetition:
        m[4*i+i] = cj
        m[4*i+j] = sj*si
        m[4*i+k] = sj*ci
        m[4*j+i] = sj*sk
        m[4*j+j] = -cj*ss+cc
        m[4*j+k] = -cj*cs-sc
        m[4*k+i] = -sj*ck
        m[4*k+j] = cj*sc+cs
        m[4*k+k] = cj*cc-ss
    else:
        m[4*i+i] = cj*ck
        m[4*i+j] = sj*sc-cs
        m[4*i+k] = sj*cc+ss
        m[4*j+i] = cj*sk
        m[4*j+j] = sj*ss+cc
        m[4*j+k] = sj*cs-sc
        m[4*k+i] = -sj
        m[4*k+j] = cj*si
        m[4*k+k] = cj*ci
    return m
  
#converts a matrix (flat, homogenous row-major transformation)  to rotation, position
# where rotation is  an euler xyz rotation and position is the position in the form [x,y,z]
def makeXYZPos(m):
  position = [m[3], m[7], m[11]]
  #xyzRot = 
  pass
 
def centroid(vertsList):
    """
    This function returns the baricenter of a set of coordinate vectors
    [[x1,y1,z1],[x2,y2,z2],...,[xn,yn,zn]], returning a coordinate vector
    formatted as a double list [double X,double Y, double Z].
    This is the sum of all of the vectors divided by the number of vectors.

    :rtype:       double iterable
    :return:      the centroid of the convex hull of all the vertices in M(vertsList)
    :type  vertsList: list of double lists
    :param vertsList: each vector in the list is in the format [x,y,z]
        (or [x,y,z,0] for affine transformations in an homogeneous space).
    """

    nVerts = len(vertsList)
    xTot = 0.0
    yTot = 0.0
    zTot = 0.0
    for v in vertsList:
        xTot += v[0]
        yTot += v[1]
        zTot += v[2]
    if nVerts != 0:
        centrX = xTot / nVerts
        centrY = yTot / nVerts
        centrZ = zTot / nVerts
    else:
        print 'Warning: no verts to calc centroid'
        return 0
    return [centrX, centrY, centrZ]

    
def rotatePoint(center, vect, rotMatrix):
    """
    This function returns the 3D vector coordinates of a
    vector rotated around a specified centre point using a
    3x3 rotation matrix.

    :rtype:       Double iterable
    :return:      Rotated point
    :type  center: Double iterable
    :param center: A 3D vector - in the format[x,y,z] containing the
        coordinates of the center of rotation.
    :type vect: Double iterable
    :param vect: A 3D vector - in the format[x,y,z] containing the
        coordinates of the point to be rotated.
    :type rotMatrix: Double iterable
    :param rotMatrix: A 3x3 rotation matrix.

    """
    # subtract rotation point
    tv = vsub(vect, center)

    # rotate
    nv = mulmatvec3x3(rotMatrix, tv)

    # add the rotation point back again
    nv = vadd(nv, center)
    return nv


def scalePoint(center, vect, scale, axis=None):
    """
    This function returns the 3D vector coordinates of a
    coordinate vector scaled relative to a specified centre point using a
    scalar value.

    :rtype:       Double iterable
    :return:      Scaled point
    :type  center: Double iterable
    :param center: A 3D vector - in the format[x,y,z] containing the
        coordinates of the center of rotation.
    :type vect: Double iterable
    :param vect: A 3D vector - in the format[x,y,z] containing the
        coordinates of the point to be scaled.
    :type scale: Float
    :param scale: Scale factor.
    :type axis: String
    :param axis: An optional axis to constrain scaling ("X", "Y", "Z" or None).
        If an axis is specified then no scaling takes place along that axis.

    """

    # subtract centre point
    tv = vsub(vect, center)

    # scale
    if axis == 'X':
        nv = [tv[0], tv[1] * scale, tv[2] * scale]
    elif axis == 'Y':
        nv = [tv[0] * scale, tv[1], tv[2] * scale]
    elif axis == 'Z':
        nv = [tv[0] * scale, tv[1] * scale, tv[2]]
    else:
        nv = [tv[0] * scale, tv[1] * scale, tv[2] * scale]

    # add the centre point back again
    nv = vadd(nv, center)
    return nv


def planeNorm(vect1, vect2, vect3):
    """
    This function returns the vector of the normal to a plane, where the
    plane is defined by the vector coordinates of 3 points on that plane.
    This function calculates two direction vectors eminating from
    vect2, and calculates the normalized cross product to derive
    a vector at right angles to both direction vectors.
    This function assumes that the input coordinate vectors
    do not all lie on the same straight line.

    :rtype:       Double iterable
    :return:      The plane normal
    :type  vect1: Double iterable
    :param vect1: A 3D vector - in the format[x,y,z] containing the
        coordinates of the first point.
        (or [x,y,z,0] for affine transformations in an homogeneous space)
    :type  vect2: Double iterable
    :param vect2: A 3D vector - in the format[x,y,z] containing the
        coordinates of the first point.
        (or [x,y,z,0] for affine transformations in an homogeneous space)
    :type  vect3: Double iterable
    :param vect3: A 3D vector - in the format[x,y,z] containing the
        coordinates of the first point.
        (or [x,y,z,0] for affine transformations in an homogeneous space)
    
    """

    # Calculate two vectors from the three points

    v1 = [vect1[0] - vect2[0], vect1[1] - vect2[1], vect1[2] - vect2[2]]
    v2 = [vect2[0] - vect3[0], vect2[1] - vect3[1], vect2[2] - vect3[2]]

    # Take the cross product

    normal = [v1[1] * v2[2] - v1[2] * v2[1], v1[2] * v2[0] - v1[0] * v2[2], v1[0] * v2[1] - v1[1] * v2[0]]
    
    return normal
    '''
    # Normalize
    length = sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2])
    if length < machine_epsilon: #should not happen but it does for badly formed meshes
      return [0.0,0.0,0.0]
    else:
      return [normal[0] / length, normal[1] / length, normal[2] / length]
    '''

def focalToFov(dimension, focal):
    if focal == 0:
        return 0
    else:
        return 2 * atan2(dimension * 0.5, focal)


def fovToFocal(dimension, fov):
    return dimension / (2 * tan(fov / 2))


def in2pts(point1, point2, t):
    """
    This function returns a vector that lies on the directed line between points, given
    a parameter t. The paraemeter t is between 0 and 1 and it parametrizes our directed line
    between point1 and point2
    
    :rtype:       Double iterable
    :return:      A vector that lies on the directed line between points
    :type  point1: Double iterable
    :param point1: A 3D vector - in the format[x,y,z] containing the
        coordinates of the first point (i.e. starting point) of a directed line.
    :type  point2: Double iterable
    :param point2: A 3D vector - in the format[x,y,z] containing the
        coordinates of the first point (i.e. starting point) of a directed line.
    :type  t: Float
    :param t: A real number between 0 and 1, that linearly parametrizes
        the directed line between point1 and point2. In other words, when t is 0 the
        return value for this function is point1 and when t is 1 the return value
        for this function is point2.

    """

    return vadd(vmul(point1, 1 - t), vmul(point2, t))

def vectorsToRotMatrix(v1,v2):
    """
    Given two points v1 and v2 in 3D space. Suppose furthermore that v2 is generated by v1 
    multiplied by a rotation matrix over the origin. We will determine this rotation matrix using this method
    
    :rtype:    float list
    :return:   Rotation matrix (unhomegenized) in 3D space that rotated v1 to v2 (center at origin)
    :type  v1: float list
    :param v1: original point in 3D 
    :type  v2: float list
    :param v2: the point for which v1 has rotated to
    """
    normal = vcross(v1,v2)
    normal = vnorm(normal)
    angle = acos(vdot(v1,v2)/(vlen(v1)*vlen(v2)))
    q = axisAngleToQuaternion(normal, angle)
    return quaternion2Matrix(q)

def randomPointFromNormal(v):
    """
    Suppose you have a vector v. Then this vector defines a plane that passes the origin and for which the vector v is a normal to it
    This method gets a random point in 3D that lies on this plane. Beware: v should not be of length 0!
    """
    x = random()
    y = random()
    z = random()
    if (v[0] != 0):
      x = -(y*v[1]+z*v[2])/v[0]
    elif (v[1] !=0):
      y = -(x*v[0]+z*v[2])/v[1]
    else:
      z = -(x*v[0]+y*v[1])/v[2]
    return [x,y,z]
    
def convexQuadrilateralArea(v1,v2,v3,v4):
    """
    This function returns the area of a Quadrilateral. See U{http://mathworld.wolfram.com/Quadrilateral.html}.

    :rtype:    float
    :return:   The area of a Quadrilateral determined by v1 to v4 (clockwise or counterclockwise order)
    :type  v1: float list
    :param v1: first vertex of a parallelogram - in the format [x,y,z]
    :type  v2: float list
    :param v2: second vertex of a parallelogram - in the format [x,y,z]
    :type  v3: float list
    :param v3: third vertex of a parallelogram - in the format [x,y,z]
    :type  v4: float list
    :param v4: fourth vertex of a parallelogram - in the format [x,y,z]
    """
    #a=vdist(v2,v1)
    #b=vdist(v3,v2)
    #c=vdist(v4,v3)
    #d=vdist(v1,v4)
    p=vdist(v2,v4)
    q=vdist(v1,v3)
    pq = vdot(vsub(v3,v1),vsub(v4,v2))
    return 0.5*sqrt(p*p*q*q - pq*pq)
    #return sqrt(4*p*p*q*q - pow((b*b+d*d-a*a-c*c),2))/4
        
def calcBBox(verts, indices=None):
    
    if indices:
        bbox = [verts[indices[0]].co[:], verts[indices[0]].co[:]]
        for i in indices:
            if verts[i].co[0] < bbox[0][0]: #minX
                bbox[0][0] = verts[i].co[0]
            if verts[i].co[0] > bbox[1][0]: #maxX
                bbox[1][0] = verts[i].co[0]
            if verts[i].co[1] < bbox[0][1]: #minY
                bbox[0][1] = verts[i].co[1]
            if verts[i].co[1] > bbox[1][1]: #maxY
                bbox[1][1] = verts[i].co[1]
            if verts[i].co[2] < bbox[0][2]: #minZ
                bbox[0][2] = verts[i].co[2]
            if verts[i].co[2] > bbox[1][2]: #maxX
                bbox[1][2] = verts[i].co[2]
        return bbox
    else:
        bbox =  [verts[0].co[:],verts[0].co[:]]
        for v in verts:
            if v.co[0] < bbox[0][0]: #minX
                bbox[0][0] = v.co[0]
            if v.co[0] > bbox[1][0]: #maxX
                bbox[1][0] = v.co[0]
            if v.co[1] < bbox[0][1]: #minY
                bbox[0][1] = v.co[1]
            if v.co[1] > bbox[1][1]: #maxY
                bbox[1][1] = v.co[1]
            if v.co[2] < bbox[0][2]: #minZ
                bbox[0][2] = v.co[2]
            if v.co[2] > bbox[1][2]: #maxX
                bbox[1][2] = v.co[2]
        return bbox


def quadPrismPseudoVol(p1,p2,p3,p4):
  """
  This function returns the unit "volume" for a prism with rectangular base, as described in  
  U{http://www-ljk.imag.fr/membres/Stefanie.Hahmann/PUBLICATIONS/RHC09/}.
  
  
  :return:   the value of the unit "volume"
  :type  pi:  list of doubles
  :param pi:  coordinate values of the base of the rectangular prism
  """
  M = [4.0, 2.0, 2.0, 1.0, 2.0, 4.0, 1.0, 2.0, 2.0, 1.0, 4.0, 2.0, 1.0, 2.0, 2.0, 4.0]
  M = vmul(M, 1.0/36)
  z = [p1.co[0],p2,p3[0],p4[0]]
  k = [(p2[0]-p1[0])*(p3[1]-p1[1]) - (p2[1]-p1[1])*(p3[0]-p1[0]),
       (p2[0]-p1[0])*(p4[1]-p2[1]) - (p2[1]-p1[1])*(p4[0]-p2[0]),
       (p4[0]-p3[0])*(p3[1]-p1[1]) - (p4[1]-p3[1])*(p3[0]-p1[0]),
       (p4[0]-p3[0])*(p4[1]-p2[1]) - (p4[1]-p3[1])*(p4[0]-p2[0])]
  return _vmulv(z,_mmul(M, k, 4,4,1))
  

#Various Functions

    
def bump(x, width=1.0):
    """
    This is the bump function (see U<Wikipedia - Bump Function>{http://en.wikipedia.org/wiki/Bump_function}). Height is always 1, if we 
    need higher bump we just scale it

    :rtype:    double
    :return:   the bump function scaked to have height of 1 
    :type  x:  double
    :param x:  the value of the function at x (>= 0)
    :type  width:  double
    :param width:  radius of the bump
    """
    if (x < width):
      return exp(-width*width/(width*width - x*x) + 1.0)
    else:
      return 0.0
      
def sign(x):
    """
    The sign function

    :rtype:    integer
    :return:   the sign of x
    :type  x:  double or integer
    :param x:  any arbitrary real number
    """
    if (x<0): return -1
    elif (x>0): return 1
    else: return 0

# u and m must be float 0<=m<=1
# returns : sn,cn,dn,phi
#TODO : Add reference: Louis V. King; Hofsommer; Salzer (after reading it yourself :P)
#pg. 9 eq'n 35 of Louis V. King for m within 1.0e-9 range
def jacobianEllipticFunction(u,m):
    """
    This function returns a triple consisting of the Jacobian elliptic functions, namely the
    Jacobian sine (sn), Jacobian cosine (cn), Jacobian *TODO.. dn*, angle (in radians)

    :rtype:    Double iterable
    :return:   A triple consisting of the Jacobian elliptic functions
    :type  u:  float
    :param u:  A 3D vector - in the format[x,y,z] containing the
        coordinates of the first point (i.e. starting point) of a directed line.
    :type  m:  float
    :param m:  A value between 0 and 1 which represent the modulus of the Jacobian elliptic function

    """
    if (m< 0) or (m >1):
        print "Coefficient for Elliptic Integral should be between 1 and 0"
        return  #error-code!
    a=[0]*9
    c=[0]*9
    if m < 1.0e-9:
        t = sin(u)
        b = cos(u)
        ai = 0.25*m*(u-t*b)
        sn = t-ai*b
        cn = b+ai*t
        ph = u-ai
        dn = 1.0 - 0.5*m*t*t
        return sn,cn,dn,ph
    if m>=1.0 - 1.0e-9:
        ai = 0.25*(1.0-m)
        b = math.cosh(u)
        t= math.tanh(u)
        phi = 1.0/b
        twon = b*math.sinh(u)
        sn = t+a*(twon-u)/(b*b)
        ph = 2.0*math.atan(math.exp(u))-math.pi/2+ai*(twon-u)/b
        ai=ai*t*phi
        cn = phi - ai*(twon-u)
        dn=phi+ai*(twon+u)
        return sn,cn,dn,ph

    a[0] = 1.0;
    b=math.sqrt(1.0-m)
    c[0]=math.sqrt(m)
    twon=1.0
    i=0

    while fabs(c[i]/a[i])>machine_epsilon:
        if i>7:
            print "Overflow in the calculation of Jacobian elliptic functions"
            break
        ai = a[i]
        i=i+1
        c[i]=0.5*(ai-b)
        t=sqrt(ai*b)
        a[i]=0.5*(ai+b)
        b=t
        twon=twon*2.0

    phi=twon*a[i]*u
    while i!=0:
        t=c[i]*sin(phi)/a[i]
        b=phi
        phi=(math.asin(t)+phi)/2.0
        i=i-1
    return sin(phi),cn,cn/cos(phi-b),phi


def newton_raphson(f, f_diff, value, start, iterations=4, epsilon = 1e-4):
    """ 
    Returns a root of the polynomial, with a starting value.
    To get both roots in a quadratic, try using with n = 1 and n = -1.
    """

    x = start - (float(f(start)-value) / f_diff(start))

    x = start
    counter = 0

    while fabs(x-start)>epsilon or counter<iterations:
        x = x - (float(f(x)-value) / f_diff(x))
        counter += 1
        return x

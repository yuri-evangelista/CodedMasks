import numpy as np
from sympy import *

def next_prime(n):
    p=n
    while (not isprime(p)):
        p+=1
    return(p)

def ura_mura(p):
    #Checking if p is prime
    if isprime(p):
        #Preparing arrays
        A = np.zeros(p, dtype=int)
        R = np.arange(p, dtype=np.int64)**2 % p
        A[R] = 1

        #Check if URA or MURA can be generated
        URA = ( round( (p-3)/4 ) - ((p-3)/4)) == 0
        MURA = ( round( (p-1)/4 ) - ((p-1)/4)) == 0

        if URA:
            print("Generating URA array")
        elif MURA:
            print("Generating MURA array")
            A[0] = 0
        return A
    else:
        raise TypeError("p must be prime")

def bura(p, modified=False):
    #Checking if p is prime
    if not isprime(p):
        raise TypeError("Number of array elements must be prime")

    #Check if BURA can be generated
    x = np.sqrt( (p-1)/4)
    if not ( (int(x) == x) & ( (x % 2) == 1) ) :
        raise TypeError("p does not fulfill the requirement of p = 4x^2 + 1 with x odd")
        
    #Preparing arrays
    A = np.zeros(p, dtype=np.int64)
    R = np.arange(p, dtype=np.int64)**4 % p
    A[R] = 1

    if modified:
        A[0] = 0
        print("Generating M-BURA array")
    else:
        print("Generating BURA array")
    
    return A

def bura33(p):
    """
    Generates a biquadratic URA with OF ~0.33.
   
    ********************************************************
    The resulting code does not seem "perfect"!!!!!
    ********************************************************
    
    From Baumert L. D. 1971, Lecture Notes in Mathematics No. 182, Cyclic Difference Sets
    Theorem 5.18 (iii)

    ----------------------------------------------------------------------------------------------
    The biquadratic residues of primes v = 4x^2 + 9, x odd, form a difference set with ~0.33 OF
    Example primes: 13, 109, 1453, 3373, 3853, 4909, 6733

    """
    #Checking if p is prime
    if not isprime(p):
        raise TypeError("Number of array elements must be prime")

    x = np.sqrt((p - 9)/4.0)
    if not ( (int(x) == x) & ( (x % 2) == 1) ) :
        raise TypeError("p does not fulfill the requirement of p = 4x^2 + 9 with x odd")
 
    #Preparing arrays
    A = np.zeros(p, dtype=np.int64)
    R = np.arange(p, dtype=np.int64)**4 % p
    A[R] = 1
    #A[0] = 0 #not sure about that

    return A

def cura(p):
    """
    Generates a cubic residues URA with OF ~0.33.
   
    ********************************************************
    The resulting code does not seem "perfect"!!!!!
    ********************************************************
    
    From Lehmer, D. H. 1962, "Machine Proof of a Theorem on Cubic Residues"
    ----------------------------------------------------------------------------------------------

    """
    #Checking if p is prime
    if not isprime(p):
        raise TypeError("Number of array elements must be prime")
    
    x = (p - 1)/6.0

    if not ( (int(x) == x) & ( (x % 2) == 1) ) :
        raise TypeError("p does not fulfill the requirement of p = 6x + 1 with x odd")
 
    #Preparing arrays
    A = np.zeros(p, dtype=np.int64)
    R = np.arange(p, dtype=np.int64)**3 % p
    A[R] = 1
    
    return A
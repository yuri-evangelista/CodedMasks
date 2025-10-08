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


def is_cyclic_difference_set(s, v, k, lambda_val):
    """
    Tests if a given set is a cyclic (v, k, lambda) difference set.

    Args:
        s: A list or set of integers representing the candidate set D.
        v: The modulus (size of the cyclic group).
        k: The expected size of the set D.
        lambda_val: The expected number of times each non-zero residue appears as a difference.

    Returns:
        A tuple (bool, str) indicating whether the set is a difference set and a message.
    """

    # 1. Check if the set has the correct number of distinct elements modulo v
    unique_elements = set(x % v for x in s)
    if len(unique_elements) != k:
        return False, f"Error: Set does not contain exactly {k} distinct elements modulo {v}. Found {len(unique_elements)}."

    # 2. Check the fundamental difference set equation: k*(k-1) = lambda*(v-1)
    if k * (k - 1) != lambda_val * (v - 1):
        return False, f"Error: Parameters do not satisfy k*(k-1) = lambda*(v-1). {k}*{k-1} = {k*(k-1)}, but {lambda_val}*{v-1} = {lambda_val*(v-1)}."

    # 3. Compute all differences (di - dj) mod v
    differences = {}
    s_list = list(unique_elements) # Work with a list of unique, mod-v elements
    n = len(s_list)

    for i in range(n):
        for j in range(n):
            if i != j:
                diff = (s_list[i] - s_list[j]) % v
                differences[diff] = differences.get(diff, 0) + 1

    # 4. Check if each non-zero residue appears exactly lambda_val times
    for r in range(1, v):  # Iterate through non-zero residues
        if r not in differences or differences[r] != lambda_val:
            return False, f"Error: Residue {r} appears {differences.get(r, 0)} times, expected {lambda_val}."

    return True, f"Success: The set {s} is a cyclic ({v}, {k}, {lambda_val}) difference set."

"""
Example usage:
p = 131
k = 66
lam = (k * (k-1))/(p-1)

v = np.arange(p, dtype=np.int64)**2 % p

is_cyclic, result_set = is_cyclic_difference_set(v, p, k, lam)
print(result_set)
"""
def reshape_ext_diagonal(array, nx, ny):
    out = np.zeros( (nx, ny), dtype=array.dtype)
    for i in range (array.size):
        out[i % nx, i % ny] = array[i]
    return out
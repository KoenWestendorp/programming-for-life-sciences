# By Koen Westendorp (s3950808)
# 2022-03-21
# Programming for Life Sciences Assignment 1

# Exercise 2: Sieve of Erastosthenes

def crude_sieve(maximum: int) -> list[int]:
    r = range(2, maximum + 1) # Create inclusive range [2, maximum]
    l = list(r)
    
    # Replace all prime multiples with 0. 
    # For every number in `r`, remove all multiples of that value, except the
    # value itself.
    for n in r:
        for i, m in enumerate(l):
            # Prevent removal of the number itself, and ignore values already
            # set to 0. 
            # If m is a multple of n, its remainder is 0, and it is set to 0 in
            # `r`.
            if not (m == 0 or m == n) and m % n == 0:
                l[i] = 0
                
    return list(l)

def remove_zeroes(number_list: list[int]) -> list[int]:
    return [n for n in number_list if n != 0]

def sieve(maximum: int) -> list[int]:
    crude = crude_sieve(maximum)
    return remove_zeroes(crude)

# Exercise 3

## 3.1

#from ../../assignment_1/prime/prime import

from math import isqrt

def prime_storing(maximum: int) -> list[int]:
    found_primes = [2]
    for n in range(2, maximum + 1):
        # First, check whether `n` is divisible by any known primes. 
        # If so, reject it as a prime number.
        is_prime = True
        for found_prime in [p for p in found_primes if p < isqrt(maximum)]:
            if n % found_prime == 0:
                # `n` is divisible by a found_prime, so it cannot be prime.
                # Set `is_prime` to Fales and break out of the inner for loop.
                is_prime = False
                break

        if is_prime:
            found_primes.append(n)

    return found_primes

if __name__ == "__main__":
    n = 18
    print(f"""
Assignment 2: Prime 2
---------------------

Exercise 2
==========
For n = {n},
crude_sieve(n) 
    -> {crude_sieve(n)}
remove_zeroes([0, 0, 1, 3, 0]) 
    -> {remove_zeroes([0, 0, 1, 3, 0])}
sieve(n) 
    -> {sieve(n)}

Exercise 3
==========
prime_storing(100)
    -> {prime_storing(100)}
""")

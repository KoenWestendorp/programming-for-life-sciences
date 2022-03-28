# By Koen Westendorp (s3950808)
# 2022-03-21
# Programming for Life Sciences Assignment 1

from math import isqrt

# Prime: Exercise 1

## 1.1

def divisible(divident: int, divider: int) -> bool:
    '''
    Returns True iff the divident is neatly divided by the divider, and
    False if this is not the case.

    Parameters
    ----------
    divident : int
        The number to be divided.
    divider : int
        The number dividing the divident.

    Returns
    -------
    bool
    '''

    return divident % divider == 0

## 1.2 

def does_it_divide(n: int):
    for i in range(2, n):
        if divisible(n, i):
            print(i, "divides", n)
        else:
            print(i, "does not divide", n)

# for i in range(1, 20):
#    does_it_divide(i)

## 1.3

def is_prime(target: int) -> bool:
    '''
    Returns True if `target` is prime, otherwise False.

    Iterates over numbers from 2 to `target` using a for loop, checking whether
    `target` is divisible by that number. If the number is divisible, the
    number is not prime and False is returned as soon as such an instance is
    encountered. If divisibility is not encountered, `target` is prime, and
    True is returned.

    Parameters
    ----------
    target : int
        Number to be checked for prime.

    Returns
    -------
    bool
    '''

    for i in range(2, target):
        if divisible(target, i):
            return False

    return True

## 1.4

def is_prime_while(target: int) -> bool:
    '''
    Returns True if `target` is prime, otherwise False.

    Iterates over numbers from 2 to `target` using a while loop, checking
    whether `target` is divisible by that number. If the number is divisible,
    the number is not prime and False is returned as soon as such an instance
    is encountered. If divisibility is not encountered, `target` is prime, and
    True is returned.

    Parameters
    ----------
    target : int
        Number to be checked for prime.

    Returns
    -------
    bool
    '''

    i = 2
    while i < target:
        if divisible(target, i):
            return False

        i += 1

    return True

def pretty_string_prime(n: int, prime_function=is_prime) -> str:
    '''
    Returns a string describing whether `n` is prime.

    Parameters
    ----------
    n : int
        Number for which to check whether it is prime.
    prime_function : (int) -> bool (default = `is_prime`)
        The function used to check whether `n` is prime.

    Returns
    -------
    str
    '''

    if prime_function(n):
        return f"{n} is prime."
    else:
        return f"{n} is not prime."

# print(pretty_string_prime(1999, prime_function=is_prime_while))

## 1.5

def primes_between(start: int, end: int):
    '''
    Iterates over all numbers between `start` and `end`, and prints the number
    iff it is prime.

    Parameters
    ----------
    start : int
    end : int
    '''

    for n in range(start, end):
        if is_prime(n):
            print(n)

# primes_between(1000, 3000)

## 1.6

def is_prime_optimized(target: int) -> bool:
    '''
    Returns True if `target` is prime, otherwise False.

    Iterates over numbers from 2 to the square root of `target` using a for
    loop, checking whether `target` is divisible by that number. If the number
    is divisible, the number is not prime and False is returned as soon as such
    an instance is encountered. If divisibility is not encountered, `target` is
    prime, and True is returned.

    Parameters
    ----------
    target : int
        Number to be checked for prime.

    Returns
    -------
    bool
    '''

    # Experiment: use `int(math.sqrt(target))`, or `math.isqrt(target)`?
    #
    # Since Python 3.8, an integer square root function `math.isqrt()` is
    # available. I wondered wether this would lend superior performance, for
    # this square root of an integer is taken numerous times. My initial
    # implementation used `math.sqrt()` on the integer `target` value, which
    # returns a float, which is then coerced to an int using `int()`. Perhaps
    # the integer square root function would be more efficient than this
    # initial implementation. These are my results.
    # 
    # Using Hyperfine, a benchmarking tool, I ran 10 benchmarks on each
    # implementation by running `largest_prime(1_000_000_000_000_000)`,
    # yielding the following overviews:
    # 
    # `int(math.sqrt(target))` (the initial implemetation)
    #   Time (mean ± σ):      3.565 s ±  0.062 s    [User: 3.528 s, System: 0.021 s]
    #   Range (min … max):    3.516 s …  3.708 s    10 runs
    # 
    # `math.isqrt(target)`
    #   Time (mean ± σ):      3.548 s ±  0.013 s    [User: 3.509 s, System: 0.022 s]
    #   Range (min … max):    3.525 s …  3.568 s    10 runs
    #
    # The difference is very small, but the `isqrt()` function does seem to
    # give a slight edge, and a more tight standard deviation, and thus a
    # possibly more consistent runtime. However, a significant difference might
    # only present itself at much higher workloads. For the slight apparent
    # advantage and especially because of the better fit for the `isqrt()`
    # function in this case, I have opted to use the integer square root
    # function.
    for i in range(2, isqrt(target)):
        # Rather than using the more semantically meaningful divisible(n, m)
        # function, we can just use its internals here, inlined. 
        if target % i == 0:
            return False

    return True

# print(pretty_string_prime(2963, prime_function=is_prime_optimized))

## 1.7

def largest_prime(end: int) -> int:
    '''
    Returns the largest possible prime that is small than `end`.

    In order to find the largest prime, it is much smarter to start looking
    from the high end. Otherwise all calculations that were necessary for
    finding the primes below the largest would have been wasted.

    For values smaller than or equal to 2, the function will return the value
    of `end`.

    Parameters
    ----------
    end : int
        Maximal value under which the largest prime is looked for. Must be
        greater than 0.

    Returns
    -------
    int
    '''

    # If end is even, it can never be prime. In that case, n = end - 1, to make
    # it an odd number, and the iteration will run over all odd numbers below
    # end until it finds a prime. 
    if end % 2 == 0:
        # end is even, decrement to an odd number
        n = end - 1
    else:
        # end is odd
        n = end

    # Check all odd numbers from n and below for whether they are prime. If n
    # is prime, return n. Otherwise decrement n by 2 to the next odd number and
    # repeat.
    while n >= 2:
        # We already have gotten rid of half of the numbers to be checked: the
        # even numbers. However, we can optimize more! We can seive out any
        # number ending on 5 too, as well as on 0, because these are divisible
        # by 5.

        if is_prime_optimized(n):
            # Note that using is_prime_optimized(n) here, rather than
            # is_prime(n), produces great time gains!
            #
            # For largest_prime(int(10e7)) running `% time python3 prime.py`
            # returns:
            #
            #           is_prime: 14.64s user 0.06s system 99% cpu 14.799 total
            # is_prime_optimized: 0.03s user 0.01s system 87% cpu 0.044 total
            #
            # In fact, the larger n is, the greater the performance gain the
            # optimized function brings, because the difference between a
            # number and its square root (as used in the is_prime_optimized
            # function) quickly grows for an increasing number. 
            # In other words, the larger n is, the larger the gain the
            # optimized function brings because of its use of the square root
            # limit.
            return n
        
        n -= 2 # Decrementing 2 to the next odd number.

    # Return the value of `end` if no prime has been found. This will only
    # occur for an `end` <= 2.
    return end

if __name__ == "__main__":
    print(largest_prime(1_000_000_000_000_000))


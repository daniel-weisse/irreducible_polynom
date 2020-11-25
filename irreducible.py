"""
Generate a list of irreducible polynoms over GF(2^n) by testing all possibles
"""
from pyGF2 import gf2_mul, gf2_div, gf2_add
from numpy import array, append, binary_repr, flip, uint8, pad
from sympy import mobius, divisors
from datetime import datetime

from multiprocessing import Pool, Lock, Value, cpu_count, current_process

def time():
    """
    Return the current time
    """
    return datetime.now().strftime('%c')


def gf2_gcd(a, b):
    """
    Performs the euclidean algorithm to find the greatest common divider of polynoms a and b
    
    Parameters
    ----------
    a : ndarray (uint8 or bool)
        Input polynomial's coefficients.
    b : ndarray (uint8 or bool)
        Input polynomial's coefficients.
    
    Returns
    -------
    out : ndarray of uint8
        greatest common divider polynomial's coefficients.

    """
    if not b.any():
        return a
    else:
        _, a = gf2_div(a, b)
        return gf2_gcd(b, a)


def main(n: int, n_factors: list, cpu_limit=None):
    """
    main function for handling multiprocessing pool

    Parameters
    ----------
    n : int
        degree of polynoms to find irreducibles for
    n_factors : list of integers
        prime factors of n used to validate irreducible polynoms
    cpu_limit : int
        Maximum number of cpu cores to use for the program, None to use all avaiable

    """
    def init(l, t):
        global lock
        global total
        lock = l
        total = t
    l = Lock()
    t = Value('i', 0)

    cpu = cpu_limit or cpu_count()
    pool = Pool(cpu, initializer=init, initargs=(l, t))
    
    args = []

    gcd_polys = []
    for val in n_factors:
        x = [0, 1]
        for _ in range(n // val):
            x = gf2_mul(x, x)
        x = gf2_add(x, array([0, 1], dtype=uint8))
        gcd_polys.append(x)

    target_number = sum([mobius(d)*pow(2, n//d) for d in divisors(n)])//n

    for i in range(cpu):
        args.append(
            (
                int(2**(n - 1) * i/cpu),
                int(2**(n - 1) * (i + 1)/cpu),
                n,
                gcd_polys,
                f'irreducibles_n-{n}.txt',
                target_number
            )
        )

    print(f'{time()} Starting polynom generation for GF(2^{n}) using {cpu} cores, looking for {target_number} polynomials')

    f = open(f'irreducibles_n-{n}.txt', 'w')
    f.write(f'{n} {target_number}\n')
    f.close()

    start = datetime.now()
    total_values = pool.starmap(find_polys, args)
    total_time = datetime.now() - start

    print(f'{time()} Found {max(total_values)} irreducible polynoms over GF(2^{n}) in {total_time}')

    #f = open(f'irreducibles_n-{n}.txt', 'w+')
    #f.write(f'{n} {sum(lengths)}\n')
    #for i in range(len(irreducibles)):
    #    for j in range(len(irreducibles[i])):
    #        out = ' '.join(map(str, irreducibles[i][j]))
    #        f.write(f'{out}\n')
    #f.close()


def find_polys(start, end, n, gcd_polys: list, filename: str, target_number: int):
    """
    finds all irreducible polynomials in a given range, polynomials are directly written to file
    Parameters
    ----------
    start : int
            starting point for polynom search
    end :   int
            end point for polynom search
    n :     int
            degree of polynoms to search
    gcd_polys : list of ndarray (uint8)
            all polynoms x^2^(n/m) - x with m being a prime factor of n, used to verify irreducible polynoms

    """

    for i in range(start, end):
        #if i == (int(end * 1/4) + start):
        #    print(f'{time()} {current_process().name} 25% done', flush=True)
        #elif i == (int(end * 2/4) + start):
        #    print(f'{time()} {current_process().name} 50% done', flush=True)
        #elif i == (int(end * 3/4) + start):
        #    print(f'{time()} {current_process().name} 75% done', flush=True)
        s = array(list(binary_repr(i)), dtype=uint8)
        s = pad(
            array=s,
            pad_width=((n-1 - len(s)), 0),
            mode='constant',
            constant_values=0,
        )
        poly = append(1, s)
        poly = append(poly, 1)
        
        T_x = flip(poly)
        p_x = array([0, 1], dtype=uint8)

        for _ in range(n):
            p_x = gf2_mul(p_x, p_x)
            _, p_x = gf2_div(p_x, T_x)
        
        # if result only has 2 values there is a good chance T_x was irreducible
        if len(p_x) == 2:

            # first criteria: p_x = x
            if (p_x == [0, 1]).all():
                irr_check = True

                # second criteria: gcd(x^2^(n/m) - x, T_x) == 1, with m := prime factors of n
                for gcd_poly in gcd_polys:
                    gcd = gf2_gcd(gcd_poly, T_x)
                    if len(gcd) > 1:
                        irr_check = False
                        break
                    elif gcd[0] != 1:
                        irr_check = False
                        break

                # if both criterias are fullfilled the polynom is irreducible
                if irr_check:
                    lock.acquire()
                    total.value += 1
                    if total.value == int(target_number * 1/10):
                        print(f'{time()} 10% done')
                    elif total.value == int(target_number * 2/10):
                        print(f'{time()} 20% done')
                    elif total.value == int(target_number * 3/10):
                        print(f'{time()} 30% done')
                    elif total.value == int(target_number * 4/10):
                        print(f'{time()} 40% done')
                    elif total.value == int(target_number * 5/10):
                        print(f'{time()} 50% done')
                    elif total.value == int(target_number * 6/10):
                        print(f'{time()} 60% done')
                    elif total.value == int(target_number * 7/10):
                        print(f'{time()} 70% done')
                    elif total.value == int(target_number * 8/10):
                        print(f'{time()} 80% done')
                    elif total.value == int(target_number * 9/10):
                        print(f'{time()} 90% done')
                    elif total.value == int(target_number * 19/20):
                        print(f'{time()} 95% done')
                    f = open(filename, 'a+')
                    f.write(' '.join(map(str, T_x)))
                    f.write('\n')
                    f.close()
                    #irreducibles.append(T_x)
                    lock.release()
    return total.value

if __name__ == '__main__':
    main(16, [2])

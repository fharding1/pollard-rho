import random
import sympy
import sys
import csv

def random_dlp_zp_prime_subgroup(d):
    p = sympy.randprime(10**8,10**16)
    b = max(sympy.factorint(p-1, multiple=True)) # largest prime factor of p-1
    a = int((p-1)/b)
    if b < 10**(d-1) or b >= 10**d:
        return
    g = pow(sympy.primitive_root(p), a, p)
    x = random.randint(0,b-1)
    h = g^x
    return (p,b,g,h,x)

if __name__ == "__main__":
    d = int(sys.argv[1])
    n = int(sys.argv[2])

    with open(sys.argv[3], mode='w') as csv_file:
        csv_writer = csv.writer(csv_file,delimiter=',',quotechar='"',quoting=csv.QUOTE_MINIMAL)

        for i in range(n):
            dlp = random_dlp_zp_prime_subgroup(d)
            while dlp is None:
                dlp = random_dlp_zp_prime_subgroup(d)

            csv_writer.writerow(dlp)

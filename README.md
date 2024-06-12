# hardin

This repository contains a Python program that counts the number of rectangular
Hardinian arrays. These are families of n x k arrays, with n >= k, that satisfy
the rules laid out in the OEIS entries [A253026](https://oeis.org/A253026) and
[A253223](https://oeis.org/A253223).

The implementation relies on the [SymPy](https://www.sympy.org/en/index.html)
package, which you can probably install using `pip` if you don't already have
it:

```
$ pip install sympy
```

Here is a brief demo:

```python
>>> from sympy import var
>>> n, x = var("n x")
>>> import hardin
>>> hardin.rectGF(x, 2, 1)
-x**3/(x - 1)**3
>>> hardin.rectGF(x, 2, 2)
-x**3*(x**2 + 6*x + 1)/(x - 1)**3
>>> hardin.rectGF(x, 2, 3)
-x**3*(2*x**4 + 13*x**3 + 48*x**2 + 16*x + 1)/(x - 1)**3

# compare the below with conjectures at https://oeis.org/A253223
>>> hardin.rectPoly(n, 2, 1) 
(n - 2)*(n - 1)/2
>>> hardin.rectPoly(n, 2, 2)
(2*n - 5)**2
>>> hardin.rectPoly(n, 2, 3)
40*n**2 - 279*n + 497
>>> hardin.rectGF(x, 3, 8)
2*x**8*(1019797*x**10 + 15878036*x**9 + 112695723*x**8 + 479179922*x**7 + 1345053211*x**6 + 2574720204*x**5 + 3309668803*x**4 + 2622816901*x**3 + 852822279*x**2 + 92194141*x + 2457863)/(x - 1)**4
>>> hardin.rectPoly(n, 3, 8)
11408506880*n**3/3 - 116708851712*n**2 + 3632793111640*n/3 - 4248252440142
```

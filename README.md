## Introduction

There is a multiplicity of applications requiring the computation of averages of various functions over the unitary groups in dimension $d$ or the behaviour of such functions in the limit of large $d$ [1-4].  The computation of this type of integral is considerably simplified using the so-called *Weingarten functions*, a nomenclature introduced in Ref.[5].

Weingarten functions depend only on a class of symmetric group $S_p$ and on the dimension $d$ of the unitaries that are averaged.  A convenient closed form expression has been given in Ref. [5]:

$$
\int U_{i_1j_1}\ldots U_{i_pj_p} \left(U_{i^\prime_1j^\prime_1}\ldots U_{i^\prime_p,j^\prime_p}\right)^{\ast} dU   =\sum_{\sigma,\tau\in S_p}\text{Wg}([\sigma\tau^{-1}];d)\, ,\tag{1} 
$$

with just checking

$$
\text{Wg}(\sigma\tau^{-1};d)=\frac{1}{(p!)^2} \sum_{\lambda} \frac{\chi^\lambda(1)^2 \chi^\lambda([\sigma\tau^{-1}])}{s_{\lambda,d}}\, ,\tag{2}
$$
where $U$ is a Haar-random $d\times d$ unitary matrix, $dU$ is the Haar measure over $U(d)$, and  $[\sigma]$ is the class of element $\sigma$.  The sum in Eq. (1) is a sum over all $\sigma\in S_p$ and all the $\tau\in S_p$ so that

$$   (i^\prime_{\sigma(1)},\ldots,i^\prime_{\sigma(p)})=(i_1,\ldots,i_p)\, ,\\
(j^\prime_{\tau(1)},\ldots,j^\prime_{\tau(p)})=(j_1,\ldots,j_p)\, ,
$$

with the integral $0$ is the $i',i$, $j'$ or $j$ strings have different lengths.  In other words, the integral on the left of Eq. (1) is a sum of Weingarten functions.  In the expression for $\text{Wg}$, $\chi^\lambda(\mathbf{1})$ is the dimension of irrep $\lambda$ of $S_p$, $\chi^\lambda(\sigma)$ is the character of element $\sigma$ in the irrep $\lambda$, $s_{\lambda,d}$ is the dimension of irrep $\lambda$ of $U(d)$, and the sum over $\lambda$ is a sum over all partitions of $p$.  Alternate derivations and properties, including generalizations to the orthogonal and the symplectic groups, can be found in Refs. [6-9].

This note describes the workings of a Mathematica code to quickly evaluate Weingarten functions as given in Eq. (2).  The code is built on an implementation of the Murnaghan-Nakayama rule for the characters irreducible representations of the symmetric group $S_p$, provided to me  by Dr. Justin Kulp [10].  The dimensionality factor $s_{\lambda,d}$ can be computed in the usual way using the hook-rule [11].  To speed up calculations it was found useful to write a dedicated function to evaluate the dimension of the irrep $\lambda$ of $S_p$, also using the hook-rule.

## The main functions: eWg and cWg

The code contains two main functions: eWg and cWg.  The functions differ in their arguments, as explained in details a little later.  The LHS of Eq. (2) is seen to depend on the product $\sigma\tau^{-1}$ of elements in $S_p$ but in fact the RHS depends only on the class of $\sigma\tau^{-1}$.  Furthermore, the sum on the RHS of Eq. (1) is a sum over products of elements in $S_p$ so it is convenient to define  eWg when group elements are used as inputs, and cWg for use when classes of elements are used as inputs.

The cWg functions takes as input a class and a dimension parameter $d$:

    cWg[class_,d_]

and returns the RHS of Eq. (2).  Classes are partitions inside curly brackets.  
Thus for instance: 

>  In[1]:= cWg[\{3,1\},d]
> 
>  Out[1]:= $\frac{-3+2 d^2}{(-3+d) (-2+d) (-1+d) {d}^{2}(1+d) (2+d) (3+d)}$
>  
is the Weingarten function for the integration of $U(d)$ functions, for the class $\{3,1\}$. 

The function  eWg has slightly different inputs:

    eWg[sigma_,p_,d_]


Mathematica constructs elements of $S_p$ in the form of cycles, such as Cycles[\{\{2,3,4\}\}].  However, it is not possible to full determine from the information in Cycles if this is an element of $S_4, S_5$ *etc* hence the need to supply the additional parameter $p$ to indicate this is an element in $S_p$.  The output is of course the same as  cWg  if the order $p$ and an element in the class are given as inputs:

> In[2]:= eWg[Cycles[\{\{2,3,4\}\}],4,d]
> 
> Out[2]:= $\frac{-3+2d^2}{(d-3)(d-2)(d-1)d^2 (1+d) (2+d) (3+d)}$


The dimension $d$ of the unitaries can be passed directly to either functions:

> In[3]:= cWg[Cycles[\{\{3,1\}\}],6]
> 
> Out[3]:= $\frac{23}{362880}$


## Auxiliary functions

### murnNaka

This function comes from the code of [10].  The inputs are a partition and a class:

      murnNaka[partition_,class_]

The output is the character of the elements in "class" of the irrep of the symmetric group labelled by "partition".  For example:

> In[4]:= murnNaka[\{3,1\},\{1,1,1,1\}]
> 
> Out[4]:= $3$

which is also the dimension of irrep $(3,1)$.

### getClass

This function inputs an element in $S_p$ and returns the class of this element:

    getClass[p_,cycles_]

For example:

> In[5]:= getClass[4,Cycles[\{\{1,2,3\}\}]]
> 
> Out[5]:= $\{3,1\}$

### snDimension

This function inputs a "partition"


    snDimension[partition_]

and returns the dimension of the irrep of $S_p$ labelled by a partition.  It uses the hook-length formula.  For instance 

> In[6]:= snDimension[\{5,4,2,1,1,1\}]
> Out[6]:= $63063$

snDimension gives the same result as murnNaka when the class is $(1,\ldots,1)$ but is considerably faster than murnNaka when $p$ is large and the partition contains many parts.  

### udDimension

This functions outputs the dimension of the $U(d)$ labelled by "ppartition":

    udDimension[ppartition_, d_].

Thus, for the partition $(5,4,2,1,1,1)$, we have

> In[7]:= udDimension[\{5,4,2,1,1,1\},d]
> 
> Out[7]:=  $\frac{(-5+d) (-4+d) (-3+d) (-2+d) (-1+d)^{2}{d}^{2} (1+d)^{2} (2+d)^{2} (3+d) (4+d)}{1382400}$
>
> In[8]:= udDimension[\{5,4,2,1,1,1\},8]
> 
> Out[8]:   $873180$

The dimension $d$ of the unitary must be greater than or equal to the number of parts in the partition, or alternatively the $d$ must be greater than or equal to  the number of rows in the Young diagram associated with the partition, else there is no irrep for this partition and the function returns $0$:

> In[9]:=  udDimension[\{5,4,2,1,1,1\},4]
> Out[9]:   $0$
> 

### gClass

This functions outputs the number of elements in the class of $S_p$ labelled by "partition":

     gClass[partition_]

Thus, for the partition $(5,1,1)$ of $S_7$, we have

> In[10]:=  gClass[\{5,1,1\},4]
> 
> Out[10]:   $504$


## Tables of Weingarten functions for $n\le 6$

 $S_2$: integrals of the type $\displaystyle\int dU U_{i\alpha}U_{j\beta} U^*_{k\eta} U^*_{\ell \nu}\, .$


|Class| Wg  |
|--|--|
|$\{2\}$ |$\displaystyle\frac{1}{(d-1) d (d+1)}$|
| $\{1,1\}$ | $\displaystyle\frac{1}{(d-1)(d+1)}$ |



 $S_3$: integrals of the type $\displaystyle\int dU U_{i\alpha}U_{j\beta}U_{k\eta} 
U^*_{m\mu} U^*_{\ell \nu}U^*_{p\kappa}$.

|Class  |Wg  |
|--|--|
|  \{3\} | $\displaystyle\frac{2}{(d-2) (d-1) d (d+1) (d+2)}$ |
|\{2,1\}|$-\displaystyle\frac{1}{(d-2) (d-1) (d+1) (d+2)}$|
| \{1,1,1\}  |$\displaystyle\frac{d^2-2}{(d-2) (d-1) d (d+1) (d+2)}$ |





$S_4$:  

|Class|Wg  |
|--|--|
| \{4\}   |$-\displaystyle\frac{5}{(d-3) (d-2) (d-1) d (d+1) (d+2) (d+3)}$  |
|\{3,1\} |$\displaystyle\frac{2 d^2-3}{(d-3) (d-2) (d-1) d^2 (d+1) (d+2) (d+3)}$ |
|\{2,2\}|$\displaystyle\frac{d^2+6}{(d-3) (d-2) (d-1) d^2 (d+1) (d+2) (d+3)}$|
| \{2,1,1\}| $-\displaystyle\frac{1}{(d-3) (d-1) d (d+1) (d+3)}$|
|\{1,1,1,1\} |$\displaystyle\frac{d^4-8 d^2+6}{(d-3) (d-2) (d-1)d^2 (d+1) (d+2)(d+3)}$|




$S_5$

|Class|Wg  |
|--|--|
|  \{5\}    |$\displaystyle\frac{14}{(d-4) (d-3) (d-2) (d-1) d (d+1) (d+2) (d+3)(d+4)}$  |
|\{4,1\}  |$\displaystyle\frac{24-5 d^2}{(d-4) (d-3) (d-2) (d-1) d^2 (d+1) (d+2)(d+3)(d+4)}$ |
| \{3,2\}|$-\displaystyle\frac{2 \left(d^2+12\right)}{(d-4) (d-3) (d-2) (d-1)d^2(d+1)(d+2)(d+3)(d+4)}$|
|  \{3,1,1\} | $\displaystyle\frac{2}{(d-4) (d-2) (d-1) d (d+1) (d+2) (d+4)}$|
| \{2,2,1\}  |$\displaystyle\frac{-d^4+14 d^2-24}{(d-4) (d-3) (d-2) (d-1) d^2  (d+1) (d+2) (d+3) (d+4)}$|
| \{1,1,1,1,1\}|$\displaystyle\frac{d^4-20 d^2+78}{(d-4) (d-3) (d-2) (d-1) d (d+1) (d+2) (d+3) (d+4)}$|

 $S_6$: 
| Class |Wg  |
|--|--|
|  \{6\} | $-\displaystyle\frac{42}{(d-5) (d-4) (d-3) (d-2) (d-1) d (d+1) (d+2)(d+3) (d+4) (d+5)}$ |
| \{5,1\} | $\displaystyle\frac{14 \left(d^2-10\right)}{(d-5) (d-4) (d-3) (d-2)(d-1) d^2 (d+1) (d+2) (d+3) (d+4) (d+5)}$|
|  \{4,2\}  | $\displaystyle\frac{5 \left(d^4+15 d^2+8\right)}{(d-5)(d-4)(d-3)(d-2)(d-1)^2 d^2(d+1)^2 (d+2) (d+3) (d+4) (d+5)}$|
|\{4,1,1\}|$\displaystyle\frac{13-5 d^2}{(d-5) (d-3) (d-2) (d-1)^2 d (d+1)^2   (d+2) (d+3) (d+5)}$|
| \{3,3\} |$\displaystyle\frac{4 \left(d^4+29 d^2-90\right)}{(d-5) (d-4) (d-3)(d-2) (d-1)^2 d^2 (d+1)^2 (d+2) (d+3) (d+4) (d+5)}$|
| \{3,2,1\}| $\displaystyle\frac{-2 d^2-13}{(d-5) (d-4) (d-2) (d-1)^2 d (d+1)^2(d+2) (d+4) (d+5)}$ |
| \{3,1,1,1\} | $\displaystyle\frac{2 d^6-51 d^4+229 d^2-60}{(d-5) (d-4) (d-3)(d-2) (d-1)^2 d^2 (d+1)^2 (d+2) (d+3) (d+4) (d+5)}$|
| \{2,2,2\}|$\displaystyle\frac{-d^4-d^2-358}{(d-5) (d-4) (d-3) (d-2) (d-1)^2 d   (d+1)^2 (d+2) (d+3) (d+4) (d+5)}$|
 |\{2,2,1,1\} |$\displaystyle\frac{d^4-3 d^2+10}{(d-5) (d-3) (d-2) (d-1)^2 d^2   (d+1)^2 (d+2) (d+3) (d+5)}$|
 | \{2,1,1,1,1\} | $\displaystyle\frac{-d^4+24 d^2-38}{(d-5) (d-4) (d-2) (d-1)^2 d(d+1)^2 (d+2) (d+4) (d+5)}$|
 \{1,1,1,1,1,1\} | $\displaystyle\frac{d^8-41 d^6+458 d^4-1258 d^2+240}{(d-5)   (d-4) (d-3) (d-2) (d-1)^2 d^2 (d+1)^2 (d+2) (d+3) (d+4) (d+5)}$

## Acknowledgements

This work was supported by [ACFAS](https://www.acfas.ca/) through their Programme de coopération en recherche dans la francophonie canadienne.  I would like to  thank  Dr. Nicolás Quesada and his group for their hospitality during two visits at École Polytechnique in Montréal and for testing the package, and David Amaro-Alcala for checking the source code of the $\beta$-version of the package and for helpful comments.

Please cite as: 

```
@software{deGuiseWeingit,
  author = {Hubert de Guise},
  title = {{weingarten_package: A Mathematica package for unitary Weingarten functions}},
  url = {https://github.com/hdeguise/Weingarten_calculus},
  year = {2024},
}
```

## Bibliography

1. P A Mello. “Averages on the unitary group and applications to the problem of disordered conductors”. In: Journal of Physics A: Mathematical and General 23.18 (1990), p. 4061.
    
2.  Stuart Samuel. “U(N) Integrals, 1/N, and the De Wit–’t Hooft anomalies”. In: Journal of Mathematical Physics 21.12 (1980), pp. 2695–2703.
    
3.  Sho Matsumoto. “Moments of a single entry of circular orthogonal ensembles and Weingarten calculus”. In: Letters in Mathematical Physics 103 (2013), pp. 113–130.
    
4.  Don Weingarten. “Asymptotic behavior of group integrals in the limit of infinite rank”. In: Journal of Mathematical Physics 19.5 (1978), pp. 999–1001.
    
5.  Benoît Collins and Piotr Śniady. “Integration with respect to the Haar measure on unitary, orthogonal and symplectic group”. In: Communications in Mathematical Physics 264.3 (2006), pp. 773–795.
    
6.  Marcel Novaes. “Elementary derivation of Weingarten functions of classical Lie groups”. In: arXiv preprint arXiv:1406.2182 (2014).
    
7.  Thomas Gorin. “Integrals of monomials over the orthogonal group”. In: Journal of Mathematical Physics 43.6 (2002), pp. 3342–3351.
    
8.  T Gorin and GV López. “Monomial integrals on the classical groups”. In: Journal of Mathematical Physics 49.1 (2008).
    
9.  Benoít Collins and Sho Matsumoto. “On some properties of orthogonal Weingarten functions”. In: Journal of Mathematical Physics 50.11 (2009).
    
10.  Justin Kulp. Murnahan-Nakayama. 2014. url: https://github.com/justinkulp/MurnaghanNakayama.
    
11.  Rutwig Campoamor-Stursberg and Michel Rausch De Traubenberg. Group theory in physics: a practitioner’s guide. World Scientific.
   

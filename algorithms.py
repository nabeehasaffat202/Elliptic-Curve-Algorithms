import math

def pollards_rho(n, x_0):
    i = 0
    xarray=[]
    
    def poly(x_i,n):
        return (x_i**2 + 1) %n

    xarray.append(x_0)
    
    #calculating x_1 and x_2
    x = poly(x_0,n)
    xarray.append(x)
    x = poly(x,n)
    xarray.append(x)
    i+=1
    
    while (math.gcd( int(xarray[2*i]- xarray[i]) , n ) == 1 ):

        #calculating up to new x_2i
        x = poly(x,n)
        xarray.append(x)
        x = poly(x,n)
        xarray.append(x)
        i+=1

    print(f'N = {n}, i = {i}, p = {math.gcd(xarray[2*i] - xarray[i], n )} ')
    print(f'i/sqrtN = {i/math.sqrt(n)}\n')
    
#pollards_rho(2201,2)
#pollards_rho(9409613,2)
#pollards_rho(1782886219,2)


def discrete_log_solver(g, h ,p):
    i = 0
    x = 1
    y = 1
    a = 0
    b = 0
    ga = 0
    de = 0

    print(f'g,h,p = {g},{h},{p}')

    def f(x, g, h, p):
        x_i1 = 0
        if (0 <= x and x < p/3):
            x_i1 = g*x
        if (p/3 <= x  and x < 2*p/3 ):
            x_i1 = x**2
        if (2*p/3 <= x  and x < p ):
            x_i1 = h*x
        return x_i1 % p

    def alpha(a, x, g, h, p):
        a_i1 = 0
        if ( 0 <= x and  x < p/3):
            a_i1 = a+1
        if (p/3 <= x  and x < 2*p/3 ):
            a_i1 = 2*a
        if (2*p/3 <= x  and x < p ):
            a_i1 = a
        return a_i1 % (p-1)
    
    def beta(b,x, g, h, p):
        b_i1 = 0
        if ( 0 <= x and x < p/3):
            b_i1 = b
        if (p/3 <= x  and x < 2*p/3 ):
            b_i1 = 2*b
        if (2*p/3 <= x  and x < p ):
            b_i1 = b+1
        return b_i1 % (p-1)


    while (i==0 or x != y ):
        a = alpha( a,x, g,h,p)
        b = beta( b, x, g, h, p)
        
        ga = alpha( alpha( ga , y , g,h,p) , f(y,g,h,p), g,h,p)
        
        de = beta(  beta(  de , y , g, h, p) , f(y,g,h,p) , g, h, p)
        
        x = f(x,g,h,p)
        y = f( f(y,g,h,p) ,g,h,p)
        i += 1

    print(f'i = {i}, x = {x}, a = {a} , b = {b} , y= {y},  ga = {ga}, de = {de} ')

    u = (a - ga ) % (p-1)
    v = ( de - b) % (p-1)

    d =  math.gcd(v,p-1)

    def euc(a, b):
            if a == 0:
                return 0,1
            x1,y1 = euc( b%a, a)
            x = y1 - (b//a) * x1
            y = x1
            return x,y
    
    if ( d == 1):
        sol = u* euc( v, p-1)[0] % (p-1)
        print(f'sol is {sol}')
        
    else:

        s = euc( v, p-1)[0]
        w = s*u
        solset=[]
        k = 0
        print(f'w = {w}, d= {d}, s={s}')
        while (k<d):
            solset.append( w/d + k*(p-1)/d )
            k += 1

        k -= 1
    
        while ( k >= 0):
            if (h == pow(g, int(solset[k]),  p) ):
                print(f'sol is {solset[k]}')
            k -= 1

        print(f'i/sqrt(p*pi/2) = {i/math.sqrt(p*math.pi/2)}\n')

        
#discrete_log_solver(2,2495,5011)
#discrete_log_solver(17,14226,17959)
#discrete_log_solver(29,5953042, 15239131)

def elliptic_curve_factorization(s):
    discrete_log_solver(3, 12*s**2 , 12*s**2+6*s+1 )

    print(f's = {s}, 2*s^2 + s = {2*s**2+s}\n')    

#elliptic_curve_factorization(3)
#elliptic_curve_factorization(5)
#elliptic_curve_factorization(7)


def q3h3(x,y, N, a):
    i=1
    Px = x
    Py = y

    #a = a, b = N
    def euc(a,b):
        if (math.gcd(a,b) != 1):
            print(f'ERROR')
        if a == 0:
            return 0,1
        x1,y1 = euc( b%a, a)
        x = y1 - (b//a) * x1
        y = x1
        return x,y
        
    l = ( (3*Px**2+a)*( euc(2*Py, N)[0] ) ) % N
    i+=1
    Pxold = Px
    Px = (l**2 - 2*Px ) %N
    Py = (l*(Pxold-Px) - Py) %N


    def lam(x,y, Px, Py, N,a):
        lam = (Py-y)* (euc( (Px-x), N)[0] )
        return lam % N 


    while ( math.gcd( Px-x ,N) == 1 ):
        l = lam(x,y,Px,Py, N,a)
        Pxold = Px
        Px = ( l**2 - 2*Px ) % N
        Py = ( l*(Pxold-Px) - Py ) % N
        i += 1

    print(f'i = {i-1}, P = ({Px}, {Py})')
    print(f'N = {N}, Px-x = {Px-x} ')
    print(f'factors = {math.gcd( Px-x ,N) } , { N /math.gcd( Px-x ,N) } \n')

q3h3(2,12, 26167, 4)
q3h3(1,1, 1386493, 3)
q3h3(7,4, 28102844557, 18)


def q1h3(x,y, N, a):
    i=1
    Px = x
    Py = y

    #a = a, b = N
    def euc(a,b):
        if (math.gcd(a,b) != 1):
            print(f'ERROR')
        if a == 0:
            return 0,1
        x1,y1 = euc( b%a, a)
        x = y1 - (b//a) * x1
        y = x1
        return x,y
        
    l = ( (3*Px**2+a)*( euc(2*Py, N)[0] ) ) % N
    i+=1
    Pxold = Px
    Px = (l**2 - 2*Px ) %N
    Py = (l*(Pxold-Px) - Py) %N


    def lam(x,y, Px, Py, N,a):
        lam = (Py-y)* (euc( (Px-x), N)[0] )
        return lam % N 


    while ( i < 876):
        l = lam(x,y,Px,Py, N,a)
        Pxold = Px
        Px = ( l**2 - 2*Px ) % N
        Py = ( l*(Pxold-Px) - Py ) % N
        i += 1

    print(f'i = {i-1}, P = ({Px}, {Py})')
    
#q1h3(2,2575, 2671, 171)
#q1h3(2,96, 2671, 171)


        

    
    




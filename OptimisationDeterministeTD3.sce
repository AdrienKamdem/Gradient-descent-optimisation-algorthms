//////////////////// Adrien Kamdem ING2 GMI2 ///////////////////////

/////////////////////////////////////////EXO3

// qts1
function POFQ = PasOptimalFonctionQuadr(A,b,x0,itermax,epsilon)
    
    // initialisations
    k=0;
    dk = -(A*x0+b)
    pk = (norm(dk)^2)/(dk'*A*dk)
    xk = x0
    xk1 = xk + pk*dk
    
    
    //stop conditions
    while norm(xk1 - xk)> epsilon*norm(xk)& k < itermax
        xktmp = xk1
        xk = xktmp
        dk = -(A*xktmp+b)
        pk = (norm(dk)^2)/(dk'*A*dk)
        xk1 = xktmp + pk*dk
        k = k+1
          
    end
    
    POFQ = pk
    
endfunction

// qts2
function POFMC = PasOptimalFonctionMoindr(A,b,x0,kmax,epsilon)
    
    // initialisations
    xk = x0
    dk = -A'*(A*xk-b)
    pk = (norm(dk)^2)/(norm(A*dk)^2)
    xk1 = xk+pk*dk
    
    // stop conditions
    while k < kmax | norm(xk1-xk) > epsilon*norm(xk)
        xtmp = xk1
        xk = xktmp
        dk = -A'*(A*xktmp-b) 
        pk = (norm(dk)^2)/(norm(A*dk)^2)
        xk1 = xktmp + pk*dk
    end
    
    POFMC = pk
    
endfunction

// qts3
function armijo = PasOptimalArmijo(A,b,c,x0,epsilon,T,B1,w1)
    
    // cas ou f est une fonction quadratique
    // initialisations
    T = 10^(-2)// ]0,1[
    w1 = 10^(-4) // < 0.5 && ]0,1[
    k=0
    xk = x0
    dk = -(A*xk+b)
    
    //stop conditions
    while (1/2)*(xk+pk*dk)'*A*(xk+pk*dk)+b'*(xk+pk*dk) > (1/2)*(xk'*A*xk)+b'*xk + c + w1*pk *(dk'*(A*xk*dk))
        xktmp = xk1
        xk = xktmp
        dk = -(A*xk+b)
        pk = norm(dk)^2/(dk'*A*dk)
        xk1 = xktmp + pk*dk
        k = k+1
    end
    
    armijo = pk
endfunction


/////////////////////////////////////////EXO4

// qts 1 
// pk = norm(dk)^2/norm(A*dk)^2
// on remarque que f est quadratique avec c=0, b= (-3,-3) et A = ([10 0],[0 1])

// qts 2
function [POFQbis] = PasOptimalFonctionQbis(A,b,x0,itermax,epsilon)
    
    // initialisations
    k=0;
    dk = -(A*x0+b);
    pk = (norm(dk)^2)/(dk'*A*dk)
    xk = x0
    xk1 = xk + pk*dk
    resd = list()
    resp = list()
    resx = list()
    
    //stop conditions
    while norm(xk1 - xk)> epsilon*norm(xk)& k < itermax
        xktmp = xk1
        xk = xktmp
        dk = -(A*xktmp+b)
        pk = (norm(dk)^2)/(dk'*A*dk)
        xk1 = xktmp + pk*dk
        k = k+1
        resd($+1)=dk
        resp($+1)=pk
        resx($+1)=xk   
    end
    
    POFQbis = list(resp,resx,resp)
    
endfunction

// qts3 (ne fonctionne pas on ne sait pas pourquoi)
function [armijobis] = PasOptimalAbis(A,b,c,x0,epsilon,B1)
    
    // cas ou f est une fonction quadratique
    
    // initialisations
    T = 0.7// ]0,1[
    w1 = 10^(-4) // < 0.5 && ]0,1[
    k=0
    xk = x0
    dk = -(A*xk+b)
    pk = 1 // = p0
    resAp = list()
    xk1 = xk + pk*dk
    
    // stop conditions
    while (1/2)*(xk+pk*dk)'*A*(xk+pk*dk)+b'*(xk+pk*dk) > (1/2)*(xk'*A*xk)+b'*xk + c + w1*pk *(dk'*(A*xk+b))
        xktmp = xk1
        xk = xktmp
        dk = -(A*xk+b)
        pk = norm(dk)^2/(dk'*A*dk)
        xk1 = xktmp + pk*dk
        k = k+1
        resAp($+1) = pk
        resXk($+1) = xk
    end
    
    armijobis = list(resAp,resXk)
endfunction

import numpy as n
import logging

lib_logger = logging.getLogger("SEBI-CF.lib")

def Bw(hr,L,z0h,z0m):
    alphah=0.12
    betah=125.
    if (hr/L) > 0:
        return -2.2*n.log(1+hr/L)
    else:
        if z0m < (alphah/betah*hr):
            return -n.log(alphah) + Psym(-alphah*hr/L) - Psym(-z0m/L)
        else:
            return n.log(hr/(betah*z0m)) + Psym(-betah*z0m/L) - Psym(-z0m/L)

def Cw(hr,L,z0h,z0m):
    alphah = 0.12
    betah = 125.
    
    res = n.zeros([len(z0m),z0m.size/len(z0m)])
    
    for i in n.arange(0,len(z0m),1):
        for j in n.arange(0,z0m.size/len(z0m),1):
            if (hr[i,j]/L[i,j]) > 0:
                res[i,j] = -7.6*n.log(hr[i,j]/L[i,j])
            elif z0m[i,j] < (alphah/betah*hr[i,j]):
                res[i,j] = -n.log(alphah) + Psyh(-alphah*hr[i,j]/L[i,j]) - Psyh(-z0h[i,j]/L[i,j])
            else:
                res[i,j] = n.log(hr[i,j]/(betah*z0m[i,j])) + Psyh(-betah*z0m[i,j]/L[i,j]) - Psyh(-z0h[i,j]/L[i,j])
    return res

def cleanup(x,label):
    try:
        if label == "ndvi":
            search = n.where(x<0.)
            x[search] = abs(x[search])
            search = n.where(x>1.)
            x[search] = n.nan
            search = n.where(x==0.)
            x[search] = 0.001

        elif label == "albedo":
            search = n.where(x<0.)
            x[search] = n.nan
            search = n.where(x>1.)
            x[search] = n.nan

        elif label == "ts":
            search = n.where(x<250.)
            x[search] = n.nan
            search = n.where(x>340.)
            x[search] = n.nan

        elif label == "swdw":
            search = n.where(x<0.)
            x[search] = n.nan
            search = n.where(x>1400.)
            x[search] = n.nan

        elif label == "lwdw":
            search = n.where(x<0.)
            x[search] = n.nan
            search = n.where(x>500.)
            x[search] = n.nan

        elif label == "RnDaily":
            search = n.where(x<0.)
            x[search] = n.nan
            search = n.where(x>1400.)
            x[search] = n.nan

        elif label == "G0_Rn":
            search = n.where(x<0.)
            x[search] = n.nan
            search = n.where(x>1.)
            x[search] = n.nan
    except :
        pass
    
    else:
        pass

    try:
        log = str(len(x)) + " pixels cleaned up for " + label
    except:
        log = "unknown cleanup criteria for " + label

    return x,log

def delta(es,t):
    return n.log(10.)*7.5*237.3/(237.3+(t-273.16))**2 * es

def deltad(cp,red,rho,Rn,G0):
    return red*(Rn-G0)/(rho*cp)

def deltaw(cp,delta,e,es,gamma,rew,rho,Rn,G0):
    return (rew*(Rn-G0)/(rho*cp)-((es-e)/gamma))/(1+(delta/gamma))

def downscaling(x,myproj):
    #try:
        res = n.zeros([myproj.gridNb[0],myproj.gridNb[1]])
        for i in n.arange(0,myproj.gridNb[0]):
            for j in n.arange(0,myproj.gridNb[1]):
                xgrid = x[i*myproj.pixPerGrid:(i+1)*myproj.pixPerGrid,
                    j*myproj.pixPerGrid:(j+1)*myproj.pixPerGrid]
                if n.nansum(xgrid) == n.nan or myproj.pixPerGrid**2 - nancount(xgrid) == 0:
                    res[i,j] = n.nan
                else:
                    res[i,j] = n.nansum(xgrid)/(myproj.pixPerGrid**2 - nancount(xgrid))
        return res
    #except Exception, err:
        #print "\nDownscaling error:"
        #print "pixels per grid=",myproj.pixPerGrid
        #try:
        #    print "Grid cell ", i,j
        #    print "nansum(xgrid)=",n.nansum(xgrid)
        #    print "pixels with data=", myproj.pixPerGrid**2 - nancount(xgrid)
        #except:
        #    pass
        #sys.stderr.write('ERROR: %s\n' % str(err))
        #return 1

def eact(p,q,rd,rv):
    '''Pelgrum 2000, eq. 2.7, p.30'''
    return p*q*rv/rd

def esat(es0,t):
    '''Tetens (cf. Guyot 1999, eq.3.13, p.109)'''
    if n.nansum(t)/(t.size-nancount(t)) > -40. and n.nansum(t)/(t.size-nancount(t)) < 100.:
        return es0*pow(10,((7.5*(t))/(237.3+t)))
    elif n.nansum(t)/(t.size-nancount(t)) > 233. and n.nansum(t)/(t.size-nancount(t)) < 373.:
        return es0*pow(10,((7.5*(t-273.16))/(237.3+t-273.16)))
    else:
        print "WARNING: unit or value inconsitency in e(sat)"

def G0(fc,Rn):
    return Rn*(0.05 + (1-fc)*(0.3-0.05))

#def getHist(x,label):
#    p.clf()
#    txtStats = label#getStats(x,label)
#    p.hist(x,normed=True)
#    p.title(txtStats)
#    p.savefig(label+'-hist.png')
    
#def getPreview(x,label):
#    p.clf()
#    p.imshow(x,interpolation='nearest',vmin=n.nanmin(x),vmax=n.nanmax(x))
#    p.colorbar()
#    p.title(label)
#    p.savefig(label+'-preview.png')

def getStats(x,label):
    lib_stats_logger = logging.getLogger("SEBI-CF.lib.stats")
    i = nancount(x)
    #return label + ' : min=' + str(n.nanmin(x)) + ' max=' + str(n.nanmax(x)
    #                            ) + ' avg=' + str(n.nansum(x)/x.size) + ' NanCount=' + str(i)
    try:
        lib_stats_logger.info('Mean ' + label + ' = ' + str("%10.4f" % (n.nansum(x)/(x.size-i))) + ' (NaN count is ' + str(i) +')')
    except ZeroDivisionError:
        lib_stats_logger.info('error on ' + label + ' (type is ' + str(x.dtype) + '): contains ' + str(n.nansum(x)) + ' Nan for a total of ' + str(x.size) + ' pixels')


def H(cp,delta_a,k,ra,rho,ustar):
    return delta_a*k*ustar*rho*cp/ra

def hv2z0m(hv):
    (i,j)=n.where(hv!=0)
    z0m = hv-hv
    z0m[i,j] = (0.136 * hv[i,j])
    return z0m

def kB(cd,ct,fc,k,hg,hr,hs,hv,lai,ndvi,p0,pr,tr,ur,z0m):
    lib_kb_logger = logging.getLogger("SEBI-CF.lib.kb")

    '''Massman kB-1 model'''
    
    fc = setMask(fc,0.,0.,1.)
    fc = substitute(fc,0.,0.001)
    #fc = substitute(fc,0.,n.NaN)
    lai = substitute(lai,0.,0.001)
    #lai = substitute(lai,0.,n.NaN)
    
    U = 0.32 - 0.264*n.exp(-15.1*cd*lai) # U=u*/u(h)
    N = cd*lai/(2.*(U**2.))
    
    # N might be zero
    (i,j)=N.nonzero()
    
    d_h = N-N
    d_h[i,j] = 1. -1/(2*N[i,j])*(1-n.exp(-2*N[i,j]))
    
    z0m_h = (1-d_h)*n.exp(-k/U)
    #z0m_h = z0m / hv
    Prdt = 0.71 # Prandtl number
    #pr = p0*((1-((0.0065*(hg+hr))/288.15))**5.256) INPUT!!
    t = tr*((pr/p0)**.286)

    # Some code added to avoid divide by zero
    nu = pr - pr # just to get an array of the same dimension
    search = n.nonzero(nu) # identify numerical
    nu[search] = (1.327*(p0/pr[search])*(t[search]/273.16)**1.81)*1.E-5
    z0m = ndvi2z0m(ndvi)
    hv = z0m_h / z0m
    d0 = z0m2d0(z0m)
    ustar = u2ustar(d0,hr,k,ur,z0m)

    # wind speed at canopy height
    uh = z0m - z0m # just to get an array of the same dimension

    # The following is done only to present calculation of log(0)
    #hflow = hv-d0
    #search = n.where(hflow=0.)
    #loghflow =


    uh = ur*((n.log(hv-d0)-n.log(z0m))/(n.log(hr-d0)-n.log(z0m)))
    #Restar = U*uh*hv/nu
    Restar = hs*ustar/nu
    Ctstar = (Prdt**(-2./3))*(Restar**(-1/2.))
   
    kb_dec_a = N-N
    kb_dec_a[i,j] = (k*cd)/(4*ct*U[i,j]*(1-n.exp(-N[i,j]/2)))
    
    kb_dec_b = (k*U*z0m_h)/Ctstar
    kb_dec_b = nan2flt(kb_dec_b,0.) # Replace NaN with 0.
    
    #kbs_1 = (2.46*(((ur*k/n.log(hr/hs))*hs/nu)**(1/4))) - n.log(7.4)
    kbs_1 = 2.46*(Restar**0.25) - n.log(7.4)
    kb_1 = kb_dec_a*fc**2 + kb_dec_b*fc*2*(1-fc) + kbs_1*(1-fc)**2
    
    ref = n.where(fc<=0.1)
    kb_1[ref] = kbs_1[ref]
    
    getStats(U,'U')
    getStats(N,'N')
    getStats(d_h,'d_h')
    getStats(z0m_h,'z0m_h')
    getStats(pr,'pr')
    getStats(t,'t')
    getStats(nu,'nu')
    getStats(d0,'d0')
    getStats(ustar,'ustar')
    getStats(uh,'uh')
    getStats(Restar,'Restar')
    getStats(Ctstar,'Ctstar')
    getStats(kb_dec_a,'kb_dec_a')
    getStats(kb_dec_b,'kb_dec_b')
    getStats(kbs_1,'kbs_1')
    
    # !!!
    #kb_1 = kb_1-kb_1+4.0
    # !!!
    
    z0h = z0m / n.exp(kb_1)
    
    return kb_1, z0h

def L(cp,delta_a,g,k,Le,rho,Rn,G0,ustar,state='none'):
    if state == 'dry':
        return -(ustar**3 *rho)/(k*g*(Rn-G0)/(delta_a*cp))
    elif state == 'wet':
        return -(ustar**3 *rho)/(k*g*0.61*(Rn-G0)/Le)

def nan2flt(x,xnew):
    '''Remplace NaN in x with xnew value and return the array'''
    ref = n.nonzero(x-x)
    x[ref] = xnew
    return x

def nancount(x):
    '''Count NaN occurences in an array. Based on the result of numpy.nonzero() on
    an array containing both zeros and NaN: only NaN are indexed'''
    
    ref = n.nonzero(x-x)
    i = x[ref]
    return len(i)

def ndvi2emi(ndvi):
    '''Following van de Griend & Owe, 1993'''
    return 1.009 + 0.047*n.log(ndvi)

def ndvi2fc(ndvi):
    '''Calculate fc = f(ndvi). Another option could be fc = 1.318*ndvi + 0.01877'''
    return 1 - ((ndvi - n.nanmax(ndvi))/(n.nanmin(ndvi) - n.nanmax(ndvi)))**0.4631

def ndvi2lai(ndvi):
    return (ndvi*(1.+ndvi)/(1.-ndvi))**0.5

def ndvi2z0m(ndvi):
    '''Following Moran. Also 0.005 + 0.5*((ndvi/ndvi.max()))**2.5'''
    return n.exp(-5.2 + 5.3*ndvi)

def ps_sea2gnd(ps,hg):
    """
    ps MUST be in Pa
    """
    return ps*((1-((0.0065*hg)/288.15))**5.256)

def Psyh(y):
    return ((1-0.057)/0.78)*n.log((0.33+(n.abs(y)**0.78))/0.33)

def Psym(y):
    a = 0.33
    b = 0.41
    psy0 = -n.log(a) + 3**(1/2.)*b*a**(1/3.)*n.pi/6
    if y > b**(-3.):
        y = b**(-3.)
    x = (1/a*n.abs(y))
    return n.log(a+n.abs(y)) - 3*b*n.abs(y)**(1/3.)+b*a**(1/3.)/2*n.log((1+x)**2 /(1-x+x**2)) + 3**(1/2.)*b*a**(1/3.)*n.arctan((2*x-1)/(3**(1/2.))) + psy0

def ra(d0,hr,z0h):
    return n.log((hr-d0)/z0h)

def re(Cw,d0,hr,k,ustar,z0h):
    re = (n.log((hr-d0)/z0h)-Cw)/(k*ustar)
    re_alt = (n.log((hr-d0)/z0h))/(k*ustar)
    search = n.where(re <= 0)
    re[search] = re_alt[search]
    return re

def rho(e,p,q,rd,t):
    '''Brutsaert 1982, eq.3.6, p.38'''
    return (p/(rd*t))*(1-(0.378*e/p))

def Rn(albedo,emi,lwdw,sigma,swdw,ts):
    return swdw*(1.-albedo) + lwdw - (1-emi)*lwdw - emi*sigma*(ts**4.)
    
def setMask(x,xnew,vmin,vmax):
    '''Replace in x values out of [vmin:vmax] with xnew'''
    ref = n.where(x<=vmin)
    x[ref] = xnew
    ref = n.where(x>=vmax)
    x[ref] = xnew
    return x

def substitute(x,v,vnew):
    '''Replace in x the value v with vnew'''
    ref = n.where(x==v)
    x[ref] = vnew
    return x

def stabFunc(x):
    """Returns the set of 3 equations, with
    x[0] = ustar
    x[1] = H
    x[2] = L
    Get the variables from one single file
     vars[0] = ur
     vars[1] = hr
     vars[2] = d0
     vars[3] = z0m
     vars[4] = Bw
     vars[5] = delta_a
     vars[6] = rho
     vars[7] = ra
     vars[8] = Ta(pot)
     vars[9] = Le_i
     vars[10] = Rn-G0"""
    try:
        vars = n.fromfile('tmp000')
        #out = [vars[0] - x[0]/0.4*n.log((vars[1]-vars[2])/vars[3]) - vars[4]]
        #out.append(x[1] - vars[5]*0.4*x[0]*vars[6]*1005./vars[7])
        #out.append(x[2] + x[0]**3 *vars[7]/(0.4*9.81*(x[1]/(vars[9]*1005.)+0.61*((vars[9]-x[1])/vars[10]))))
        #return out
        return [(vars[0] - x[0]/(0.4*n.log((vars[1]-vars[2])/vars[3]) - vars[4])),(x[1] - vars[5]*0.4*x[0]*vars[6]*1005./vars[7]),(x[2] + ((x[0]**3) * vars[7])/(0.4*9.81*(x[1]/(vars[8]*1005.)+0.61*((vars[10]-x[1])/vars[9]))))]
    except:
        print "Error in vars ", vars
        
def tpot(cp,p,p0,q,rd,t):
    kr = rd*(1-0.23*q)/cp
    return t*((p0/p)**kr)
    
def u2ustar(d0,h,k,u,z0m):
    """
    Calculate u* from the wind speed
    """
    #TODO: Unit testing for u2ustar

    ustar = z0m - z0m # to get an empty array of proper dimension
    ratio = (h-d0)/z0m
    search = n.nonzero(ratio)
    ustar[search] = (u[search]*k)/(n.log(ratio[search]))
    return ustar

def z0m2d0(z0m):
    return z0m*4.9

def z0m2hv(z0m):
    return z0m / 0.136

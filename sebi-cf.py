#! /usr/bin/env python

"""
Stand-alone SEBI-CF prototype for CEOP-AEGIS WP3

Use the demo test case:
- demo.xml file
- demo directory containing a test dataset

To run (sebi-cf.py should be +x):
./sebi-cf.py demo.xml
Have a look at sebi-cf.log for logs

# TODO:
- centralize cleanup statements in lib.cleanup (also for outputs)
- bound kB outputs
"""
__author__ = "Jerome Colin"
__author_email__ = "first-name dot last-name @ lsiit.u-strasbg.fr"
__date__ = "2011-01-20"
__version__ = "0.2.1 (sneezy)"

print 'SEBI-CF ' + __version__ + "\n Use -h for help"
import congrid
import getopt
import lib as l
import numpy as n
import os
import progressBar
import Project
#import pylab as p
import rebin
from scipy.optimize import fsolve, broyden2
import sys
import logging
import time
import warnings

# You may consider reset it to "default" or "always" for testing purposes ONLY
warnings.simplefilter("ignore")

#create logger with "spam_application"
logger = logging.getLogger("SEBI-CF")
logger.setLevel(logging.DEBUG)

# clear previous logs
try:
    os.remove("sebi-cf.log")
except OSError:
    pass

#create file handler and set level to debug
fh = logging.FileHandler("sebi-cf.log")
fh.setLevel(logging.DEBUG)

#create console handler and set level to error
ch = logging.StreamHandler()
ch.setLevel(logging.ERROR)

#create formatter
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")

#add formatter to fh
fh.setFormatter(formatter)

#add formatter to ch
ch.setFormatter(formatter)

#add fh to logger
logger.addHandler(fh)

#add ch to logger
logger.addHandler(ch)
logger.info("Launching SEBI-CF")

def main():
    """Parse command line options."""
    try:
        opts, args = getopt.getopt(sys.argv[1:], "h", ["help"])
    except getopt.error, msg:
        print msg
        print "for help use --help"
        sys.exit(2)
    # process options
    for o, a in opts:
        if o in ("-h", "--help"):
            print __doc__
            sys.exit(0)
    # process arguments
    for arg in args:
        process(arg) # process() is defined elsewhere

def process(fproject):
    """Backbone: Tasks are processed the following order:
    - Get time
    - Setting common constants
    - Instantiate a new project
    - Load images
    - Try load optional images, else use surrogate
    - Crop images
    - Generate atmospheric forcing layers
    - Input diagnostic
    - kB-1 sequence
    - Radiative budget sequence
    - Downscaling
    - Stability sequence
    - Low resolution at-limit parameters
    - Upscaling
    - External resistances and gradients
    - Saving to files
    - Get time and save logs
    
    """
    # Dedicated logger
    process_logger = logging.getLogger("SEBI-CF.Process")

    try:
        # Get time
        time0 = time.time()
        # Setting common constants
        cd = 0.2
        ct = 0.01
        cp = 1005.
        es0 = 610.7 #Pa
        g = 9.81
        hs = 0.009
        k = 0.4
        p0 = 101325.
        pdtl = 0.71
        gamma = 67.
        rd = 287.04
        rv = 461.05
        sigma = 5.678E-8
        
        # Instantiate a new project
        myproj = Project.project(fproject)
        process_logger.info("Instantiate a new project")
        
    except Exception, err:
        sys.stderr.write('ERROR: %s\n' % str(err))
      
    #try:
    # Calculate scales
    myproj.setGrids()

    # Load images
    widgets = [' Loading surface data:    ', progressBar.Percentage(), ' ', progressBar.Bar(marker='-',left='[',right=']'),
                       ' ', ' ', ' ', ' ']
    pbar = progressBar.ProgressBar(widgets=widgets, maxval=8).start()

    albedo = myproj.read(myproj.falbedo)
    pbar.update(1)
    ts = myproj.read(myproj.flst)
    pbar.update(2)
    ndvi = myproj.read(myproj.fndvi)
    pbar.update(3)

    # Cleanup
    albedo,log = l.cleanup(albedo,"albedo")
    process_logger.info(log)
    ts,log = l.cleanup(ts,"ts")
    process_logger.info(log)
    ndvi,log = l.cleanup(ndvi,"ndvi")
    process_logger.info(log)


    # Try load optional images, else use surrogate
    if myproj.femissivity != 'undef':
        emi = myproj.read(myproj.femissivity)
        process_logger.info("Emissivity image found")
    else:
        emi = l.ndvi2emi(ndvi)

    pbar.update(4)

    if myproj.ffc != 'undef':
        fc = myproj.read(myproj.ffc)
        process_logger.info("Fc image found")
    else:
        fc = l.ndvi2fc(ndvi)

    pbar.update(5)

    if myproj.fhv != 'undef':
        hv = myproj.read(myproj.fhv)
        process_logger.info("Hv image found")
        hv = l.substitute(hv, 0., 0.01)
        z0m = l.hv2z0m(hv)
    else:
        z0m = l.ndvi2z0m(ndvi)
        z0m = l.substitute(z0m, 0., 0.000001)
        hv = l.z0m2hv(z0m)
        hv = l.substitute(hv, 0., 0.01)
    pbar.update(6)

    if myproj.flai != 'undef':
        lai = myproj.read(myproj.flai)
        process_logger.info("LAI image found")
    else:
        lai = l.ndvi2lai(ndvi)
    pbar.update(7)

    if myproj.mask != 'undef':
        mask = myproj.read(myproj.mask)
        process_logger.info("Mask image found")
    else:
        mask = ndvi-ndvi

    if myproj.RnDaily != 'undef':
        RnDaily = myproj.read(myproj.RnDaily)
        process_logger.info("RnDaily image found")
    else:
        mask = ndvi-ndvi



    pbar.update(8)

    pbar.finish()

    # Crop images
    albedo = albedo[0:myproj.imgDims[0], 0:myproj.imgDims[1]]
    ts = ts[0:myproj.imgDims[0], 0:myproj.imgDims[1]]
    ndvi = ndvi[0:myproj.imgDims[0], 0:myproj.imgDims[1]]
    emi = emi[0:myproj.imgDims[0], 0:myproj.imgDims[1]]
    fc = fc[0:myproj.imgDims[0], 0:myproj.imgDims[1]]
    hv = hv[0:myproj.imgDims[0], 0:myproj.imgDims[1]]
    z0m = z0m[0:myproj.imgDims[0], 0:myproj.imgDims[1]]
    lai = lai[0:myproj.imgDims[0], 0:myproj.imgDims[1]]
    mask = mask[0:myproj.imgDims[0], 0:myproj.imgDims[1]]
    if myproj.RnDaily != 'undef':
        RnDaily = RnDaily[0:myproj.imgDims[0], 0:myproj.imgDims[1]]

    # Generate atmospheric forcing layers
    widgets = [' Loading PBL data:        ', progressBar.Percentage(), ' ', progressBar.Bar(marker='-',left='[',right=']'),
                       ' ', ' ', ' ', ' ']
    pbar = progressBar.ProgressBar(widgets=widgets, maxval=10).start()

    if myproj.atmmode == '1D':
        hg = n.zeros(myproj.imgDims,dtype=float) + myproj.hg
        pbar.update(1)
        hr = n.zeros(myproj.imgDims,dtype=float) + myproj.hr
        pbar.update(2)
        lwdw = n.zeros(myproj.imgDims,dtype=float) + myproj.lwdw
        pbar.update(3)
        pg = n.zeros(myproj.imgDims,dtype=float) + myproj.pg
        pbar.update(4)
        pr = n.zeros(myproj.imgDims,dtype=float) + myproj.pr
        pbar.update(5)
        qg = n.zeros(myproj.imgDims,dtype=float) + myproj.qg
        pbar.update(6)
        qr = n.zeros(myproj.imgDims,dtype=float) + myproj.qr
        pbar.update(7)
        swdw = n.zeros(myproj.imgDims,dtype=float) + myproj.swdw
        pbar.update(8)
        tr = n.zeros(myproj.imgDims,dtype=float) + myproj.tr
        pbar.update(9)
        ur = n.zeros(myproj.imgDims,dtype=float) + myproj.ur
        pbar.update(10)
    if myproj.atmmode == '2D':
        hg = myproj.read(myproj.fhg)
        pbar.update(1)
        hr = myproj.read(myproj.fhr)
        pbar.update(2)
        lwdw = myproj.read(myproj.flwdw)
        pbar.update(3)
        pg = myproj.read(myproj.fpg)
        pbar.update(4)
        pr = myproj.read(myproj.fpr)
        pbar.update(5)
        qg = myproj.read(myproj.fqg)
        pbar.update(6)
        qr = myproj.read(myproj.fqr)
        pbar.update(7)
        swdw = myproj.read(myproj.fswdw)
        pbar.update(8)
        tr = myproj.read(myproj.ftr)
        pbar.update(9)
        ur = myproj.read(myproj.fur)
        pbar.update(10)

        # Additional cleanup
        swdw,log = l.cleanup(swdw,"swdw")
        process_logger.info(log)
        lwdw,log = l.cleanup(lwdw,"lwdw")
        process_logger.info(log)
        if myproj.RnDaily != 'undef':
            RnDaily,log = l.cleanup(lwdw,"RnDaily")
            process_logger.info(log)

        # Crop images
        hg = hg[0:myproj.imgDims[0], 0:myproj.imgDims[1]]
        hr = hr[0:myproj.imgDims[0], 0:myproj.imgDims[1]]
        lwdw = lwdw[0:myproj.imgDims[0], 0:myproj.imgDims[1]]
        pg = pg[0:myproj.imgDims[0], 0:myproj.imgDims[1]]
        pr = pr[0:myproj.imgDims[0], 0:myproj.imgDims[1]]
        qg = qg[0:myproj.imgDims[0], 0:myproj.imgDims[1]]
        qr = qr[0:myproj.imgDims[0], 0:myproj.imgDims[1]]
        swdw = swdw[0:myproj.imgDims[0], 0:myproj.imgDims[1]]
        tr = tr[0:myproj.imgDims[0], 0:myproj.imgDims[1]]
        ur = ur[0:myproj.imgDims[0], 0:myproj.imgDims[1]]

        if myproj.pressureUnit == "hPa":
            pg = pg*100.
            pr = pr*100.
        if myproj.pressureLevel == "SL":
            pg = l.ps_sea2gnd(pg,hg)

        #TMP TMP TMP
        #search = n.where(mask == 5.5)
        #ur[search] = ur[search]*1.2

    pbar.finish()
    # Apply mask
    if myproj.mask != 'undef':
        search = n.where(mask==0)

        albedo[search] = n.nan
        ts[search] = n.nan
        ndvi[search] = n.nan
        emi[search] = n.nan
        fc[search] = n.nan
        hv[search] = n.nan
        z0m[search] = n.nan
        lai[search] = n.nan
        hg[search] = n.nan
        hr[search] = n.nan
        lwdw[search] = n.nan
        pg[search] = n.nan
        pr[search] = n.nan
        qg[search] = n.nan
        qr[search] = n.nan
        swdw[search] = n.nan
        tr[search] = n.nan
        ur[search] = n.nan
        if myproj.RnDaily != 'undef':
            RnDaily[search] = n.nan


#    # Input diagnostic
    widgets = [' Input diagnostic:        ', progressBar.Percentage(), ' ', progressBar.Bar(marker='-',left='[',right=']'),
                       ' ', ' ', ' ', ' ']
    pbar = progressBar.ProgressBar(widgets=widgets, maxval=18).start()
    l.getStats(albedo,'albedo')
    pbar.update(1)
    l.getStats(ts,'ts')
    pbar.update(2)
    l.getStats(ndvi,'ndvi')
    pbar.update(3)
    l.getStats(emi,'emi')
    pbar.update(4)
    l.getStats(fc,'fc')
    pbar.update(5)
    l.getStats(hv,'hv')
    pbar.update(6)
    l.getStats(z0m,'z0m')
    pbar.update(7)
    l.getStats(lai,'lai')
    pbar.update(8)
    l.getStats(hg,'hg')
    pbar.update(9)
    l.getStats(hr,'hr')
    pbar.update(10)
    l.getStats(lwdw,'lwdw')
    pbar.update(11)
    l.getStats(pg,'pg')
    pbar.update(12)
    l.getStats(pr,'pr')
    pbar.update(13)
    l.getStats(qg,'qg')
    pbar.update(14)
    l.getStats(qr,'qr')
    pbar.update(15)
    l.getStats(swdw,'swdw')
    pbar.update(16)
    l.getStats(tr,'tr')
    pbar.update(17)
    l.getStats(ur,'ur')
    pbar.update(18)

    pbar.finish()

    # kB-1 sequence
    widgets = [' Running kB-1 model:      ', progressBar.Percentage(), ' ', progressBar.Bar(marker='-',left='[',right=']'),
                       ' ', ' ', ' ', ' ']
    pbar = progressBar.ProgressBar(widgets=widgets, maxval=100).start()
    if myproj.kbMode == "Massman":
        process_logger.info("Launching kB-1 model")
        kB_1,z0h = l.kB(cd,ct,fc,k,hg,hr,hs,hv,lai,ndvi,p0,pr,tr,ur,z0m)
        
    else:
        kB_1 = n.zeros(myproj.imgDims,dtype=float) + 4.
        z0h = z0m / n.exp(kB_1)

    l.getStats(kB_1,'kB_1')
    l.getStats(z0h,'z0h')

    pbar.update(100)
    pbar.finish()

    # Radiative budget
    widgets = [' Radiative budget:        ', progressBar.Percentage(), ' ', progressBar.Bar(marker='-',left='[',right=']'),
                       ' ', ' ', ' ', ' ']
    pbar = progressBar.ProgressBar(widgets=widgets, maxval=3).start()
    myproj.logs += '\n\nRadiative budget:'

    Rn = l.Rn(albedo,emi,lwdw,sigma,swdw,ts)
    l.getStats(Rn,'Rn')
    pbar.update(1)

    G0 = l.G0(fc,Rn)
    l.getStats(G0,'G0')
    pbar.update(2)

    G0_Rn = G0/Rn
    l.getStats(G0_Rn,'G0_Rn')

    G0_Rn,log = l.cleanup(G0_Rn,"G0_Rn")
    myproj.logs += log

    pbar.update(3)

    pbar.finish()

    # Downscaling
    widgets = [' Downscaling:             ', progressBar.Percentage(), ' ', progressBar.Bar(marker='-',left='[',right=']'),
                       ' ', ' ', ' ', ' ']
    pbar = progressBar.ProgressBar(widgets=widgets, maxval=12).start()
    myproj.logs += '\n\nDownscaling:'

    low_z0m = l.downscaling(z0m,myproj)
    l.getStats(low_z0m,'low_z0m')
    pbar.update(1)

    low_z0h = l.downscaling(z0h,myproj)
    l.getStats(low_z0h,'low_z0h')
    pbar.update(2)

    low_ts = l.downscaling(ts,myproj)
    l.getStats(low_ts,'low_ts')
    pbar.update(3)

    low_Rn = l.downscaling(Rn,myproj)
    l.getStats(low_Rn,'low_Rn')
    pbar.update(4)

    low_G0 = l.downscaling(G0,myproj)
    l.getStats(low_G0,'low_G0')
    pbar.update(5)

    low_ur = l.downscaling(ur,myproj)
    l.getStats(low_ur,'low_ur')
    pbar.update(6)

    low_hr = l.downscaling(hr,myproj)
    l.getStats(low_hr,'low_hr')
    pbar.update(7)

    low_pr = l.downscaling(pr,myproj)
    l.getStats(low_pr,'low_pr')
    pbar.update(8)

    low_pg = l.downscaling(pg,myproj)
    l.getStats(low_pg,'low_pg')
    pbar.update(9)

    low_qr = l.downscaling(qr,myproj)
    l.getStats(low_qr,'low_qr')
    pbar.update(10)

    low_qg = l.downscaling(qg,myproj)
    l.getStats(low_qg,'low_qg')
    pbar.update(11)

    low_tr = l.downscaling(tr,myproj)
    l.getStats(low_tr,'low_tr')
    pbar.update(12)

    pbar.finish()

    # Stability sequence
    widgets = [' Stability sequence:      ', progressBar.Percentage(), ' ', progressBar.Bar(marker='-',left='[',right=']'),
                       ' ', ' ', ' ', ' ']
    pbar = progressBar.ProgressBar(widgets=widgets, maxval=myproj.gridNb[0]*myproj.gridNb[1]).start()

    myproj.logs += '\n\nStability sequence'

    low_d0 = l.z0m2d0(low_z0m)
    l.getStats(low_d0,'low_d0')
    ustar_i = l.u2ustar(low_d0,low_hr,k,low_ur,low_z0m)
    l.getStats(ustar_i,'ustar_i')
    ra_i = l.ra(low_d0,low_hr,low_z0h)
    l.getStats(ra_i,'ra_i')
    low_delta_a = l.tpot(cp,low_pg,p0,low_qg,rd,low_ts) - \
                    l.tpot(cp,low_pr,p0,low_qr,rd,low_tr)
    l.getStats(low_delta_a,'low_delta_a')
    low_es = l.esat(es0,low_ts)
    l.getStats(low_es,'low_es')
    low_e = l.eact(low_pg,low_qg,rd,rv)
    l.getStats(low_e,'low_e')
    low_rho = l.rho(low_e,low_pg,low_qg,\
                    rd,low_ts)
    l.getStats(low_rho,'low_rho')
    H_i = l.H(cp,low_delta_a,k,ra_i,low_rho,ustar_i)
    l.getStats(H_i,'H_i')
    delta = l.delta(low_es,low_ts)
    l.getStats(delta,'delta')
    Le_i = (delta*rd*(low_ts)**2)/(0.622*low_es)
    l.getStats(Le_i,'Le_i')
    L_i = (-ustar_i**3 *low_rho)/(k*g*0.61*(low_Rn-low_G0)/Le_i)
    l.getStats(L_i,'L_i')
    #pbar.update(1)

    # Modif >><<

    H_ic = low_delta_a * k * low_rho * cp
    L_ic = -low_rho *cp * (low_ts * (1 + 0.61 * low_qr)) / (k * g)
    ustar_i = k * low_ur / n.log((low_hr-low_d0) / low_z0m)
    H_i = H_ic * ustar_i / n.log((low_hr-low_d0) / low_z0h)
    H_target = H_i
    L_i = L_i-L_i

    # >><<

    # Starting the iterative sequence
    vars = n.zeros([11,1],dtype=float)
    # Variables for output
    slvUstar = n.zeros([myproj.gridNb[0],myproj.gridNb[1]],dtype=float)
    slvH = n.zeros([myproj.gridNb[0],myproj.gridNb[1]],dtype=float)
    slvL = n.zeros([myproj.gridNb[0],myproj.gridNb[1]],dtype=float)
    iterator = n.zeros([myproj.gridNb[0],myproj.gridNb[1]],dtype=float)

    if myproj.iterate == "True":
        for i in n.arange(0,myproj.gridNb[0]):
            for j in n.arange(0,myproj.gridNb[1]):

                stabDif = 10.
                stabLoops = 0

                while stabDif > 0.01 and stabLoops < 100:
                    L_i[i,j] = L_ic[i,j] * (ustar_i[i,j]**3) / H_i[i,j]
                    ustar_i[i,j] = k* low_ur[i,j] / (n.log((low_hr[i,j]-low_d0[i,j]) / low_z0m[i,j]) - l.Bw(low_hr[i,j],L_i[i,j],low_z0h[i,j],low_z0m[i,j]))
                    H_i[i,j] = H_ic[i,j] * ustar_i[i,j] / (n.log((low_hr[i,j]-low_d0[i,j]) / low_z0h[i,j]) - (-7.6*n.log(low_hr[i,j]/L_i[i,j])))
                    stabDif   = n.abs(H_target[i,j] - H_i[i,j])
                    H_target[i,j] = H_i[i,j]
                    stabLoops+=1

                slvUstar[i,j] = ustar_i[i,j]
                slvH[i,j] = H_i[i,j]
                slvL[i,j] = L_i[i,j]
                iterator[i,j] = stabLoops

                ## Grid stability functions
                #Cw = -7.6*n.log(low_hr[i,j]/L_i[i,j])
                #ra = n.log((low_hr[i,j]-low_d0[i,j])/low_z0h[i,j])
                #Bw = l.Bw(low_hr[i,j],L_i[i,j],low_z0h[i,j],low_z0m[i,j])
                ## Prepare the file to provide to l.stabFunc
                #vars[0] = low_ur[i,j]
                #vars[1] = low_hr[i,j]
                #vars[2] = low_d0[i,j]
                #vars[3] = low_z0m[i,j]
                #vars[4] = Bw
                #vars[5] = low_delta_a[i,j]
                #vars[6] = low_rho[i,j]
                #vars[7] = ra
                #vars[8] = l.tpot(cp,low_pg[i,j],p0,low_qg[i,j],rd,0.5*(low_ts[i,j]+low_tr[i,j]))
                #vars[9] = Le_i[i,j]
                #vars[10] = low_Rn[i,j] - low_G0[i,j]
                #vars.tofile('tmp000')
                #
                ##slvUstar[i,j],slvH[i,j],slvL[i,j] = fsolve(l.stabFunc,[ustar_i[i,j],H_i[i,j],L_i[i,j]],warning=False)
                #try:
                #    slvUstar[i,j],slvH[i,j],slvL[i,j] = broyden2(\
                #            l.stabFunc,[ustar_i[i,j],H_i[i,j],L_i[i,j]],iter=40,verbose=False)
                #except(OverflowError):
                #    slvUstar[i,j] = ustar_i[i,j]
                #    slvH[i,j] = H_i[i,j]
                #    slvL[i,j] = L_i[i,j]

                pbar.update(myproj.gridNb[1]*i+j)

        # add some stats
        l.getStats(slvUstar,'slvUstar')
        l.getStats(slvH,'slvH')
        l.getStats(slvL,'slvL')
        l.getStats(iterator,'iterator')
        pbar.finish()

    else:
        # 2010-02-05: TMP TMP TMP TMP TMP TMP TMP TMP TMP TMP TMP TMP TMP TMP TMP TMP TMP TMP TMP TMP TMP TMP
        slvUstar = ustar_i
        slvH = H_i
        slvL = L_i

        # TMP TMP TMP TMP TMP TMP TMP TMP TMP TMP TMP TMP TMP TMP TMP TMP TMP TMP TMP TMP TMP TMP TMP TMP TMP

    # Low resolution at-limit parameters
    low_Ld = l.L(cp,low_delta_a,g,k,Le_i,low_rho,low_Rn,low_G0,slvUstar,state='dry')
    low_Lw = l.L(cp,low_delta_a,g,k,Le_i,low_rho,low_Rn,low_G0,slvUstar,state='wet')
    low_Cwd = l.Cw(low_hr,low_Ld,low_z0h,low_z0m)
    low_Cww = l.Cw(low_hr,low_Lw,low_z0h,low_z0m)
    low_Cw = l.Cw(low_hr,slvL,low_z0h,low_z0m)
    l.getStats(low_Ld,'low_Ld')
    l.getStats(low_Lw,'low_Lw')
    l.getStats(low_Cwd,'low_Cwd')
    l.getStats(low_Cww,'low_Cww')

    # Upscaling
    widgets = [' Upscaling:               ', progressBar.Percentage(), ' ', progressBar.Bar(marker='-',left='[',right=']'),
                       ' ', ' ', ' ', ' ']
    pbar = progressBar.ProgressBar(widgets=widgets, maxval=6).start()

    # TMP TMP TMP TMP TMP TMP TMP TMP TMP TMP TMP TMP TMP TMP TMP TMP TMP TMP TMP TMP TMP TMP TMP TMP TMP
    low_rho = l.nan2flt(low_rho,n.nansum(low_rho)/low_rho.size)
    slvUstar = l.nan2flt(slvUstar,n.nansum(slvUstar)/slvUstar.size)
    low_Cwd = l.nan2flt(low_Cwd,n.nansum(low_Cwd)/low_Cwd.size)
    low_Cww = l.nan2flt(low_Cww,n.nansum(low_Cww)/low_Cww.size)
    low_Cw = l.nan2flt(low_Cw,n.nansum(low_Cw)/low_Cw.size)
    slvL = l.nan2flt(slvL,n.nansum(slvL)/slvL.size)
    # TMP TMP TMP TMP TMP TMP TMP TMP TMP TMP TMP TMP TMP TMP TMP TMP TMP TMP TMP TMP TMP TMP TMP TMP TMP

    rho = congrid.congrid(low_rho,[myproj.imgDims[0],myproj.imgDims[1]],method='spline',minusone=True)
    pbar.update(1)
    ustar = congrid.congrid(slvUstar,[myproj.imgDims[0],myproj.imgDims[1]],method='spline',minusone=True)
    pbar.update(2)
    Cwd = congrid.congrid(low_Cwd,[myproj.imgDims[0],myproj.imgDims[1]],method='spline',minusone=True)
    pbar.update(3)
    Cww = congrid.congrid(low_Cww,[myproj.imgDims[0],myproj.imgDims[1]],method='spline',minusone=True)
    pbar.update(4)
    Cw = congrid.congrid(low_Cw,[myproj.imgDims[0],myproj.imgDims[1]],method='spline',minusone=True)
    pbar.update(5)
    L = congrid.congrid(slvL,[myproj.imgDims[0],myproj.imgDims[1]],method='spline',minusone=True)
    pbar.update(6)

    pbar.finish()

    # External resistances and gradients
    widgets = [' Processing SEBI:         ', progressBar.Percentage(), ' ', progressBar.Bar(marker='-',left='[',right=']'),
                       ' ', ' ', ' ', ' ']
    pbar = progressBar.ProgressBar(widgets=widgets, maxval=4).start()

    myproj.logs += '\n\nExernal resistances and gradients:'

    d0 = l.z0m2d0(z0m)
    es = l.esat(es0,0.5*(ts+tr))
    e = l.eact(0.5*(pg+pr),0.5*(qg+qr),rd,rv)
    delta = l.delta(es,0.5*(ts+tr))
    pbar.update(1)

    # dry limit
    red = l.re(Cwd,d0,hr,k,ustar,z0h)
    deltad = l.deltad(cp,red,rho,Rn,G0)
    l.getStats(Cwd,'Cwd')
    l.getStats(red,'red')
    l.getStats(deltad,'deltad')
    pbar.update(2)

    # wet limit
    rew = l.re(Cww,d0,hr,k,ustar,z0h)
    deltaw = l.deltaw(cp,delta,e,es,gamma,rew,rho,Rn,G0)
    l.getStats(Cww,'Cww')
    l.getStats(rew,'rew')
    l.getStats(deltaw,'deltaw')
    pbar.update(3)

    # actual conditions
    Cw = l.Cw(hr,L,z0h,z0m)
    re = l.re(Cw,d0,hr,k,ustar,z0h)
    deltaa = l.tpot(cp,pg,p0,qg,rd,ts) - l.tpot(cp,pr,p0,qr,rd,tr)
    SEBI = (deltaa/re - deltaw/rew) / (deltad/red - deltaw/rew)
    ef = 1 - SEBI
    search = n.where(ef > 1.)
    ef[search] = 1.
    LE = (Rn-G0)*ef
    H = Rn - LE - G0

    # relative evap (alternate)
    Hd = Rn - G0
    Hw = ((Rn - G0) - (rho*cp / rew) * es / gamma) / (1.0 + delta / gamma)
    Ha = rho*cp * deltaa / re

    search = n.where(Ha>Hd)
    Ha[search] = Hd[search]

    search = n.where(Ha<Hw)
    Ha[search] = Hw[search]




    Le_re = 1. - ((Ha - Hw) / (Hd - Hw))

    search = n.where(Hd<=Hw)
    Le_re[search] = 1.

    Le_fr = Le_re * (Rn - G0 - Hw) / (Rn - G0)

    search = n.where((Rn-G0)<=0)
    Le_fr[search] = 1.

    if myproj.RnDaily != 'undef':
        ETdaily = Le_fr * (RnDaily * (1 - G0_Rn)) * 24*3600 / 2.45E6

    pbar.update(4)
    pbar.finish()

    widgets = [' Output statistics        ', progressBar.Percentage(), ' ', progressBar.Bar(marker='-',left='[',right=']'),
                       ' ', ' ', ' ', ' ']
    pbar = progressBar.ProgressBar(widgets=widgets, maxval=12).start()

    l.getStats(L,'L')
    pbar.update(1)
    l.getStats(Cw,'Cw')
    pbar.update(2)
    l.getStats(re,'re')
    pbar.update(3)
    l.getStats(deltaa,'deltaa')
    pbar.update(4)
    l.getStats(SEBI,'SEBI')
    pbar.update(5)
    l.getStats(ef,'ef')
    pbar.update(6)
    l.getStats(Ha,'Ha')
    pbar.update(7)
    l.getStats(Hd,'Hd')
    pbar.update(8)
    l.getStats(Hw,'Hw')
    pbar.update(9)
    l.getStats(Le_re,'Le_re')
    pbar.update(10)
    l.getStats(Le_fr,'Le_fr')
    pbar.update(11)
    l.getStats(ETdaily,'ETdaily')
    pbar.update(12)
    pbar.finish()



    # Saving to files
    widgets = [' Saving to files:         ', progressBar.Percentage(), ' ', progressBar.Bar(marker='-',left='[',right=']'),
                       ' ', ' ', ' ', ' ']
    pbar = progressBar.ProgressBar(widgets=widgets, maxval=8).start()

    # Data layers:
    #ascIO.ascWritter(myproj.path+myproj.prefix+'ef',
    #                 ef,myproj.imgDims[1],myproj.imgDims[0],xllcorner=myproj.xllcorner,yllcorner=myproj.yllcorner,cellsize=myproj.cellsize,NODATA_value='nan')
    myproj.writeRaw(ef,'ef')
    pbar.update(1)

    #ascIO.ascWritter(myproj.path+myproj.prefix+'kb',
    #                 kB_1,myproj.imgDims[1],myproj.imgDims[0],xllcorner=myproj.xllcorner,yllcorner=myproj.yllcorner,cellsize=myproj.cellsize,NODATA_value='nan')
    myproj.writeRaw(kB_1,'kb')
    pbar.update(2)

    #ascIO.ascWritter(myproj.path+myproj.prefix+'LE',
    #                 LE,myproj.imgDims[1],myproj.imgDims[0],xllcorner=myproj.xllcorner,yllcorner=myproj.yllcorner,cellsize=myproj.cellsize,NODATA_value='nan')
    myproj.writeRaw(LE,'LE')
    pbar.update(3)

    #ascIO.ascWritter(myproj.path+myproj.prefix+'H',
    #                 H,myproj.imgDims[1],myproj.imgDims[0],xllcorner=myproj.xllcorner,yllcorner=myproj.yllcorner,cellsize=myproj.cellsize,NODATA_value='nan')
    myproj.writeRaw(H,'H')
    pbar.update(4)

    #ascIO.ascWritter(myproj.path+myproj.prefix+'Rn',
    #                 Rn,myproj.imgDims[1],myproj.imgDims[0],xllcorner=myproj.xllcorner,yllcorner=myproj.yllcorner,cellsize=myproj.cellsize,NODATA_value='nan')
    myproj.writeRaw(Rn,'Rn')
    pbar.update(5)

    #ascIO.ascWritter(myproj.path+myproj.prefix+'G0',
    #                 G0,myproj.imgDims[1],myproj.imgDims[0],xllcorner=myproj.xllcorner,yllcorner=myproj.yllcorner,cellsize=myproj.cellsize,NODATA_value='nan')
    myproj.writeRaw(G0,'G0')
    pbar.update(6)

    myproj.writeRaw(Le_fr,'Le_fr')
    pbar.update(7)

    if myproj.RnDaily != 'undef':
        myproj.writeRaw(ETdaily,'ETdaily')

    pbar.update(8)

    pbar.finish()

    print 'Done...'
        
    #except Exception, err:
        #log.exception('Error from process():')
        #myproj.logs += log.exception('Error from process():')
        #process_logger.info("ERROR: process aborted")
        
    #finally:
        # Get time and save logs
    myproj.logs += '\n\nTotal processing time (s):' + str(time.time()-time0)
    myproj.logs += '\n\nSnowWhite ' + __version__
    fd=open(myproj.path + myproj.prefix + 'log.txt','w')
    fd.write(myproj.logs)
    fd.close()
    print 'Log file: '+ myproj.prefix + 'log.txt'
        
if __name__ == "__main__":
    main()

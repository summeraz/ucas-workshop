import numpy as np
import math
import random
import string
    
def md_start(nmol, rho, temp, eps, sig, mass, delt, cutoff, 
          x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4,
          px0, py0, pz0, px1, py1, pz1, px2, py2, pz2, px3, py3, pz3, px4, py4, pz4):
    # 
    # initlaize simulation
    #
    rstar = rho * sig*sig*sig / mass
    rstar = rstar * 1.0E-24
    boxlen = (nmol / rstar)**(1.0/3.0)
    vol = boxlen**3
    halfl = boxlen / 2.0
    tstar = temp / eps
    cutf2 = cutoff*cutoff
    n = int(round((nmol/4)**(1.0/3.0)))

    print rstar, boxlen, vol, halfl, tstar, cutf2, n

    nby2 = n*2
    nsby2 = 2*n*n
    delta = 2./nby2
    gap = 2.*delta
    strt = -1. + delta/2.0
    strt1 = strt + delta
    print nby2, nsby2, delta, gap, strt, strt1

    # initialize molecules in fcc lattice
    for ilay in range(1,nby2,2):
        for k in range(1,nby2+1,1):
            for j in range(1,n+1,1):
                idisk=((k/2)*2)/k
                disx=idisk*delta
                i=j+(k-1)*n+(ilay-1)*nsby2
                x0[i-1]=(strt+disx+gap*(j-1)+1.0)*halfl
                y0[i-1]=(strt+delta*(k-1)+1.0)*halfl
                z0[i-1]=(strt+delta*(ilay-1)+1.0)*halfl
    for ilay in range(2,nby2+1,2):
        for k in range(1,nby2+1,1):
            for j in range(1,n+1,1):
                idisk=-((k/2)*2)/k
                disx=idisk*delta
                i=j+(k-1)*n+(ilay-1)*nsby2
                x0[i-1]=(strt1+disx+gap*(j-1)+1.0)*halfl
                y0[i-1]=(strt+delta*(k-1)+1.0)*halfl
                z0[i-1]=(strt1+delta*(ilay-2)+1.0)*halfl
    #
    # now assign velocities
    #
    sumx = 0.0
    sumy = 0.0
    sumz = 0.0
    arg = 3.0*tstar
    rtr = delt * math.sqrt(arg)
    for i in range(nmol):
        xx = random.random()
        yy = random.random()
        zz = random.random()
        xyz = rtr/math.sqrt(xx*xx + yy*yy + zz*zz)
        x1[i] = xx * xyz
        y1[i] = yy * xyz
        z1[i] = zz * xyz
        sumx = sumx + x1[i]
        sumy = sumy + y1[i]
        sumz = sumz + z1[i]
    #
    #   subtract center of mass motion.
    #
    vcx = sumx/nmol
    vcy = sumy/nmol
    vcz = sumz/nmol
    for i in range(nmol):
        x1[i] = x1[i] - vcx
        y1[i] = y1[i] - vcy
        z1[i] = z1[i] - vcz
        px0[i] = x1[i]/delt
        py0[i] = y1[i]/delt
        pz0[i] = z1[i]/delt
    #
    # calculate  kinetic energy per molecule
    #
    tsum = 0.0
    for i in range(nmol):
        tsum = tsum + px0[i]*px0[i] + py0[i]*py0[i] + pz0[i]*pz0[i]
    ekin = 0.5*tsum/nmol
    # 
    # work out temperature of configuration
    #
    mdtemp = 2.0*ekin/3.0
    # 
    # calculate scaling factor to change temperature of configuration to required temperature
    #
    scaling_factor = math.sqrt(tstar/mdtemp)
    #
    # rescale velocities
    #
    for i in range(nmol):
        px0[i] = scaling_factor*px0[i]
        py0[i] = scaling_factor*py0[i]
        pz0[i] = scaling_factor*pz0[i]
        x1[i] = px0[i]*delt
        y1[i] = py0[i]*delt
        z1[i] = pz0[i]*delt
        
    return tstar, rstar, cutf2, boxlen

def loop(nstep,nblock,nvt,blkavr,tstar,nmol,rho,temp,eps,sig,mass,delt,cutoff,cutf2,boxlen,time,
         x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4,
         px0, py0, pz0, px1, py1, pz1, px2, py2, pz2, px3, py3, pz3, px4, py4, pz4,
         traj_time, traj_ke, traj_pe, traj_tot):
    #
    # constants for 4th order predictor-corrextor
    #
    f01, f21, f31, f41 = 0.348611111, 0.916666667,0.333333333, 0.041666667
    #
    #  the following variables are used to accumulate variables to form
    #   block averages.  iblk controls which block average is being performed
    #
    eksum, ek2sum, epsum, ep2sum, psum, p2sum = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    halfl=boxlen/2.0
    vol=boxlen**3
    for istep in range(1,nstep+1):
        time = time + delt
        traj_time[istep-1]=time
        for i in range(nmol):
            #
            # predictor step
            #
            x0[i] = x0[i] + x1[i] + x2[i] + x3[i] + x4[i]
            y0[i] = y0[i] + y1[i] + y2[i] + y3[i] + y4[i]
            z0[i] = z0[i] + z1[i] + z2[i] + z3[i] + z4[i]
            x1[i] = x1[i] + 2*x2[i] + 3*x3[i] + 4*x4[i]
            y1[i] = y1[i] + 2*y2[i] + 3*y3[i] + 4*y4[i]
            z1[i] = z1[i] + 2*z2[i] + 3*z3[i] + 4*z4[i]
            x2[i] = x2[i] + 3*x3[i] + 6*x4[i]
            y2[i] = y2[i] + 3*y3[i] + 6*y4[i]
            z2[i] = z2[i] + 3*z3[i] + 6*z4[i]
            x3[i] = x3[i] + 4*x4[i]
            y3[i] = y3[i] + 4*y4[i]
            z3[i] = z3[i] + 4*z4[i]
            px0[i] = px0[i] + px1[i] + px2[i] + px3[i] + px4[i]
            py0[i] = py0[i] + py1[i] + py2[i] + py3[i] + py4[i]
            pz0[i] = pz0[i] + pz1[i] + pz2[i] + pz3[i] + pz4[i]
            px1[i] = px1[i] + 2*px2[i] + 3*px3[i] + 4*px4[i]
            py1[i] = py1[i] + 2*py2[i] + 3*py3[i] + 4*py4[i]
            pz1[i] = pz1[i] + 2*pz2[i] + 3*pz3[i] + 4*pz4[i]
            px2[i] = px2[i] + 3*px3[i] + 6*px4[i]
            py2[i] = py2[i] + 3*py3[i] + 6*py4[i]
            pz2[i] = pz2[i] + 3*pz3[i] + 6*pz4[i]
            px3[i] = px3[i] + 4*px4[i]
            py3[i] = py3[i] + 4*py4[i]
            pz3[i] = pz3[i] + 4*pz4[i]
        #
        #   this code calculates the forces on each molecule and the
        #   internal energy per molecule at time t.
        #
        esum = 0.0
        tsum = 0.0
        ptf = np.zeros((3,3)) # ptf accumulates the f*r part of pressure tensor
        ptp = np.zeros((3,3)) # ptp accumulates the momentum part of pressure tensor
        fx, fy, fz = np.zeros(nmol), np.zeros(nmol), np.zeros(nmol)  
        #
        # fx, fy, fz will accumulate x,y,z components of force for each molecule
        #
        for i in range(nmol-1):
            for j in range(i+1,nmol):
                xij = x0[i] - x0[j]
                yij = y0[i] - y0[j]
                zij = z0[i] - z0[j]
                xij = xij - boxlen*round(xij/boxlen)
                yij = yij - boxlen*round(yij/boxlen)
                zij = zij - boxlen*round(zij/boxlen)
                rijsqr = xij*xij + yij*yij + zij*zij
                if (rijsqr < cutf2):
                    r6 = 1/rijsqr**3
                    r12 = r6*r6
                    fij = (r12 + r12 - r6)/rijsqr
                    fxij = xij*fij
                    fyij = yij*fij
                    fzij = zij*fij
                    fx[i] = fx[i] + fxij
                    fy[i] = fy[i] + fyij
                    fz[i] = fz[i] + fzij
                    fx[j] = fx[j] - fxij
                    fy[j] = fy[j] - fyij
                    fz[j] = fz[j] - fzij
                    esum = esum + r12 - r6
                    ptf[0,0] = ptf[0,0] + xij*fxij
                    ptf[0,1] = ptf[0,1] + xij*fyij
                    ptf[0,2] = ptf[0,2] + xij*fzij
                    ptf[1,0] = ptf[1,0] + yij*fxij
                    ptf[1,1] = ptf[1,1] + yij*fyij
                    ptf[1,2] = ptf[1,2] + yij*fzij
                    ptf[2,0] = ptf[2,0] + zij*fxij
                    ptf[2,1] = ptf[2,1] + zij*fyij
                    ptf[2,2] = ptf[2,2] + zij*fzij
                # end of if statement
            # end of for j
        # end of for i
        esum = 4.0*esum
        epot = esum/nmol
        anum, aden, alpha = 0.0, 0.0, 0.0
        for i in range(nmol):
            fx[i] = 24.0*fx[i]
            fy[i] = 24.0*fy[i]
            fz[i] = 24.0*fz[i]
            aden = aden + px0[i]*px0[i] + py0[i]*py0[i] + pz0[i]*pz0[i]
            anum = anum + fx[i]*px0[i] + fy[i]*py0[i] + fz[i]*pz0[i]
        alpha = anum/aden
        if (nvt == 0):
            alpha = 0
        ekin = 0.5*aden/nmol
        etot = ekin + epot
        traj_ke[istep-1]=ekin
        traj_pe[istep-1]=epot
        traj_tot[istep-1]=etot
        #
        # this code performs the corrector cycle
        #
        for i in range(nmol):
            xcorr = x1[i] - px0[i]*delt
            ycorr = y1[i] - py0[i]*delt
            zcorr = z1[i] - pz0[i]*delt
            pxcorr = px1[i] - (fx[i] - alpha*px0[i])*delt
            pycorr = py1[i] - (fy[i] - alpha*py0[i])*delt
            pzcorr = pz1[i] - (fz[i] - alpha*pz0[i])*delt
            x0[i] = x0[i] - xcorr*f01
            y0[i] = y0[i] - ycorr*f01
            z0[i] = z0[i] - zcorr*f01
            px0[i] = px0[i] - pxcorr*f01
            py0[i] = py0[i] - pycorr*f01
            pz0[i] = pz0[i] - pzcorr*f01
            x1[i] = x1[i] - xcorr
            y1[i] = y1[i] - ycorr
            z1[i] = z1[i] - zcorr
            px1[i] = px1[i] - pxcorr
            py1[i] = py1[i] - pycorr
            pz1[i] = pz1[i] - pzcorr
            x2[i] = x2[i] - xcorr*f21
            y2[i] = y2[i] - ycorr*f21
            z2[i] = z2[i] - zcorr*f21
            px2[i] = px2[i] - pxcorr*f21
            py2[i] = py2[i] - pycorr*f21
            pz2[i] = pz2[i] - pzcorr*f21
            x3[i] = x3[i] - xcorr*f31
            y3[i] = y3[i] - ycorr*f31
            z3[i] = z3[i] - zcorr*f31
            px3[i] = px3[i] - pxcorr*f31
            py3[i] = py3[i] - pycorr*f31
            pz3[i] = pz3[i] - pzcorr*f31
            x4[i] = x4[i] - xcorr*f41
            y4[i] = y4[i] - ycorr*f41
            z4[i] = z4[i] - zcorr*f41
            px4[i] = px4[i] - pxcorr*f41
            py4[i] = py4[i] - pycorr*f41
            pz4[i] = pz4[i] - pzcorr*f41
            #
            # apply periodic boundary conditions.
            #
            rx = x0[i] - halfl
            ry = y0[i] - halfl
            rz = z0[i] - halfl
            #
            rx = rx - boxlen*round(rx/boxlen)
            ry = ry - boxlen*round(ry/boxlen)
            rz = rz - boxlen*round(rz/boxlen)
            #
            x0[i] = rx + halfl
            y0[i] = ry + halfl
            z0[i] = rz + halfl
            #
            #   put momentum part of pressure tensor in.
            #
            ptp[0,0] = ptp[0,0] + px0[i]*px0[i]
            ptp[0,1] = ptp[0,1] + px0[i]*py0[i]
            ptp[0,2] = ptp[0,2] + px0[i]*pz0[i]
            ptp[1,0] = ptp[1,0] + py0[i]*px0[i]
            ptp[1,1] = ptp[1,1] + py0[i]*py0[i]
            ptp[1,2] = ptp[1,2] + py0[i]*pz0[i]
            ptp[2,0] = ptp[2,0] + pz0[i]*px0[i]
            ptp[2,1] = ptp[2,1] + pz0[i]*py0[i]
            ptp[2,2] = ptp[2,2] + pz0[i]*pz0[i]
        #
        # divide ptp, ptf by voluime to obtain pressure tensor elements
        #
        for ix in range(3):
            for iy in range(3):
                ptp[ix,iy] = ptp[ix,iy]/vol
                ptf[ix,iy] = ptf[ix,iy]/vol
        eksum = eksum + ekin
        ek2sum = ek2sum + ekin*ekin
        epsum = epsum + epot
        ep2sum = ep2sum + epot*epot
        pressure = (ptp[0,0]+ptp[1,1]+ptp[2,2]+ptf[0,0]+ptf[1,1]+ptf[2,2])/3.0
        psum = psum + pressure
        p2sum = p2sum + pressure*pressure
        #
        # accumulate block average
        #
        j = istep - (istep/nblock)*nblock
        # 
        # if istep is a multiple of nblock, j = 0
        if (j == 0):
            iblk = istep/nblock-1
            ekmean = eksum/nblock
            blkavr[iblk,0] = ekmean
            print ekmean, ek2sum/nblock-ekmean*ekmean
            blkavr[iblk,1] = math.sqrt(ek2sum/nblock - ekmean*ekmean)
            epmean = epsum/nblock
            blkavr[iblk,2] = epmean
            blkavr[iblk,3] = math.sqrt(ep2sum/nblock - epmean*epmean)
            pmean = psum/nblock
            blkavr[iblk,4] = pmean
            blkavr[iblk,5] = math.sqrt(p2sum/nblock - pmean*pmean)
            blkavr[iblk,6] = time
            #
            #   zero out the accumulator variables before going to another block.
            #
            eksum = 0.0
            ek2sum = 0.0
            epsum = 0.0
            ep2sum = 0.0
            psum = 0.0
            p2sum = 0.0
            #
            #   rescale momenta and velocities
            #
            if (nvt == 1):
                mdtemp = 2.0*ekin/3.0
                scaling_factor = math.sqrt(tstar/mdtemp)
                for i in range(nmol):
                    px0[i] = scaling_factor*px0[i]
                    py0[i] = scaling_factor*py0[i]
                    pz0[i] = scaling_factor*pz0[i]
                    x1[i] = px0[i]*delt
                    y1[i] = py0[i]*delt
                    z1[i] = pz0[i]*delt
                    # end for i
                # end if j=0
            # end if nvt
        # end for istep
    return time

def read_prior_configuration(nmol,x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4,
          px0, py0, pz0, px1, py1, pz1, px2, py2, pz2, px3, py3, pz3, px4, py4, pz4):
    prior_configuration = open('mdeq.con','r')
    for i in range(nmol):
        numbers = prior_configuration.readline()
        x = numbers.split()
        x0[i] = string.atof(x[0])
        y0[i] = string.atof(x[1])
        z0[i] = string.atof(x[2])
        numbers = prior_configuration.readline()
        x = numbers.split()
        x1[i] = string.atof(x[0])
        y1[i] = string.atof(x[1])
        z1[i] = string.atof(x[2])
        numbers = prior_configuration.readline()
        x = numbers.split()
        x2[i] = string.atof(x[0])
        y2[i] = string.atof(x[1])
        z2[i] = string.atof(x[2])
        numbers = prior_configuration.readline()
        x = numbers.split()
        x3[i] = string.atof(x[0])
        y3[i] = string.atof(x[1])
        z3[i] = string.atof(x[2])
        numbers = prior_configuration.readline()
        x = numbers.split()
        x4[i] = string.atof(x[0])
        y4[i] = string.atof(x[1])
        z4[i] = string.atof(x[2])
        numbers = prior_configuration.readline()
        x = numbers.split()
        px0[i] = string.atof(x[0])
        py0[i] = string.atof(x[1])
        pz0[i] = string.atof(x[2])
        numbers = prior_configuration.readline()
        x = numbers.split()
        px1[i] = string.atof(x[0])
        py1[i] = string.atof(x[1])
        pz1[i] = string.atof(x[2])
        numbers = prior_configuration.readline()
        x = numbers.split()
        px2[i] = string.atof(x[0])
        py2[i] = string.atof(x[1])
        pz2[i] = string.atof(x[2])
        numbers = prior_configuration.readline()
        x = numbers.split()
        px3[i] = string.atof(x[0])
        py3[i] = string.atof(x[1])
        pz3[i] = string.atof(x[2])
        numbers = prior_configuration.readline()
        x = numbers.split()
        px4[i] = string.atof(x[0])
        py4[i] = string.atof(x[1])
        pz4[i] = string.atof(x[2])
    # end for i
    prior_configuration.close()
    return

def read_mdeq_dat():
    data_in = open('mdeq.dat','r') 
    numbers = data_in.readline() # read first line of mdeq.dat
    x = numbers.split()
    print x
    irun = string.atoi(x[0])
    nmol = string.atoi(x[1])
    nstep = string.atoi(x[2])
    nblock = string.atoi(x[3])
    nvt = string.atoi(x[4])
    temp = string.atof(x[5])
    rho = string.atof(x[6])
    eps = string.atof(x[7])
    sig = string.atof(x[8])
    mass = string.atof(x[9])
    numbers = data_in.readline() # read second line of mdeq.dat
    x = numbers.split() 
    cutoff = string.atof(x[0]) 
    delt = string.atof(x[1])
    time = string.atof(x[2]) 
    data_in.close()
    return irun,nmol,nstep,nblock,nvt,temp,rho,eps,sig,mass,cutoff,delt,time

def write_final_configuration(nmol,x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4,
          px0, py0, pz0, px1, py1, pz1, px2, py2, pz2, px3, py3, pz3, px4, py4, pz4):
    final_configuration = open('mdeq.con','w')
    for i in range(nmol):
        final_configuration.write('{}'.format(x0[i])+"   "+'{}'.format(y0[i])+\
                                  "   "+'{}'.format(z0[i])+"\n")
        final_configuration.write('{}'.format(x1[i])+"   "+'{}'.format(y1[i])+\
                                  "   "+'{}'.format(z1[i])+"\n")
        final_configuration.write('{}'.format(x2[i])+"   "+'{}'.format(y2[i])+\
                                  "   "+'{}'.format(z2[i])+"\n")
        final_configuration.write('{}'.format(x3[i])+"   "+'{}'.format(y3[i])+\
                                  "   "+'{}'.format(z3[i])+"\n")
        final_configuration.write('{}'.format(x4[i])+"   "+'{}'.format(y4[i])+\
                                  "   "+'{}'.format(z4[i])+"\n")
        final_configuration.write('{}'.format(px0[i])+"   "+'{}'.format(py0[i])+\
                                  "   "+'{}'.format(pz0[i])+"\n")
        final_configuration.write('{}'.format(px1[i])+"   "+'{}'.format(py1[i])+\
                                  "   "+'{}'.format(pz1[i])+"\n")
        final_configuration.write('{}'.format(px2[i])+"   "+'{}'.format(py2[i])+\
                                  "   "+'{}'.format(pz2[i])+"\n")
        final_configuration.write('{}'.format(px3[i])+"   "+'{}'.format(py3[i])+\
                                  "   "+'{}'.format(pz3[i])+"\n")
        final_configuration.write('{}'.format(px4[i])+"   "+'{}'.format(py4[i])+\
                                  "   "+'{}'.format(pz4[i])+"\n")
    # end for i
    final_configuration.close()
    return





    
    

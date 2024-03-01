import rebound
import numpy as np
import csv
import h5py
import time
filename = "NBody_MCMC_Posteriors.hdf5"
f = h5py.File(filename, "r")
system = 'Kepler-23' ##

start_time = time.process_time()
##Retrieve epoch and star mass data
file = open("Epochs.txt", 'r') 
KOI = "1236.01" ##
Epoch = 0
mstar = 1
lines = []
for line in file:
    l=line.split()
    if (str(l[0])) == KOI:
        Epoch = float(l[3])
        mstar = float(l[11])
        break
file.close()

file = open(system+'_output.txt', 'w') ##output file
file.write("System\tIndex\tSimulation\tMass Ecc alternating \tIn. Semi-major and Final Semi-major and percentage difference alternating \tUnstable \n")
num = 1 ##number of simulations
length = len(f[system]['DefaultPriors']['PosteriorSample'])
partition = length//num
for index in range(num):
    file.write(system+"\t")
    print("Simulation " + str(index))
    ##first sim
    sim = rebound.Simulation()
    sim.add(m=mstar)    
    n = np.random.randint(index*partition, (index+1)*partition)
    file.write(str(n) + ", " + str(n+1) + "\t")
    file.write(str(index)+"\t\t")
    smas = np.zeros(8)
    data = list(f[system]['DefaultPriors']['PosteriorSample'][n])
    nplanets=len(data)//5
    nsystems = 8//nplanets
    pmin = data[1]/365.25*2*np.pi

    if nplanets<3:    
        npersim=2
        for i in range(nplanets):
            m = data[i*5]*mstar
            file.write(str(m)+"\t")
            P = data[i*5+1]/365.*2*np.pi
            if P< pmin:
               pmin = P
            h = data[i*5+2]
            k = data[i*5+3]
            l = np.mod((Epoch-data[i*5+4])/data[i*5+1]*2*np.pi, 2*np.pi)
            inc = np.random.random() * 0.005/180*np.pi
            Omega = np.random.random()*2*np.pi
            sim.add(m = m, P = P, h = k, k = h, ix = 2*np.sin(inc/2)*np.cos(Omega), iy = 2*np.sin(inc/2)*np.sin(Omega), l = l)
            o = sim.particles[-1].orbit(primary=sim.particles[0])
            smas[i] = o.a
            file.write(str(o.e)+"\t")
        for i in range(2-nplanets):
            h = np.random.random()*data[i*5+2]
            k = np.random.random()*data[i*5+3]
            l = np.mod((Epoch-data[i*5+4])/data[i*5+1]*2*np.pi, 2*np.pi)
            inc = np.random.random() *0.005/180*np.pi
            Omega = np.random.random()*2*np.pi
            sim.add(m = 0, P = P, h = k, k = h, ix = 2*np.sin(inc/2)*np.cos(Omega), iy = 2*np.sin(inc/2)*np.sin(Omega), l = l)        
        sim.move_to_com()
        
        ##second sim
        data = list(f[system]['DefaultPriors']['PosteriorSample'][n+1])
        sim2 = rebound.Simulation()
        sim2.add(m=mstar)
        for i in range(nplanets):
            m = data[i*5]*mstar
            file.write(str(m)+"\t")
            P = data[i*5+1]/365.*2*np.pi
            if P< pmin:
               pmin = P
            h = data[i*5+2]
            k = data[i*5+3]
            l = np.mod((Epoch-data[i*5+4])/data[i*5+1]*2*np.pi, 2*np.pi)
            inc = np.random.random() * 0.005/180*np.pi
            Omega = np.random.random()*2*np.pi
            sim2.add(m = m, P = P, h = k, k = h, ix = 2*np.sin(inc/2)*np.cos(Omega), iy = 2*np.sin(inc/2)*np.sin(Omega), l = l)
            o = sim2.particles[-1].orbit(primary=sim2.particles[0])
            smas[i+2] = o.a
            file.write(str(o.e)+"\t")
        for i in range(2-nplanets):
            h = np.random.random()*data[i*5+2]
            k = np.random.random()*data[i*5+3]
            l = np.mod((Epoch-data[i*5+4])/data[i*5+1]*2*np.pi, 2*np.pi)
            inc = np.random.random() *0.005/180*np.pi
            Omega = np.random.random()*2*np.pi
            sim2.add(m = 0, P = P, h = k, k = h, ix = 2*np.sin(inc/2)*np.cos(Omega), iy = 2*np.sin(inc/2)*np.sin(Omega), l = l)   
        sim2.move_to_com()
        
        ##combine simulations
        for particle in sim2.particles:
            sim.add(particle)

        #third sim
        data = list(f[system]['DefaultPriors']['PosteriorSample'][n+2])
        sim3 = rebound.Simulation()
        sim3.add(m=mstar)
        for i in range(nplanets):
            m = data[i*5]*mstar
            file.write(str(m)+"\t")
            P = data[i*5+1]/365.*2*np.pi
            if P< pmin:
               pmin = P
            h = data[i*5+2]
            k = data[i*5+3]
            l = np.mod((Epoch-data[i*5+4])/data[i*5+1]*2*np.pi, 2*np.pi)
            inc = np.random.random() * 0.005/180*np.pi
            Omega = np.random.random()*2*np.pi
            sim3.add(m = m, P = P, h = k, k = h, ix = 2*np.sin(inc/2)*np.cos(Omega), iy = 2*np.sin(inc/2)*np.sin(Omega), l = l)
            o = sim3.particles[-1].orbit(primary=sim3.particles[0])
            smas[i+4] = o.a
            file.write(str(o.e)+"\t")
        for i in range(2-nplanets):
            h = np.random.random()*data[i*5+2]
            k = np.random.random()*data[i*5+3]
            l = np.mod((Epoch-data[i*5+4])/data[i*5+1]*2*np.pi, 2*np.pi)
            inc = np.random.random() *0.005/180*np.pi
            Omega = np.random.random()*2*np.pi
            sim3.add(m = 0, P = P, h = k, k = h, ix = 2*np.sin(inc/2)*np.cos(Omega), iy = 2*np.sin(inc/2)*np.sin(Omega), l = l)
        sim3.move_to_com()
        ##combine simulations
        for particle in sim3.particles:
            sim.add(particle)

        #fourth sim
        data = list(f[system]['DefaultPriors']['PosteriorSample'][n+4])
        sim4 = rebound.Simulation()
        sim4.add(m=mstar)
        for i in range(nplanets):
            m = data[i*5]*mstar
            file.write(str(m)+"\t")
            P = data[i*5+1]/365.*2*np.pi
            if P< pmin:
               pmin = P
            h = data[i*5+2]
            k = data[i*5+3]
            l = np.mod((Epoch-data[i*5+4])/data[i*5+1]*2*np.pi, 2*np.pi)
            inc = np.random.random() * 0.005/180*np.pi
            Omega = np.random.random()*2*np.pi
            sim4.add(m = m, P = P, h = k, k = h, ix = 2*np.sin(inc/2)*np.cos(Omega), iy = 2*np.sin(inc/2)*np.sin(Omega), l = l)
            o = sim4.particles[-1].orbit(primary=sim4.particles[0])
            smas[i+6] = o.a
            file.write(str(o.e)+"\t")
        for i in range(2-nplanets):
            h = np.random.random()*data[i*5+2]
            k = np.random.random()*data[i*5+3]
            l = np.mod((Epoch-data[i*5+4])/data[i*5+1]*2*np.pi, 2*np.pi)
            inc = np.random.random() *0.005/180*np.pi
            Omega = np.random.random()*2*np.pi
            sim4.add(m = 0, P = P, h = k, k = h, ix = 2*np.sin(inc/2)*np.cos(Omega), iy = 2*np.sin(inc/2)*np.sin(Omega), l = l)
           
        sim4.move_to_com()
        ##combine simulations
        for particle in sim4.particles:
            sim.add(particle)

    elif nplanets<5:
        npersim=4
        for i in range(nplanets):
            m = data[i*5]*mstar
            file.write(str(m)+"\t")
            P = data[i*5+1]/365.*2*np.pi
            if P< pmin:
               pmin = P
            h = data[i*5+2]
            k = data[i*5+3]
            l = np.mod((Epoch-data[i*5+4])/data[i*5+1]*2*np.pi, 2*np.pi)
            inc = np.random.random() * 0.005/180*np.pi
            Omega = np.random.random()*2*np.pi
            sim.add(m = m, P = P, h = k, k = h, ix = 2*np.sin(inc/2)*np.cos(Omega), iy = 2*np.sin(inc/2)*np.sin(Omega), l = l)
            o = sim.particles[-1].orbit(primary=sim.particles[0])
            smas[i] = o.a
            file.write(str(o.e)+"\t")
        for i in range(4-nplanets):
            h = np.random.random()*data[i*5+2]
            k = np.random.random()*data[i*5+3]
            l = np.mod((Epoch-data[i*5+4])/data[i*5+1]*2*np.pi, 2*np.pi)
            inc = np.random.random() *0.005/180*np.pi
            Omega = np.random.random()*2*np.pi
            sim.add(m = 0, P = P, h = k, k = h, ix = 2*np.sin(inc/2)*np.cos(Omega), iy = 2*np.sin(inc/2)*np.sin(Omega), l = l)
        sim.move_to_com()

        ##second sim
        data = list(f[system]['DefaultPriors']['PosteriorSample'][n+1])
        sim2 = rebound.Simulation()
        sim2.add(m=mstar)
        for i in range(nplanets):
            m = data[i*5]*mstar
            file.write(str(m)+"\t")
            P = data[i*5+1]/365.*2*np.pi
            if P< pmin:
               pmin = P
            h = data[i*5+2]
            k = data[i*5+3]
            l = np.mod((Epoch-data[i*5+4])/data[i*5+1]*2*np.pi, 2*np.pi)
            inc = np.random.random() * 0.005/180*np.pi
            Omega = np.random.random()*2*np.pi
            sim2.add(m = m, P = P, h = k, k = h, ix = 2*np.sin(inc/2)*np.cos(Omega), iy = 2*np.sin(inc/2)*np.sin(Omega), l = l)
            o = sim2.particles[-1].orbit(primary=sim2.particles[0])
            smas[i+4] = o.a
            file.write(str(o.e)+"\t")
        for i in range(4-nplanets):
            h = np.random.random()*data[i*5+2]
            k = np.random.random()*data[i*5+3]
            l = np.mod((Epoch-data[i*5+4])/data[i*5+1]*2*np.pi, 2*np.pi)
            inc = np.random.random() *0.005/180*np.pi
            Omega = np.random.random()*2*np.pi
            sim2.add(m = 0, P = P, h = k, k = h, ix = 2*np.sin(inc/2)*np.cos(Omega), iy = 2*np.sin(inc/2)*np.sin(Omega), l = l)
           
        sim2.move_to_com()
        ##combine simulations
        for particle in sim2.particles:
            sim.add(particle)
            
    else:
        npersim=8
        for i in range(nplanets):
            m = data[i*5]*mstar
            file.write(str(m)+"\t")
            P = data[i*5+1]/365.*2*np.pi
            if P< pmin:
               pmin = P
            h = data[i*5+2]
            k = data[i*5+3]
            l = np.mod((Epoch-data[i*5+4])/data[i*5+1]*2*np.pi, 2*np.pi)
            inc = np.random.random() *0.005/180*np.pi
            Omega = np.random.random()*2*np.pi
            sim.add(m = m, P = P, h = k, k = h, ix = 2*np.sin(inc/2)*np.cos(Omega), iy = 2*np.sin(inc/2)*np.sin(Omega), l = l)
            o = sim.particles[-1].orbit(primary=sim.particles[0])
            smas[i] = o.a
            file.write(str(o.e)+"\t")
        for i in range(8-nplanets):
            h = np.random.random()*data[i*5+2]
            k = np.random.random()*data[i*5+3]
            l = np.mod((Epoch-data[i*5+4])/data[i*5+1]*2*np.pi, 2*np.pi)
            inc = np.random.random() *0.005/180*np.pi
            Omega = np.random.random()*2*np.pi
            sim.add(m = 0, P = P, h = k, k = h, ix = 2*np.sin(inc/2)*np.cos(Omega), iy = 2*np.sin(inc/2)*np.sin(Omega), l =l)
        sim.move_to_com()
    
    #print(sim.particles[:])
    ##integrate
    sim.dt = 0.09 *pmin
    sim.t = Epoch/365.25*2*np.pi
    sim.ri_whfast512.N_systems = nsystems
    sim.G = 1
    sim.integrator = "whfast512"
    int_time = pmin*10000000 #interval of time
    file.write("\n")
    for t in range (10): # repetition
        sim.integrate(sim.t+int_time, exact_finish_time=0)
        p = sim.particles[1:1+npersim]
        for sys in range(1,nsystems):
            for planet in range(npersim):
                p.append(sim.particles[sys*(npersim+1)+1+planet])
        unstable = np.zeros(nsystems)
        #print(len(p))
        #print(smas)
        for ind, particle in enumerate(p): #check for instability
            sys =  ind//npersim
            print(sys*(npersim+1))
            o = particle.orbit(primary=sim.particles[sys*(npersim+1)])
            if smas[ind]!=0 and (o.a < 0.9*smas[ind] or o.a > 1.1*smas[ind] or o.a <0 or np.isnan(o.a)):
                unstable[sys] = 1
        #if (all(element == 1 for element in unstable)):
        #   break
    ##final output
    p = sim.particles[1:1+npersim]
    for sys in range(1,nsystems):
        for planet in range(npersim):
            p.append(sim.particles[sys*(npersim+1)+1+planet])
    #print(p)
    for ind, particle in enumerate(p): 
        if smas[ind]!=0:
            sys =  ind//npersim
            o = particle.orbit(primary=sim.particles[sys*(npersim+1)])
            file.write(str(smas[ind])+"\t")
            file.write(str(o.a)+"\t")
            per_chan = abs(o.a-smas[ind])/smas[ind]*100
            file.write(str(per_chan)+"\t")
            file.write(str(sim.t)+"\t")
    file.write("\t"+str(unstable)+"\n")
    print("Stability: "+str(unstable)+"\n")
file.close()

end_time = time.process_time()
timetaken = end_time - start_time
print("Time taken: "+str(timetaken))

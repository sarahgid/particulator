class Particle {
  // To calculate trajectories by some other method, overload this class.
  // Eventually, it might be better to make this class overload a more generic Particle.

  ParticleExpt expt;
  ROMSRun run;
  HashMap current, initial, prev;
  boolean surfaceTrapped = false; // retained for back-compatibility; same as sigmaTrapped with trapLevel = 0
  boolean zTrapped = false, sigmaTrapped = false;
  float trapLevel = 0;
  boolean diffusive = true;
  boolean midpointStepping = true;
  int step = 0; // number of steps taken
  boolean active = true;
  float dt;
  
  Particle(float x0, float y0, float z0, float t0, ParticleExpt expt) {
    this.expt = expt;
    run = expt.run;
    dt = expt.dt; // a change from previous behavior!
    current = new HashMap();
    current.put("x",x0);
    current.put("y",y0);
    current.put("z",z0);
    current.put("t",t0);
    current.put("wdiff",0); // if the diffusion code is never called, wdiff and dksdz will retain 0
    current.put("dkskz",0); 
    interpEverything();
    initial = (HashMap)current.clone();
    prev = (HashMap)current.clone();
  }
  
  void trapToZLevel(float z) {
    zTrapped = true;
    sigmaTrapped = false;
    surfaceTrapped = false;
    trapLevel = z;
  }
  
  void trapToSigmaLevel(float cs) {
    sigmaTrapped = true;
    zTrapped = false;
    trapLevel = constrain(cs, -1, 0);
    surfaceTrapped = (cs==0);
  }
  
  float x() {return (Float)current.get("x");}
  float y() {return (Float)current.get("y");}
  float z() {return (Float)current.get("z");}
  float cs() {return (Float)current.get("cs");}
  float t() {return (Float)current.get("t");}
  float u() {return (Float)current.get("u");}
  float v() {return (Float)current.get("v");}
  float w() {return (Float)current.get("w");}
  float wdiff() {return (Float)current.get("wdiff");}
  float dksdz() {return (Float)current.get("dksdz");}
  
  
  void interpEverything() { // define velocities and other derived quantities that go with the current (x,y,z,t)
    // synchronize cs and z, make sure both are good values
    if (surfaceTrapped) {
      current.put("cs",0.0);
      current.put("z",run.interpZeta(t(), y(), x()));
    } else if (sigmaTrapped) {
      current.put("cs",trapLevel);
      current.put("z",run.interpZeta(t(), y(), x()));
    } else if (zTrapped) {
      current.put("z",trapLevel);
      float cs = run.z2cs(t(), z(), y(), x());
      if (cs != constrain(cs, -1, 0)) {
        cs = constrain(cs, -1, 0);
        current.put("z",run.cs2z(t(), cs, y(), x()));
      }
      current.put("cs",cs);      
    } else { // regular 3-d motion
      float cs = run.z2cs(t(), z(), y(), x());
      if (cs != constrain(cs, -1, 0)) {
        cs = constrain(cs, -1, 0);
        current.put("z",run.cs2z(t(), cs, y(), x()));
      }
      current.put("cs",cs);
    }
    // alternate coordinates
    current.put("i",findIndexBefore(run.Xu, x()));
    current.put("j",findIndexBefore(run.Yv, y()));
    current.put("j",findIndexBefore(run.Csw, cs()));
    current.put("lon",run.meters2lon(x()));
    current.put("lat",run.meters2lat(y()));    
    // ancillary variables
    current.put("H",run.interpH(y(), x()));
    current.put("zeta",run.interpZeta(t(), y(), x()));
    current.put("mask",run.interpMask(t(), y(), x()));
    // velocity and diffusion
    current.put("u",run.interpU(t(), cs(), y(), x()));
    current.put("v",run.interpV(t(), cs(), y(), x()));
    current.put("w",run.interpW(t(), cs(), y(), x()));
    if (diffusive) {
      current.put("dksdz",diffusionGradient());
      current.put("wdiff",diffusionVelocity());
    }
    // additional tracers
    for (int i=0; i<expt.tracerNames.length; i++) {
      float cs = cs(); // by default, interpolate the tracer at the actual position of the particle
      if (expt.tracerInterpMode[i].equals("cs")) { // interpolate the tracer at the (t,y,x) of the particle and a given cs
        cs = expt.tracerInterpDepth[i]; 
      } else if (expt.tracerInterpMode[i].equals("z")) { // likewise but for a given z
        cs = run.z2cs(t(), expt.tracerInterpDepth[i], y(), x());
      }
      current.put(expt.tracerSaveName(i), run.interpTracer(expt.tracerNames[i], t(), cs, y(), x()));        
    }
  }
  
  
  float diffusionGradient() {
    float dz = 1; // ****half-span to take d(Ks)/dz over****
    float cstop = constrain(run.z2cs(t(), z()+dz, y(), x()), -1, 0);
    float csbot = constrain(run.z2cs(t(), z()-dz, y(), x()), -1, 0);
    float kstop = run.interpKs(t(), cstop, y(), x());
    float ksbot = run.interpKs(t(), csbot, y(), x());
    float H = (Float)current.get("H");
    float zeta = (Float)current.get("zeta");
    return (kstop-ksbot) / (cstop - csbot) / (zeta+H);
  }
  
  
  float diffusionVelocity() {
    // Ks at the appropriate z
    float cs0 = constrain(run.z2cs(t(), z() + 0.5*dksdz()*abs(dt), y(), x()), -1, 0); // Batchelder et al 2000
    float ks0 = run.interpKs(t(), cs0, y(), x());
    // diffusion velocity
    float randn1 = sqrt(-2*log(random(1)))*cos(TWO_PI*random(1)); // normally distributed random variable with std dev 1
    return sqrt(2*ks0/abs(dt)) * randn1;
  }
  
  
  void storePrevious() {
    prev = (HashMap)current.clone();
  }
  
  void restoreFromPrevious() { // everything except t
    float tcurr = t();
    current = (HashMap)prev.clone();
    current.put("t",tcurr);
  }
  
  
  boolean inDomain() {
    int i = (Integer)current.get("i");
    int j = (Integer)current.get("j");
    if (!isfinite(x()+y()+z()+u()+v()+w()+i+j)) return false;
    if (i != constrain(i,0,run.Iu-2) || j != constrain(j,0,run.Jv-2)) return false;
    return true;
  }
  
  
  void takeStep() {
    if (!inDomain()) {
      current.put("t",t()+dt);
      return;
    }
    
    float x = x();
    float y = y();
    float z = z();
    float t = t();
    
    if (!midpointStepping) { // simple Euler step
      storePrevious();
      x += u()*dt;
      y += v()*dt;
      z += (w()+wdiff()+dksdz())*dt;
      t += dt;
    } else { // midpoint method
      storePrevious();
      x += u()*0.5*dt;
      y += v()*0.5*dt;
      z += w()*0.5*dt;
      t += 0.5*dt;
      interpEverything(); // at the midpoint
      float umid = u();
      float vmid = v();
      float wmid = w();
      restoreFromPrevious();
      x += umid*dt; // take the full step with the velocities from the midpoint (but wdiff, dksdz as calculated at the start point)
      y += vmid*dt;
      z += (w()+wdiff()+dksdz())*dt;
      t += 0.5*dt;
    }
    
    current.put("x",x);
    current.put("y",y);
    current.put("z",z);
    current.put("t",t);

    interpEverything();   
    if (u()==0 && v()==0) restoreFromPrevious();    
    step++;
  }
  
}

class Particle {
  // To calculate trajectories by some other method, overload this class.
  // Eventually, it might be better to make this class overload a more generic Particle.

  ParticleExpt expt;
  ROMSRun run;
  boolean surfaceTrapped = false;
  boolean diffusive = true;
  boolean midpointStepping = true;
  float lon,lat,x,y,z,cs,t,u,v,w,wdiff,dksdz,mask,H,zeta; // current position and velocity
  int i,j,k; // current grid location (i..i+1, j..j+1, k..k+1)
  float lon0,lat0,x0,y0,z0,cs0,u0,v0,w0,wdiff0,dksdz0,mask0,H0,zeta0; // initial  
  float lon1,lat1,x1,y1,z1,cs1,u1,v1,w1,wdiff1,dksdz1,mask1,H1,zeta1; // previous
  int i0,j0,k0; // initial
  int i1,j1,k1; // previous
  int step = 0; // number of steps taken
  boolean active = true;
  float dt = 400; // 200-400 seems about right for Puget Sound
  
  Particle(float x0, float y0, float z0, float t0, ParticleExpt expt) {
    x = x0;
    y = y0;
    z = z0;
    t = t0;
    this.expt = expt;
    run = expt.run;
    wdiff = 0;
    dksdz = 0;
    interpEverything();
    // save initial values
    lon0 = lon;
    lat0 = lat;
    x0 = x;
    y0 = y;
    z0 = z;
    cs0 = cs;
    u0 = u;
    v0 = v;
    w0 = w;
    wdiff0 = wdiff;
    dksdz0 = dksdz;
    i0 = i;
    j0 = j;
    k0 = k;
    mask0 = mask;
    H0 = H;
    zeta0 = zeta;
  }
  
  void interpEverything() {
    // define velocities and other derived quantities that go with the current (x,y,z,t)
    if (surfaceTrapped) {
      cs = 0;
      z = run.interpZeta(t, y, x);      
    } else {
      cs = run.z2cs(t, z, y, x);
      if (cs != constrain(cs, -1, 0)) {
        cs = constrain(cs, -1, 0);
        z = run.cs2z(t, cs, y, x);
      }
    }
    u = run.interpU(t, cs, y, x);
    v = run.interpV(t, cs, y, x);
    w = run.interpW(t, cs, y, x);
    i = findIndexBefore(run.Xu, x);
    j = findIndexBefore(run.Yv, y);
    k = findIndexBefore(run.Csw, cs);
    lon = run.meters2lon(x);
    lat = run.meters2lat(y);
    H = run.interpH(y, x);
    zeta = run.interpZeta(t,y,x);
    mask = run.interpMask(t, y, x);
    if (diffusive) {
      dksdz = diffusionGradient();
      wdiff = diffusionVelocity();
    }
  }
  
  float diffusionGradient() {
    float dz = 1; // ****half-span to take d(Ks)/dz over****
    float cstop = constrain(run.z2cs(t, z+dz, y, x), -1, 0);
    float csbot = constrain(run.z2cs(t, z-dz, y, x), -1, 0);
    float kstop = run.interpKs(t, cstop, y, x);
    float ksbot = run.interpKs(t, csbot, y, x);
    return (kstop-ksbot) / (cstop - csbot) / (zeta+H);
  }
  
  float diffusionVelocity() {
    // Ks at the appropriate z
    float cs0 = constrain(run.z2cs(t, z + 0.5*dksdz*abs(dt), y, x), -1, 0); // Batchelder et al 2000
    float ks0 = run.interpKs(t, cs0, y, x);
    // diffusion velocity
    float randn1 = sqrt(-2*log(random(1)))*cos(TWO_PI*random(1)); // normally distributed random variable with std dev 1
    return sqrt(2*ks0/abs(dt)) * randn1;
  }
  
  void storePrevious() { // everything except t
    lon1 = lon;
    lat1 = lat;
    x1 = x;
    y1 = y;
    z1 = z;
    cs1 = cs;
    u1 = u;
    v1 = v;
    w1 = w;
    wdiff1 = wdiff;
    dksdz1 = dksdz;
    i1 = i;
    j1 = j;
    k1 = k;
    mask1 = mask;
    H1 = H;
    zeta1 = zeta;
  }
  
  void restoreFromPrevious() { // everything except t
    lon = lon1;
    lat = lat1;
    x = x1;
    y = y1;
    z = z1;
    cs = cs1;
    u = u1;
    v = v1;
    w = w1;
    wdiff = wdiff1;
    dksdz = dksdz1;
    i = i1;
    j = j1;
    k = k1;
    mask = mask1;
    H = H1;
    zeta = zeta1;
  }
  
  void takeStep() {
    // make sure it's still in the domain
    if (!isfinite(i+j+k+x+y+z+u+v+w+cs) || i!=constrain(i,0,run.Iu-2) || j!=constrain(j,0,run.Jv-2)) active = false;
    if (!active) return;
    
    if (!midpointStepping) { // simple Euler step
      storePrevious();
      x += u*dt;
      y += v*dt;
      z += (w+wdiff+dksdz)*dt;
      t += dt;
    } else { // midpoint method
      storePrevious();
      x += u*0.5*dt;
      y += v*0.5*dt;
      z += w*0.5*dt;
      t += 0.5*dt;
      interpEverything(); // at the midpoint
      float umid = u;
      float vmid = v;
      float wmid = w;
      restoreFromPrevious();
      x += umid*dt; // take the full step with the velocities from the midpoint (but wdiff, dksdz as calculated at the start point)
      y += vmid*dt;
      z += (w+wdiff+dksdz)*dt;
      t += 0.5*dt;
    }

    interpEverything();    
    if (u==0 && v==0) restoreFromPrevious();    
    step++;
  }
  
}

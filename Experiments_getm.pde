class GETM_2DExpt extends ParticleExpt {
  GETM_2DExpt() {}  
  String[] saveNames() {return new String[] {"x","y","t","H","zeta","u","v"};}  
  float[] saveValues(Particle P) {return new float[] {P.x, P.y, P.t, P.H, P.zeta, P.u, P.v};}
}


  
void getmMap2D(String filename, String outname) {
  ParticleExpt expt = new GETM_2DExpt();
  expt.run = new GETM_2DRun(filename);
  expt.ncname = outname;
  float startTime = expt.run.firstTime() + 44712*4;
  float endTime = startTime + 44712*4;
  float[] lon = new float[expt.run.I-2];
  for (int i=0; i<lon.length; i++) lon[i] = expt.run.lon[i+1];
  float[] lat = new float[expt.run.J-2];
  for (int j=0; j<lat.length; j++) lat[i] = expt.run.lat[j+1];
  expt.seedParticles(lon, lat, new float[] {0}, startTime, 1);
  for (int i=0; i<expt.particles.length; i++) {
    expt.particles[i].dt = 108;
    expt.particles[i].surfaceTrapped = true;
    expt.particles[i].diffusive = false;
    expt.particles[i].midpointStepping = false;
  }
  expt.saveInterval = 1;
  expt.calcToTime(endTime);
}


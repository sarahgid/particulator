void getmMap2D_set(String basedir, int firstRun, int lastRun) {
  float[] starts = new float[4]; // 4 start times over the tidal cycle
  for (int i=0; i<starts.length; i++) starts[i] = 44712*4 + (44712/starts.length) * i;
  float[] durations = {44712}; // just track for a single cycle  
  for (int ri=firstRun; ri<=lastRun; ri++) {
    String filename = basedir + "run_" + ri + "/model_output.2d.nc";
    int rel = 0;
    for (int di=0; di<durations.length; di++) {
      for (int si=0; si<starts.length; si++) {
        rel++;
        String outname = "run_" + ri + "/particles.run" + ri + ".release" + rel + ".nc";
        println(outname + "-------------------------------------------------------");
        getmMap2D(filename, basedir + outname, starts[si], durations[di]);
      }
    }
  }
}




void getmMap2D(String filename, String outname, float startTimeOffset, float duration) {
  ParticleExpt expt = new ParticleExpt();
  expt.saveNames = new String[] {"x","y","t","H","zeta","u","v","mask"};
  expt.run = new GETM_2DRun(filename);
  expt.ncname = outname;
  float startTime = expt.run.firstTime() + startTimeOffset;
  float endTime = startTime + duration;
  float[] lon = new float[expt.run.I-2]; for (int i=0; i<lon.length; i++) lon[i] = expt.run.lon[i+1];
  float[] lat = new float[expt.run.J-2]; for (int j=0; j<lat.length; j++) lat[j] = expt.run.lat[j+1];
  expt.preallocSteps = round(duration / 27 / 9);
  expt.seedParticles(lon, lat, new float[] {0}, startTime, 1);
  for (int i=0; i<expt.particles.length; i++) {
    expt.particles[i].dt = 207;
    expt.particles[i].surfaceTrapped = true;
    expt.particles[i].diffusive = false;
    expt.particles[i].midpointStepping = false;
  }
  expt.saveInterval = 9;
  expt.calcToTime(endTime);
}

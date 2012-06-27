class Integrator implements Runnable {

  ArrayList queue = new ArrayList();
  ParticleRelease release;
  NetcdfFileWriteable nc;
  float endTime = 0;
  
  Integrator(ParticleRelease release, NetcdfFileWriteable nc) {
    this.release = release;
    this.nc = nc;
  }
  
  void add(Particle P) {
    queue.add(P);
  }
  
  void run() {
    while (queue.size() > 0) {
      Particle P = (Particle)queue.get(queue.size()-1);
      calcOne(P);
      queue.remove(P);
    }
  }
  
  void calcOne(Particle P) {
    while (P.active && P.t() < endTime) { // ...integrate each particle to the end of the loaded frame (or tend, whichever comes first)
      P.takeStep();
      /* autoSave doesn't work with multithreading, but if it did:
      if (release.autoSave && P.step % release.saveInterval == 0) {
        release.savePosition(nc, P.step / release.saveInterval, P.id, P); // save to netcdf file
      }
      */
    }
  }

}

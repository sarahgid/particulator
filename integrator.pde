class Integrator implements Runnable {

  ArrayList queue = new ArrayList();
  ParticleRelease parent;
  NetcdfFileWriteable nc;
  float endTime = 0;
  
  Integrator(ParticleRelease parent, NetcdfFileWriteable nc) {
    this.parent = parent;
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
      if (parent.autoSave && P.step % parent.saveInterval == 0) {
        parent.savePosition(nc, P.step / parent.saveInterval, P.id, P); // save to netcdf file
      }
    }
  }

}

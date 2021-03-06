import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.ArrayList;


class ParticleRelease {

  ROMSRun run;
  Particle[] particles = new Particle[0]; // a flat list of all particles: mandatory. This can be assembled however you want.
  Particle[][][][][] particlesRNKJI; // this is an optional, more organized, alternate indexing of the particles
                                     // (reps x release time x release depth x release lat x release lon).
                                     // it's populated by some versions of seedParticles but not used anywhere in the
                                     // basic ParticleRelease class: it's useful in specialized cases like ReturnMap.
  String ncname = "myrelease.nc"; // this and the other default values below are just placeholders. Don't set them here: that's what ParticleRelease is for.
  int saveInterval = 1;
  int preallocSteps = 1;
  boolean saveFirstLastOnly = false;
  boolean autoSave = true; // the normal way of saving, when multithreading is off. Ignore this flag until the multithreaded saving options are rethought!
  String[] saveNames = {"lon","lat","z","cs","t","H","mask"};
  // the following variables are copied from the ParticleRelease into each Particle in the Particle constructor
  float dt = 400;
  String[] tracerNames = new String[0];
  String[] tracerInterpMode = new String[0];
  float[] tracerInterpDepth = new float[0];
  
  boolean multithread = false;
  int nThreads = 6;
  Integrator[] integrators;
  int saveStep;
  
  
  ParticleRelease() {}
  
  void linkToRun(String basename, int ncn0, int ncn1) {
    run = new ROMSRun(basename, ncn0, ncn1);
  }
    
  void addTracer(String name) {addTracer(name,"3D",0);}
  void addTracer(String name, String interpMode, float interpDepth) {
    tracerNames = (String[]) append(tracerNames, name);
    tracerInterpMode = (String[]) append(tracerInterpMode, interpMode);
    tracerInterpDepth = (float[]) append(tracerInterpDepth, interpDepth);
    int ii = tracerNames.length-1;
    boolean found = false;
    for (int i=0; i<saveNames.length && !found; i++) found = saveNames[i].equals(tracerSaveName(ii));
    if (!found) saveNames = (String[]) append(saveNames, tracerSaveName(ii));
    run.tracerNames = tracerNames; // tell the ROMSRun object that it needs to load more fields when it reads a new ROMS file
    run.reload();
  }
  
  String tracerSaveName(int i) {
    String name = tracerNames[i];
    if (tracerInterpMode[i] !="3D") name += "_" + tracerInterpMode[i] + "_" + tracerInterpDepth[i];
    name = strrep(name,'.',"");
    return name;
  }



  void seedParticles(float   lon0, float   lat0, float   cs0, float   t0, int Nreps) {seedParticles(new float[] {lon0}, new float[] {lat0}, new float[] {cs0}, new float[] {t0}, Nreps);}
  void seedParticles(float   lon0, float   lat0, float   cs0, float[] t0, int Nreps) {seedParticles(new float[] {lon0}, new float[] {lat0}, new float[] {cs0}, t0, Nreps);}
  void seedParticles(float   lon0, float   lat0, float[] cs0, float   t0, int Nreps) {seedParticles(new float[] {lon0}, new float[] {lat0}, cs0, new float[] {t0}, Nreps);}
  void seedParticles(float   lon0, float   lat0, float[] cs0, float[] t0, int Nreps) {seedParticles(new float[] {lon0}, new float[] {lat0}, cs0, t0, Nreps);}
  void seedParticles(float   lon0, float[] lat0, float   cs0, float   t0, int Nreps) {seedParticles(new float[] {lon0}, lat0, new float[] {cs0}, new float[] {t0}, Nreps);}
  void seedParticles(float   lon0, float[] lat0, float   cs0, float[] t0, int Nreps) {seedParticles(new float[] {lon0}, lat0, new float[] {cs0}, t0, Nreps);}
  void seedParticles(float   lon0, float[] lat0, float[] cs0, float   t0, int Nreps) {seedParticles(new float[] {lon0}, lat0, cs0, new float[] {t0}, Nreps);}
  void seedParticles(float   lon0, float[] lat0, float[] cs0, float[] t0, int Nreps) {seedParticles(new float[] {lon0}, lat0, cs0, t0, Nreps);}
  void seedParticles(float[] lon0, float   lat0, float   cs0, float   t0, int Nreps) {seedParticles(lon0, new float[] {lat0}, new float[] {cs0}, new float[] {t0}, Nreps);}
  void seedParticles(float[] lon0, float   lat0, float   cs0, float[] t0, int Nreps) {seedParticles(lon0, new float[] {lat0}, new float[] {cs0}, t0, Nreps);}
  void seedParticles(float[] lon0, float   lat0, float[] cs0, float   t0, int Nreps) {seedParticles(lon0, new float[] {lat0}, cs0, new float[] {t0}, Nreps);}
  void seedParticles(float[] lon0, float   lat0, float[] cs0, float[] t0, int Nreps) {seedParticles(lon0, new float[] {lat0}, cs0, t0, Nreps);}
  void seedParticles(float[] lon0, float[] lat0, float   cs0, float   t0, int Nreps) {seedParticles(lon0, lat0, new float[] {cs0}, new float[] {t0}, Nreps);}
  void seedParticles(float[] lon0, float[] lat0, float   cs0, float[] t0, int Nreps) {seedParticles(lon0, lat0, new float[] {cs0}, t0, Nreps);}
  void seedParticles(float[] lon0, float[] lat0, float[] cs0, float   t0, int Nreps) {seedParticles(lon0, lat0, cs0, new float[] {t0}, Nreps);}

  void seedParticles(float[] lon0, float[] lat0, float[] cs0, float[] t0, int Nreps) {
    if (debug) println("seeding " + (lon0.length*lat0.length*cs0.length*t0.length*Nreps) + " particles");
    particles = new Particle[cs0.length * lat0.length * lon0.length * t0.length * Nreps];
    particlesRNKJI = new Particle[Nreps][t0.length][cs0.length][lat0.length][lon0.length];
    int m=0;
    for (int n=0; n<t0.length; n++) {
      for (int k=0; k<cs0.length; k++) {
        for (int j=0; j<lat0.length; j++) {
          for (int i=0; i<lon0.length; i++) {
            for (int r=0; r<Nreps; r++) {
              float y0j = run.lat2meters(lat0[j]);
              float x0i = run.lon2meters(lon0[i]);
              float z0 = run.cs2z(t0[n], cs0[k], y0j, x0i);
              particles[m] = new Particle(x0i, y0j, z0, t0[n], this);
              particlesRNKJI[r][n][k][j][i] = particles[m];
              particles[m].id = m;
              m++;
            }
          }
        }
      }
    }
    if (debug) println("seedParticles: done");
  }
  
  void seedParticles(String specialMode, float[] lon0, float[] lat0, float   cs0, float   t0, int Nreps) {seedParticles(specialMode, lon0, lat0, new float[] {cs0}, new float[] {t0}, Nreps);}
  void seedParticles(String specialMode, float[] lon0, float[] lat0, float   cs0, float[] t0, int Nreps) {seedParticles(specialMode, lon0, lat0, new float[] {cs0}, t0, Nreps);}
  void seedParticles(String specialMode, float[] lon0, float[] lat0, float[] cs0, float   t0, int Nreps) {seedParticles(specialMode, lon0, lat0, cs0, new float[] {t0}, Nreps);}
  void seedParticles(String specialMode, float[] lon0, float[] lat0, float[] cs0, float[] t0, int Nreps) {
    // lonLatList --------
    if (specialMode.equals("lonLatList") || specialMode.equals("latLonList")) {
      if (lon0.length==lat0.length) {
        if (debug) println("seeding " + (lon0.length*cs0.length*t0.length*Nreps) + " particles");
        particles = new Particle[cs0.length * lat0.length * t0.length * Nreps];
        int m=0;
        for (int n=0; n<t0.length; n++) {
          for (int k=0; k<cs0.length; k++) {
            for (int ji=0; ji<lat0.length; ji++) {
              for (int r=0; r<Nreps; r++) {
              float y0ji = run.lat2meters(lat0[ji]);
              float x0ji = run.lon2meters(lon0[ji]);
              float z0 = run.cs2z(t0[n], cs0[k], y0ji, x0ji);
              particles[m] = new Particle(x0ji, y0ji, z0, t0[n], this);
              particles[m].id = m;
              m++;
              }
            }
          }
        }
        if (debug) println("seedParticles: done");   
      } else {
        println("error: " + specialMode + " requires that lon0 and lat0 be the same length.");
      }   
   // other specialModes would go here
   } else {
      println("error: don't recognize" + specialMode);
    }
  }
  
  
  void seedParticles(String specialMode, String filename, float   cs0, float   t0, int Nreps) {seedParticles(specialMode, filename, new float[] {cs0}, new float[] {t0}, Nreps);}
  void seedParticles(String specialMode, String filename, float   cs0, float[] t0, int Nreps) {seedParticles(specialMode, filename, new float[] {cs0}, t0, Nreps);}
  void seedParticles(String specialMode, String filename, float[] cs0, float   t0, int Nreps) {seedParticles(specialMode, filename, cs0, new float[] {t0}, Nreps);}
  void seedParticles(String specialMode, String filename, float[] cs0, float[] t0, int Nreps) {
    if (specialMode.equals("lonLatList") || specialMode.equals("latLonList")) {
      String[] s = loadStrings(filename); // this is how to read x and y from an external file...
      float[] x = new float[s.length];
      float[] y = new float[s.length];
      for (int i=0; i<s.length; i++) {
        float[] xy = float(splitTokens(s[i]));
        x[i] = xy[0];
        y[i] = xy[1];
      }
      seedParticles(specialMode, x, y, cs0, t0, Nreps);
    } else {
      println("error: don't recognize" + specialMode);
    }
  }

  
  
  void surfaceTrap() {
    trapToSigmaLevel(0);
  }
  
  void trapToSigmaLevel(float cs) {
    for (int i=0; i<particles.length; i++) particles[i].trapToSigmaLevel(cs);
  }
  
  void trapToZLevel(float z) {
    for (int i=0; i<particles.length; i++) particles[i].trapToZLevel(z);
  }
  
  
  
  
  void createNetcdf(String ncname) {createNetcdf(ncname,1);}
  
  void createNetcdf(String ncname, int preallocSteps) {
    // expects particles[] to be initialized already
    try {
      NetcdfFileWriteable nc = NetcdfFileWriteable.createNew(ncname, false);
      if (debug) print("creating " + nc.getLocation() + "...");
      ucar.nc2.Dimension particleDim = nc.addDimension("particle", particles.length, true, false, false);
      ucar.nc2.Dimension stepDim = nc.addDimension("step", preallocSteps, true, true, false);
      for (int vi=0; vi<saveNames.length; vi++) {
        nc.addVariable(saveNames[vi], ucar.ma2.DataType.FLOAT, new ucar.nc2.Dimension[] {stepDim, particleDim});
      }
      nc.create();
      nc_close(nc);
      if (debug) println("done");
    } catch (IOException ioe) {
      println("trouble creating netcdf file " + ncname);
      println(ioe.toString());
    }
    // save initial positions
    if (debug) print("saving initial positions...");
    saveAllPositions(0); 
    if (debug) println("done");
  }

  
  
  void calcToTime(float tend) {
    NetcdfFileWriteable nc = null;
    if (autoSave) {
      try {
        nc = NetcdfFileWriteable.openExisting(ncname, false);
        if (debug) println(ncname + " found, appending");
      } catch(IOException ex) {
        if (debug) println(ncname + " not found, creating");
        createNetcdf(ncname);
        try {
          nc = NetcdfFileWriteable.openExisting(ncname, false);
        } catch(IOException ex2) {
          if (debug) println(ncname + " still not found");
        }
      }
    }
    if (multithread) {
      integrators = new Integrator[nThreads];
      for (int i=0; i<nThreads; i++) integrators[i] = new Integrator(this, nc);
      saveStep = 0;
    }
    tend = min(tend, run.lastTime());
    if (tend < run.firstTime()) return;
    float tpmin = Inf; // earliest current particle time
    for (int m=0; m<particles.length; m++) if (particles[m].active) tpmin = min(tpmin, particles[m].t());
    if (tpmin >= tend) return;
    if (tpmin < run.firstLoadedTime() || tpmin > run.lastLoadedTime()) {
      run.loadFramesAtTime(tpmin);
    }
    boolean anyActive = true;
    while (anyActive && tpmin < tend) { // until the whole integration is done...
      if (!multithread) {
        for (int m=0; m<particles.length; m++) {
          while (particles[m].t() < min(tend, run.lastLoadedTime()) && particles[m].active) { // ...integrate each particle to the end of the loaded frame (or tend, whichever comes first)
            particles[m].takeStep();
            if (autoSave && particles[m].step % saveInterval == 0) {
              savePosition(nc, particles[m].step / saveInterval, m, particles[m]); // save to netcdf file, one particle at a time
            }
          }
        }
      } else { // multithreaded version
         // send all the particles (including the out-of-range ones) to the integrators
        int curr = 0;
        for (int m=0; m<particles.length; m++) {
          integrators[curr].add(particles[m]);
          curr = (curr + 1) % nThreads;
        }
        // set them working
        if (debug) print("    starting integration...");
        ExecutorService pool = Executors.newFixedThreadPool(nThreads);
        for (int i=0; i<nThreads; i++) {
          integrators[i].endTime = min(tend, run.lastLoadedTime());
          pool.submit(integrators[i]);
        }
        if (debug) print("now...");
        // wait for them to finish
        pool.shutdown();
        float tic = millis();
        while (!pool.isTerminated()) {
          if (debug) {
            if (millis()-tic > 1000) {
              int c = 0;
              for (int i=0; i<nThreads; i++) {
                c += integrators[i].queue.size();
              }
              println("    " + c + " / " + particles.length);
              tic = millis();
            }
          }
        }
        if (debug) println("done");
        if (autoSave) {
          if (debug) print("    saving...");
          saveStep++;
          saveAllPositions(nc, saveStep);
          if (debug) println("done");
        }
      }
      
      try {
        nc.flush(); // is this necessary or helpful? not sure.
      } catch (IOException ex) {
        println("nc.flush: " + ex.toString());
      }
      
      // redetermine anyActive and tend, to decide whether to keep going
      tpmin = Inf;
      for (int m=0; m<particles.length; m++) {
        if (particles[m].active) {
          tpmin = min(tpmin, particles[m].t());   
          anyActive = true;
        }
      }   
      if (debug) println("t = " + tpmin + " (" + (tpmin/86400-run.fileTimes[0]/86400) + " d)");
      run.advance();
    }
    try {
      if (autoSave) nc.close();
    } catch (IOException ex) {
      println("nc.close: " + ex.toString());
    }
    if (debug) println("done!");
  }


  void savePosition(NetcdfFileWriteable nc, int row, int m, Particle P) {
    int[] pos = {row, m};
    try {
      for (int i=0; i<saveNames.length; i++) {
        ArrayFloat.D2 data = new ArrayFloat.D2(1,1);
        data.set(0,0,(Float)P.current.get(saveNames[i]));
        nc.write(saveNames[i], pos, data);  
      }
    } catch (IOException ioe) {
      if (debug) print("trouble saving particle " + m + " at step " + P.step + ": ");
      println(ioe.toString());
    } catch (InvalidRangeException ire) {
      if (debug) print("trouble saving particle " + m + " at step " + P.step + ": ");
      println(ire.toString());
    }
  }
  
  void savePosition(int row, int m, Particle P) {
    try {
      NetcdfFileWriteable nc = NetcdfFileWriteable.openExisting(ncname, false);
      savePosition(nc, row, m, P);
      nc.close();
    } catch (IOException ioe) {
      if (debug) print("trouble saving particle " + m + " at step " + P.step + ": ");
      println(ioe.toString());
    }
  }
  
  void saveAllPositions(int row) {
    int[] pos = {row,0};
    try {
      NetcdfFileWriteable nc = NetcdfFileWriteable.openExisting(ncname, false);
      saveAllPositions(nc, row);
      nc.close();
    } catch (IOException ioe) {
      if (debug) print("trouble saving particles at row " + row + ": ");
      println(ioe.toString());
    }
  }
  
  void saveAllPositions(NetcdfFileWriteable nc, int row) {
    int[] pos = {row,0};
    try {
      for (int i=0; i<saveNames.length; i++) {
        ArrayFloat.D2 data = new ArrayFloat.D2(1,particles.length);
        for (int j=0; j<particles.length; j++) data.set(0,j,(Float)particles[j].current.get(saveNames[i]));
        nc.write(saveNames[i], pos, data);
      }
    } catch (IOException ioe) {
      if (debug) print("trouble saving particles at row " + row + ": ");
      println(ioe.toString());  
    } catch (InvalidRangeException ire) {
      if (debug) print("trouble saving particles at row " + row + ": ");
      println(ire.toString());        
    }
  }

  
}

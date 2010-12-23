class ParticleExpt {

  ROMSRun run;
  Particle[] particles;
  String ncname = "myexpt.nc";
  int saveInterval = 1;
  int preallocSteps = 1;
  boolean saveFirstLastOnly = false;
  boolean autoSave = true; // the normal way of saving
  String[] saveNames = {"lon","lat","z","cs","t","H","mask"};
  // the following variables are copied from the experiment in the Particle constructor
  float dt = 400;
  String[] tracerNames = new String[0];
  String[] tracerInterpMode = new String[0];
  float[] tracerInterpDepth = new float[0];
  
  ParticleExpt() {}
  
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
    run.tracerNames = tracerNames;
    run.reload();
  }
  
  String tracerSaveName(int i) {
    String name = tracerNames[i];
    if (tracerInterpMode[i] !="3D") name += "_" + tracerInterpMode[i] + "_" + tracerInterpDepth[i];
    name = strrep(name,'.',"");
    return name;
  }

  void seedParticles(float[] lon0, float[] lat0, float[] cs0, float t0, int Nreps) {
    particles = new Particle[cs0.length * lat0.length * lon0.length * Nreps];
    int m=0;
    for (int k=0; k<cs0.length; k++) {
      for (int j=0; j<lat0.length; j++) {
        for (int i=0; i<lon0.length; i++) {
          for (int r=0; r<Nreps; r++) {
            float y0j = run.lat2meters(lat0[j]);
            float x0i = run.lon2meters(lon0[i]);
            float z0 = run.cs2z(t0, cs0[k], y0j, x0i);
            particles[m] = new Particle(x0i, y0j, z0, t0, this);
            m++;
          }
        }
      }
    }
    if (autoSave) createNetcdf(ncname);
  }
  
  void seedParticles(float[] lon0, float[] lat0, float[] cs0, float[] t0, int Nreps) {
    particles = new Particle[cs0.length * lat0.length * lon0.length * t0.length * Nreps];
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
              m++;
            }
          }
        }
      }
    }
    if (autoSave) createNetcdf(ncname);
  }
  
  void createNetcdf(String ncname) {createNetcdf(ncname,1);}
  
  void createNetcdf(String ncname, int preallocSteps) {
    if (debug) print("creating output file " + ncname + "...");
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
    try {
      NetcdfFileWriteable nc = null;
      if (autoSave) nc = NetcdfFileWriteable.openExisting(ncname, false);    
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
        for (int m=0; m<particles.length; m++) {
          Particle P = particles[m];
          while (P.t() < min(tend, run.lastLoadedTime()) && P.active) { // ...integrate each particle to the end of the loaded frame (or tend, whichever comes first)
            P.takeStep();
            if (autoSave && P.step % saveInterval == 0) {
              savePosition(nc, P.step / saveInterval, m, P); // save to netcdf file, one particle at a time
            }
          }
        }
        // redetermine anyActive and tend, to decide whether to keep going
        tpmin = Inf;
        for (int m=0; m<particles.length; m++) {
          if (particles[m].active) {
            tpmin = min(tpmin, particles[m].t());   
            anyActive = true;
          }
        }   
        println("tpmin = " + tpmin + " (" + (tpmin/86400-run.fileTimes[0]/86400) + " d)");
        run.advance();
      }
      if (autoSave) nc.close();
      if (debug) println("done");
    } catch (IOException ioe) {
      if (debug) print("trouble saving particles: " + ioe.toString());
    }
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
      for (int i=0; i<saveNames.length; i++) {
        ArrayFloat.D2 data = new ArrayFloat.D2(1,particles.length);
        for (int j=0; j<particles.length; j++) data.set(0,j,(Float)particles[j].current.get(saveNames[i]));
        nc.write(saveNames[i], pos, data);
      }
      nc.close();
    } catch (IOException ioe) {
      if (debug) print("trouble saving particles at row " + row + ": ");
      println(ioe.toString());  
    } catch (InvalidRangeException ire) {
      if (debug) print("trouble saving particles at row " + row + ": ");
      println(ire.toString());        
    }
  }
  
}

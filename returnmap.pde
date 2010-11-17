class ReturnMapMaker extends ParticleExpt {
  
  int xySubsample = 1;
  String ncbasename;
  float[] X, Y;
  boolean surfaceTrapped = false;
    
  ReturnMapMaker(String dirname, int ncn0, int ncn1, String ncbasename, int xySubsample) {
    linkToRun(dirname + "ocean_his_", ncn0, ncn1);
    this.ncbasename = ncbasename;
    this.xySubsample = xySubsample;
    X = new float[run.I / xySubsample];
    for (int i=0; i<X.length; i++) X[i] = run.lon[i*xySubsample];
    Y = new float[run.J / xySubsample];
    for (int j=0; j<Y.length; j++) Y[j] = run.lat[j*xySubsample];
    autoSave = false; // the file output gets special handling
  }  
  
  
  String[] saveNames() {
    return new String[] {"lon","lat","t"};
  }
  
  float[] saveValues(Particle P) {
    return new float[] {P.lon, P.lat, P.t};
  }
  
  void makeFirstLast(float tstart, float dt, float[] cslevels, int Nreps, float internalTimestep, String filesuffix) {
    // calc particle trajectories, saving first and last positions to a normal output file
    ncname = ncbasename + "." + filesuffix + "firstlast.nc";
    seedParticles(X, Y, cslevels, tstart, Nreps);
    createNetcdf(ncname); // create the file and save initial positions. If autoSave were true, we wouldn't have to do this manually
    for (int i=0; i<particles.length; i++) particles[i].dt = internalTimestep;
    calcToTime(tstart+dt);   
    saveAllPositions(1); // save final particle positions 
  }
  
 
  void makeMap(float tstart, float dt, float[] cslevels, int Nreps, float internalTimestep, String filesuffix) {
    makeMap(tstart, dt, cslevels, Nreps, internalTimestep, filesuffix, true);
  }

  void makeMap(float tstart, float dt, float[] cslevels, int Nreps, float internalTimestep, String filesuffix, boolean includeIndices) {
    // calc particle trajectories, but save output in the form used by flowWeaver 0.4
    
    // calc final positions, but don't save them anywhere yet
    seedParticles(X, Y, cslevels, tstart, Nreps);
    for (int i=0; i<particles.length; i++) particles[i].dt = internalTimestep;
    if (surfaceTrapped) {
      for (int i=0; i<particles.length; i++) {
        particles[i].surfaceTrapped = true;
        particles[i].diffusive = false;
      }
    }
    calcToTime(tstart+dt);

    // assemble xnext, ynext for the wet cells
    int J = Y.length, I = X.length, K = cslevels.length;
    float[][] xn = new float[K*Nreps][J*I]; // treat multiple cs levels as another kind of replicate. J*I is the max size this array might be, if there are no land cells
    float[][] yn = new float[K*Nreps][J*I];
    int[] ind1 = new int[J*I];
    float small = 1e-4; // about 10 meters
    int m=-1, mgood=-1;
    for (int k=0; k<cslevels.length; k++) { // the order of looping has to match the way the particles are created in seedParticles()
      for (int j=0; j<Y.length; j++) {
        for (int i=0; i<X.length; i++) {
          for (int n=0; n<Nreps; n++) {
            m++;
            if (particles[m].mask0 > small) {
              mgood++;
              xn[k*Nreps+n][mgood] = particles[m].lon;
              yn[k*Nreps+n][mgood] = particles[m].lat;
              ind1[mgood] = i * J + j;
            }
          }
        }
      }
    }
    int Ncells = mgood+1;

    // create the file
    ncname = ncbasename + "." + filesuffix + ".map.nc";
    try {
      NetcdfFileWriteable nc = NetcdfFileWriteable.createNew(ncname, false);
      if (debug) print("writing " + nc.getLocation() + "...");
      ucar.nc2.Dimension JDim = nc.addDimension("J", J, true, false, false);
      ucar.nc2.Dimension IDim = nc.addDimension("I", I, true, false, false);
      ucar.nc2.Dimension NrepsDim = nc.addDimension("Nreps", Nreps, true, false, false);
      ucar.nc2.Dimension oneDim = nc.addDimension("one", 1, true, false, false);
      ucar.nc2.Dimension cellDim = nc.addDimension("cell", Ncells, true, false, false);      
      nc.addVariable("x", ucar.ma2.DataType.FLOAT, new ucar.nc2.Dimension[] {IDim});
      nc.addVariable("y", ucar.ma2.DataType.FLOAT, new ucar.nc2.Dimension[] {JDim});
      if (includeIndices) nc.addVariable("ind", ucar.ma2.DataType.INT, new ucar.nc2.Dimension[] {cellDim});
      nc.addVariable("xnext", ucar.ma2.DataType.SHORT, new ucar.nc2.Dimension[] {NrepsDim, cellDim});
      nc.addVariable("ynext", ucar.ma2.DataType.SHORT, new ucar.nc2.Dimension[] {NrepsDim, cellDim});
      nc.addVariable("timestep", ucar.ma2.DataType.FLOAT, new ucar.nc2.Dimension[] {oneDim});
      nc.addVariable("numAlternates", ucar.ma2.DataType.FLOAT, new ucar.nc2.Dimension[] {oneDim});
      nc.addVariable("numCells", ucar.ma2.DataType.FLOAT, new ucar.nc2.Dimension[] {oneDim});
      nc.create();
      
      // write coordinates and such
      ArrayFloat.D1 data = new ArrayFloat.D1(I);
      for (int i=0; i<I; i++) data.set(i,X[i]);
      nc.write("x", new int[] {0}, data);
      data = new ArrayFloat.D1(J);
      for (int j=0; j<J; j++) data.set(j,Y[j]);
      nc.write("y", new int[] {0}, data);
      data = new ArrayFloat.D1(1);
      data.set(0,dt);
      nc.write("timestep", new int[] {0}, data);
      data.set(0,K*Nreps);
      nc.write("numAlternates", new int[] {0}, data);
      data.set(0,Ncells);
      nc.write("numCells", new int[] {0}, data); 
      
      // rescale xn,yn into short ints (cuts the file size in half) and save them
      ArrayShort.D2 xnext = new ArrayShort.D2(K*Nreps, Ncells);
      ArrayShort.D2 ynext = new ArrayShort.D2(K*Nreps, Ncells);
      ArrayInt.D1 ind = new ArrayInt.D1(Ncells);
      for (m=0; m<Ncells; m++) {
        if (includeIndices) ind.set(m, ind1[m]);
        for (int r=0; r<K*Nreps; r++) {
          xnext.set(r,m,(short) round(map(constrain(xn[r][m], X[0], X[I-1]), X[0], X[I-1], -32768, 32767)));
          ynext.set(r,m,(short) round(map(constrain(yn[r][m], Y[0], Y[J-1]), Y[0], Y[J-1], -32768, 32767)));
        }
      }
      if (includeIndices) nc.write("ind", new int[] {0}, ind);   
      nc.write("xnext", new int[] {0,0}, xnext);      
      nc.write("ynext", new int[] {0,0}, ynext);           
      nc_close(nc);
      if (debug) println("done");
    } catch (IOException ioe) {
      println("trouble writing netcdf file " + ncname);
      println(ioe.toString());
    } catch (InvalidRangeException ire) {
      println("trouble writing netcdf file " + ncname);
      println(ire.toString());
    }

  }
  
}

class ReturnMapMaker extends ParticleExpt {
  
  int xySubsample = 1;
  String ncbasename;
  float[] X, Y;
  
    
  ReturnMapMaker(String dirname, int ncn0, int ncn1, String ncbasename, int xySubsample) {
    linkToRun(dirname + "ocean_his_", ncn0, ncn1);
    this.ncbasename = ncbasename;
    this.xySubsample = xySubsample;
    X = new float[run.I / xySubsample];
    for (int i=0; i<X.length; i++) X[i] = run.lon[i*xySubsample];
    Y = new float[run.J / xySubsample];
    for (int j=0; j<Y.length; j++) Y[j] = run.lat[j*xySubsample];
    autoSave = false; // the file output gets special handling
    saveNames = new String[] {"lon","lat","t"};
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
    // calc particle trajectories, but save output in the form used by flowWeaver 0.5
    // (this removes the optimizations used in flowWeaver 0.4: i.e., entire grid saved, not just wet points, and values are not rescaled to short ints)
    
    // calc final positions, but don't save them anywhere yet
    if (debug) println("  calculating particles starting at " + tstart + " (" + (tstart/86400.) + " d) for " + dt);
    int nAlt = cslevels.length * Nreps;
    dt = internalTimestep;
    seedParticles(X, Y, cslevels, tstart, Nreps);
    for (int i=0; i<particles.length; i++) {
      particles[i].trapToSigmaLevel(particles[i].cs());
    }
    calcToTime(tstart+dt);


    // create the file
    ncname = ncbasename + "." + filesuffix + ".map.nc";
    int I = X.length;
    int J = Y.length;
    try {
      NetcdfFileWriteable nc = NetcdfFileWriteable.createNew(ncname, false);
      if (debug) print("writing " + nc.getLocation() + "...");
      ucar.nc2.Dimension JDim = nc.addDimension("J", J, true, false, false);
      ucar.nc2.Dimension IDim = nc.addDimension("I", I, true, false, false);
      ucar.nc2.Dimension nAltDim = nc.addDimension("nAlt", nAlt, true, false, false); // replicates and multiple depths
      ucar.nc2.Dimension oneDim = nc.addDimension("one", 1, true, false, false);
      nc.addVariable("x", ucar.ma2.DataType.FLOAT, new ucar.nc2.Dimension[] {IDim});
      nc.addVariable("y", ucar.ma2.DataType.FLOAT, new ucar.nc2.Dimension[] {JDim});
      nc.addVariable("xnext", ucar.ma2.DataType.FLOAT, new ucar.nc2.Dimension[] {nAltDim, JDim, IDim});
      nc.addVariable("ynext", ucar.ma2.DataType.FLOAT, new ucar.nc2.Dimension[] {nAltDim, JDim, IDim});
      nc.addVariable("timestep", ucar.ma2.DataType.FLOAT, new ucar.nc2.Dimension[] {oneDim});
      nc.addVariable("numAlternates", ucar.ma2.DataType.FLOAT, new ucar.nc2.Dimension[] {oneDim});
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
      data.set(0,nAlt);
      nc.write("numAlternates", new int[] {0}, data);
      
     // assemble and write xnext, ynext
      if (debug) print("xnext, ynext...");
      ArrayFloat.D3 xnext = new ArrayFloat.D3(nAlt, J, I);
      ArrayFloat.D3 ynext = new ArrayFloat.D3(nAlt, J, I);
      int m=-1;
      for (int r=0; r<Nreps; r++) {
        for (int k=0; k<cslevels.length; k++) {
          for (int j=0; j<Y.length; j++) {
            for (int i=0; i<X.length; i++) {
              m++;
              xnext.set(k*Nreps+r,j,i,particlesRNKJI[r][0][k][j][i].lon());
              ynext.set(k*Nreps+r,j,i,particlesRNKJI[r][0][k][j][i].lat());
            }
          }
        }
      }
      nc.write("xnext", new int[] {0,0,0}, xnext);      
      nc.write("ynext", new int[] {0,0,0}, ynext); 
      
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

class GETM_2DRun extends ROMSRun {
  
  String filename;
  float[][] Ubar_n0, Ubar_n1, Vbar_n0, Vbar_n1;
  // don't use filenames[] or [U,V]_n[0,1]
  float[][][] Ubar_full, Vbar_full, zeta_full, mask_full;
  
  GETM_2DRun(String runname) {
    print("loading " + runname + "...");
    filename = runname;
    NetcdfFile nc = nc_open(filename);
    fileTimes = nc_read1D(nc,"time");
    H = nc_read2D(nc,"bathymetry");
    X = firstRow(nc_read2D(nc,"xc"));
    Y = firstCol(nc_read2D(nc,"yc"));
    I = X.length;
    J = Y.length;
    Xu = X;
    Yu = Y;
    Iu = I;
    Ju = J;
    Xv = X;
    Yv = Y;
    Iv = I;
    Jv = J;
    meanLon = 0;
    meanLat = 0;
    cosMeanLat = 1;
    lon = meters2lon(X);
    lat = meters2lat(Y);
    lonu = lon;
    latu = lat;
    lonv = lon;
    latv = lat;
    Cs = new float[] {-1, 0}; // these shouldn't be used, but someone (like particle.interpEverything) might try to access them
    Csw = new float[] {-1, 0};
    
    // preload all the data. The business about loading one pair of frames at a time is thus unnecessary,
    // but the structure is retained in order to make it easy to extend this to massive 3D getm runs too big to preload.
    // now, loadFrame just updates t and ncn.
    print("...loading u,v,zeta...");
    Ubar_full = nc_read3D(nc,"u");
    Vbar_full = nc_read3D(nc,"v");
    zeta_full = nc_read3D(nc,"elev");

    // construct the mask and eliminate bad velocity points
    mask_full = arrayFill(fileTimes.length, J, I, 1);
    for (int n=0; n<fileTimes.length; n++) {
      for (int j=0; j<J; j++) {
        for (int i=0; i<I; i++) {
          if (Ubar_full[n][j][i] == -9999 || Vbar_full[n][j][i] == -9999 || zeta_full[n][j][i] == -9999) {
            mask_full[n][j][i] = 0;
            Ubar_full[n][j][i] = 0;
            Vbar_full[n][j][i] = 0;
            zeta_full[n][j][i] = 0;
          }
        }
      }
    }
    
    nc_close(nc);
    loadFrame(0);
    advance();
    println("done");
  }
  
  void loadFrame(int ncn) {
    ncn_n1 = ncn;
    t_n1 = fileTimes[ncn];
  }
  
  void advance() {
    if (ncn_n1 < fileTimes.length-1) {
      ncn_n0 = ncn_n1;
      t_n0 = t_n1;
      loadFrame(ncn_n1 + 1);
    } else {
      if (debug) println("last frames loaded: can't advance");
    }
  }

  
  
  float cs2z(float t, float cs, float y, float x) {return 0;}
  float z2cs(float t, float z, float y, float x) {return 0;}

  float W(float t, int k, int j, int i) {return 0;}
  float Ks(float t, int k, int j, int i) {return 0;}  

  float interpW(float t, float cs, float y, float x) {return 0;}
  float interpKs(float t, float cs, float y, float x) {return 0;}

  float zeta(float t, int j, int i) {
    float f = constrain((t - t_n0) / (t_n1 - t_n0), 0, 1);
    return (1-f) * zeta_full[ncn_n0][j][i] + f * zeta_full[ncn_n1][j][i];        
  }
  
  float mask(float t, int j, int i) {
    float f = constrain((t - t_n0) / (t_n1 - t_n0), 0, 1);
    return (1-f) * mask_full[ncn_n0][j][i] + f * mask_full[ncn_n1][j][i];        
  }

  float Ubar(float t, int j, int i) {
    float f = constrain((t - t_n0) / (t_n1 - t_n0), 0, 1);
    return (1-f) *  Ubar_full[ncn_n0][j][i] + f * Ubar_full[ncn_n1][j][i];
  }
  
  float Vbar(float t, int j, int i) {
    float f = constrain((t - t_n0) / (t_n1 - t_n0), 0, 1);
    return (1-f) * Vbar_full[ncn_n0][j][i] + f * Vbar_full[ncn_n1][j][i];
  }
    
  float interpU(float t, float cs, float y, float x) {
    int ii = findIndexBefore(Xu, x);
    int jj = findIndexBefore(Yv, y);
    if (ii<0 || jj<0) return 0;
    float f = constrain((x - Xu[ii]) / (Xu[ii+1] - Xu[ii]), 0, 1);
    float ui = f * Ubar(t, jj+1, ii+1) + (1-f) * Ubar(t, jj+1, ii);
    return ui;
  }
  
  float interpV(float t, float cs, float y, float x) {
    int ii = findIndexBefore(Xu, x);
    int jj = findIndexBefore(Yv, y);
    if (ii<0 || jj<0) return 0;
    float f = constrain((y - Yv[jj]) / (Yv[jj+1] - Yv[jj]), 0, 1);
    float vi = f * Vbar(t, jj+1, ii+1) + (1-f) * Vbar(t, jj, ii+1);
    return vi;
  }
  


}

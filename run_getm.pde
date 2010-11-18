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
    // but it's retained in order to make it easy to extend this to massive 3D getm runs too big to preload
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
    Ubar_n1 = new float[Ju][Iu];
    for (int j=0; j<Ju; j++) {
      for (int i=0; i<Iu; i++) {
        Ubar_n1[j][i] = Ubar_full[ncn][j][i];
      }
    }
    Vbar_n1 = new float[Jv][Iv];
    for (int j=0; j<Jv; j++) {
      for (int i=0; i<Iv; i++) {
        Vbar_n1[j][i] = Vbar_full[ncn][j][i];
      }
    }
    zeta_n1 = new float[J][I];
    mask_n1 = new float[J][I];
    for (int j=0; j<J; j++) {
      for (int i=0; i<I; i++) {
        zeta_n1[j][i] = zeta_full[ncn][j][i];
        mask_n1[j][i] = mask_full[ncn][j][i];
      }
    }
  }
  
  void advance() {
    if (ncn_n1 < fileTimes.length-1) {
      Ubar_n0 = Ubar_n1;
      Vbar_n0 = Vbar_n1;
      zeta_n0 = zeta_n1;
      ncn_n0 = ncn_n1;
      t_n0 = t_n1;
      mask_n0 = mask_n1;
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

  float U(float t, int k, int j, int i) {
    float f = constrain((t - t_n0) / (t_n1 - t_n0), 0, 1);
    return (1-f) *  Ubar_n0[j][i] + f * Ubar_n1[j][i];
  }
  
  float V(float t, int k, int j, int i) {
    float f = constrain((t - t_n0) / (t_n1 - t_n0), 0, 1);
    return (1-f) * Vbar_n0[j][i] + f * Vbar_n1[j][i];
  }
    
  float interpU(float t, float cs, float y, float x) {
    int ii = findIndexBefore(Xu, x);
    int jj = findIndexBefore(Yv, y);
    if (ii<0 || jj<0) return 0;
    float f = constrain((x - Xu[ii]) / (Xu[ii+1] - Xu[ii]), 0, 1);
    float ui = f * U(t, 0, jj+1, ii+1) + (1-f) * U(t, 0, jj+1, ii); // dummy argument for cs
    return ui;
  }
  
  float interpV(float t, float cs, float y, float x) {
    int ii = findIndexBefore(Xu, x);
    int jj = findIndexBefore(Yv, y);
    if (ii<0 || jj<0) return 0;
    float f = constrain((y - Yv[jj]) / (Yv[jj+1] - Yv[jj]), 0, 1);
    float vi = f * V(t, 0, jj+1, ii+1) + (1-f) * V(t, 0, jj, ii+1); // dummy argument for cs
    return vi;
  }
  


}

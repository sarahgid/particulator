class ROMSRun {

  String[] filenames;
  float[] fileTimes;
  
  // grid
  int I,J,K, Iu,Ju, Iv,Jv, Kw;
  float[] lon,lat,X,Y,  lonu,latu,Xu,Yu,  lonv,latv,Xv,Yv,  Cs,Csw;
  float[][] H;
  float meanLon, meanLat, cosMeanLat;

  // flow data from two frames
  HashMap frame0, frame1;
  float t_n0, t_n1;
  int ncn_n0, ncn_n1;
  
  // additional tracers to load
  String[] tracerNames = new String[0];

  ROMSRun() {}

  ROMSRun(String basename, int ncn0, int ncn1) {
    // locate files and load timebase
    print("reading timebase");
    filenames = new String[ncn1-ncn0+1];
    fileTimes = new float[filenames.length];
    for (int i=0; i<filenames.length; i++) {
      filenames[i] = basename + nf(ncn0+i,4) + ".nc";
      NetcdfFile nc = nc_open(filenames[i]);
      fileTimes[i] = nc_readOne(nc, "ocean_time");
      nc_close(nc);
      if (i % 100 == 0) print(".");
    }
    println();
    println("" + filenames.length + " files found from " + filenames[0] + " to " + filenames[filenames.length-1]);
    println("times from " + fileTimes[0] + " to " + fileTimes[fileTimes.length-1]);
    
    // read the grid from the first file
    NetcdfFile nc = nc_open(filenames[0]);

    lon = firstRow(nc_read2D(nc,"lon_rho"));
    lat = firstCol(nc_read2D(nc,"lat_rho"));
    meanLon = 0; for (int i=0; i<lon.length; i++) meanLat += lon[i]; meanLon /= lon.length;
    meanLat = 0; for (int j=0; j<lat.length; j++) meanLat += lat[j]; meanLat /= lat.length;
    cosMeanLat = cos(meanLat / 180 * PI);
    X = lon2meters(lon);
    Y = lat2meters(lat);
    I = X.length;
    J = Y.length;

    lonu = firstRow(nc_read2D(nc,"lon_u"));  
    latu = firstCol(nc_read2D(nc,"lat_u"));  
    Xu = lon2meters(lonu);
    Yu = lat2meters(latu);
    Iu = Xu.length;
    Ju = Yu.length;

    lonv = firstRow(nc_read2D(nc,"lon_v"));  
    latv = firstCol(nc_read2D(nc,"lat_v"));  
    Xv = lon2meters(lonv);
    Yv = lat2meters(latv);
    Iv = Xv.length;
    Jv = Yv.length;

    Cs = nc_read1D(nc,"Cs_r");                           
    Csw = nc_read1D(nc,"Cs_w");                       
    K = Cs.length;
    Kw = Csw.length;
    H = nc_read2D(nc,"h");

    nc_close(nc);
        
    // load first two frames of velocity data
    frame1 = new HashMap();
    loadFrame(0);
    advance();
    
    if (debug) println("first two frames loaded");
  }
  
  
  
  float firstLoadedTime() {return t_n0;}
  float lastLoadedTime() {return t_n1;}
  float firstTime() {return fileTimes[0];}
  float lastTime() {return fileTimes[fileTimes.length-1];}
  
  void loadFrame(int ncn) {
    if (debug) print("loading frame " + ncn + "...");
    NetcdfFile nc = nc_open(filenames[ncn]);
    ncn_n1 = ncn;
    t_n1 = fileTimes[ncn];
    float[][][] data = nc_read3D(nc, "u");
      if (data==null) data = zeros(K,Ju,Iu);
      frame1.put("U",finitize(data));
    data = nc_read3D(nc, "v");
      if (data==null) data = zeros(K,Jv,Iv);
      frame1.put("V",finitize(data));
    data = nc_read3D(nc, "w");
      if (data==null) data = zeros(Kw,J,I);
      frame1.put("W",finitize(data));
    data = nc_read3D(nc, "AKs");
      if (data==null) data = zeros(Kw,J,I);
      frame1.put("Ks",finitize(data));
    float[][] data2 = nc_read2D(nc, "zeta");
      if (data2==null) data2 = zeros(J,I);
      frame1.put("zeta",finitize(data2));
    data2 = nc_read2D(nc, "mask_rho");
      if (data2==null) data2 = arrayFill(J,I,1);
      frame1.put("mask",finitize(data2));
    for (int i=0; i<tracerNames.length; i++) {
      data = nc_read3D(nc, tracerNames[i]);
      if (data==null) data = zeros(K,J,I);
      frame1.put(tracerNames[i],finitize(data));
    }
    nc_close(nc);
    if (debug) println("done");
  }
  
  void loadFramesAtTime(float ti) {
    int ncni = findIndexBefore(fileTimes, ti);
    if (ncni == -1) {
     if (debug) println("time " + ti + " is out of range");
     return;
    }
    if (ncni == ncn_n1) {
      advance();
    } else {
      loadFrame(ncni);
      advance();
    }
  }
     
  void advance() {
    if (ncn_n1 < fileTimes.length-1) {
      frame0 = (HashMap)frame1.clone();
      ncn_n0 = ncn_n1;
      t_n0 = t_n1;
      loadFrame(ncn_n1 + 1);
    } else {
      if (debug) println("last frames loaded");
    }
  }
  
  void reload() {
    loadFrame(ncn_n0);
    advance();
  }
    
  
  float H(int j, int i) {return H[j][i];}
  float zeta(float t, int j, int i) {return get_tji("zeta",t,j,i);}
  float mask(float t, int j, int i) {return get_tji("mask",t,j,i);}       

  float U(float t, int k, int j, int i) {return get_tkji("U",t,k,j,i);}
  float V(float t, int k, int j, int i) {return get_tkji("V",t,k,j,i);}
  float W(float t, int k, int j, int i) {return get_tkji("W",t,k,j,i);}
  float Ks(float t, int k, int j, int i) {return get_tkji("Ks",t,k,j,i);}
  
  float get_tji(String name, float t, int j, int i) {
    float[][] data_n0 = (float[][])frame0.get(name);
    float[][] data_n1 = (float[][])frame1.get(name);
    float f = constrain((t - t_n0) / (t_n1 - t_n0), 0, 1);
    return (1-f) * data_n0[j][i] + f * data_n1[j][i];            
  }
  
  float get_tkji(String name, float t, int k, int j, int i) {
    float[][][] data_n0 = (float[][][])frame0.get(name);
    float[][][] data_n1 = (float[][][])frame1.get(name);
    float f = constrain((t - t_n0) / (t_n1 - t_n0), 0, 1);
    return (1-f) * data_n0[k][j][i] + f * data_n1[k][j][i];            
  }
  
  
  
  float interpH(float y, float x) {
    int ii = findIndexBefore(X, x);
    int jj = findIndexBefore(Y, y);
    if (ii<0 || jj<0) return 0;
    float fx = constrain((x - X[ii]) / (X[ii+1] - X[ii]), 0, 1);
    float fy = constrain((y - Y[jj]) / (Y[jj+1] - Y[jj]), 0, 1);
    float Hi = fx * fy * H(jj+1,ii+1)  +  (1-fx) * fy * H(jj+1,ii)  +  fx * (1-fy) * H(jj,ii+1)  +  (1-fx) * (1-fy) * H(jj,ii);
    return Hi;
  }
  
  float interpZeta(float t, float y, float x) {
    int ii = findIndexBefore(X, x);
    int jj = findIndexBefore(Y, y);
    if (ii<0 || jj<0) return 0;
    float fx = constrain((x - X[ii]) / (X[ii+1] - X[ii]), 0, 1);
    float fy = constrain((y - Y[jj]) / (Y[jj+1] - Y[jj]), 0, 1);
    float zetai = fx * fy * zeta(t,jj+1,ii+1)  +  (1-fx) * fy * zeta(t,jj+1,ii)  +  fx * (1-fy) * zeta(t,jj,ii+1)  +  (1-fx) * (1-fy) * zeta(t,jj,ii);
    return zetai;    
  }
  
  float interpMask(float t, float y, float x) {
    int ii = findIndexBefore(X, x);
    int jj = findIndexBefore(Y, y);
    if (ii<0 || jj<0) return 0;
    float fx = constrain((x - X[ii]) / (X[ii+1] - X[ii]), 0, 1);
    float fy = constrain((y - Y[jj]) / (Y[jj+1] - Y[jj]), 0, 1);
    float maski = fx * fy * mask(t,jj+1,ii+1)  +  (1-fx) * fy * mask(t,jj+1,ii)  +  fx * (1-fy) * mask(t,jj,ii+1)  +  (1-fx) * (1-fy) * mask(t,jj,ii);
    return maski;    
  }

  float interpU(float t, float cs, float y, float x) {
    int ii = findIndexBefore(Xu, x);
    int jj = findIndexBefore(Yv, y);
    int kk = findIndexBefore(Csw, cs);
    if (ii<0 || jj<0 || kk<0) return 0;
    float f = constrain((x - Xu[ii]) / (Xu[ii+1] - Xu[ii]), 0, 1);
    float ui = f * U(t, kk, jj+1, ii+1) + (1-f) * U(t, kk, jj+1, ii); // Y[jj+1] is between Yv[jj] and Yv[jj+1], but Cs[kk] is between Csw[kk] and Csw[kk+1]
    return ui;
  }
  
  float interpV(float t, float cs, float y, float x) {
    int ii = findIndexBefore(Xu, x);
    int jj = findIndexBefore(Yv, y);
    int kk = findIndexBefore(Csw, cs);
    if (ii<0 || jj<0 || kk<0) return 0;
    float f = constrain((y - Yv[jj]) / (Yv[jj+1] - Yv[jj]), 0, 1);
    float vi = f * V(t, kk, jj+1, ii+1) + (1-f) * V(t, kk, jj, ii+1);  
    return vi;
  }
  
  float interpW(float t, float cs, float y, float x) {
    int ii = findIndexBefore(Xu, x);
    int jj = findIndexBefore(Yv, y);
    int kk = findIndexBefore(Csw, cs);
    if (ii<0 || jj<0 || kk<0) return 0;
    float f = constrain((cs - Csw[kk]) / (Csw[kk+1] - Csw[kk]), 0, 1);
    float wi = f * W(t, kk+1, jj+1, ii+1) + (1-f) * W(t, kk, jj+1, ii+1);   
    return wi;
  }
  
  float interpKs(float t, float cs, float y, float x) { // diffusivity
    int ii = findIndexBefore(Xu, x);
    int jj = findIndexBefore(Yv, y);
    int kk = findIndexBefore(Csw, cs);
    if (ii<0 || jj<0 || kk<0) return 0;
    float f = constrain((cs - Csw[kk]) / (Csw[kk+1] - Csw[kk]), 0, 1);
    float ksi = f * Ks(t, kk+1, jj+1, ii+1) + (1-f) * Ks(t, kk, jj+1, ii+1);   
    return ksi;
  }
  
  float interpTracer(String name, float t, float cs, float y, float x) { // some other tracer on the rho grid
    int ii = findIndexBefore(X, x);
    int jj = findIndexBefore(Y, y);
    int kk = findIndexBefore(Cs, constrain(cs, Cs[0], Cs[K-1]));
    if (ii<0 || jj<0 || kk<0) return 0;
    float fx = constrain((x - X[ii]) / (X[ii+1] - X[ii]), 0, 1);
    float fy = constrain((y - Y[jj]) / (Y[jj+1] - Y[jj]), 0, 1);
    float c0 = fx * fy * get_tkji(name,t,kk,jj+1,ii+1)    +  (1-fx) * fy * get_tkji(name,t,kk,jj+1,ii)    +  fx * (1-fy) * get_tkji(name,t,kk,jj,ii+1)    +  (1-fx) * (1-fy) * get_tkji(name,t,kk,jj,ii);
    float c1 = fx * fy * get_tkji(name,t,kk+1,jj+1,ii+1)  +  (1-fx) * fy * get_tkji(name,t,kk+1,jj+1,ii)  +  fx * (1-fy) * get_tkji(name,t,kk+1,jj,ii+1)  +  (1-fx) * (1-fy) * get_tkji(name,t,kk+1,jj,ii);
    float f = constrain((cs - Cs[kk]) / (Cs[kk+1] - Cs[kk]), 0, 1);
    float c = f * c1 + (1-f) * c0;   
    return c;
  }

  
  
  float cs2z(float t, float cs, float y, float x) {
    float zetai = interpZeta(t, y, x);
    float Hi = interpH(y, x);
    float zi = (Hi + zetai) * cs + zetai;
    return zi;
  }
  
  float z2cs(float t, float z, float y, float x) {
    float zetai = interpZeta(t, y, x);
    float Hi = interpH(y, x);
    if (!isfinite(Hi+zetai)) return 0;
    float csi = (z - zetai) / (Hi + zetai);
    return csi;
  }
  
  
  
  // in a large grid, these might be off by a lot as absolute measurements.
  // but they aren't used in particle integration, only in I/O.
  
  float lon2meters(float lon) {
    return (lon - meanLon) * 111325 * cosMeanLat;
  }
  float[] lon2meters(float[] lon) {
    float[] R = new float[lon.length];
    for (int i=0; i<R.length; i++) R[i] = lon2meters(lon[i]);
    return R;
  }
  
  float lat2meters(float lat) {
    return (lat - meanLat) * 111325;
  }
  float[] lat2meters(float[] lat) {
    float[] R = new float[lat.length];
    for (int i=0; i<R.length; i++) R[i] = lat2meters(lat[i]);
    return R;
  }
  
  float meters2lon(float x) {
    return x / 111325 / cosMeanLat + meanLon;
  }
  float[] meters2lon(float[] x) {
    float[] R = new float[x.length];
    for (int i=0; i<R.length; i++) R[i] = meters2lon(x[i]);
    return R;
  }
  
  float meters2lat(float y) {
    return y / 111325 + meanLat;
  }
  float[] meters2lat(float[] y) {
    float[] R = new float[y.length];
    for (int i=0; i<R.length; i++) R[i] = meters2lat(y[i]);
    return R;
  }
  
}

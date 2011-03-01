void pugetSoundTest() {
  ParticleExpt expt = new ParticleExpt();
  expt.linkToRun("/Users/neil/pnwtox/psmid_jul8_18/ocean_his_", 1633, 1873);
  expt.ncname = "/Users/neil/Desktop/particles_ps.nc";
  expt.addTracer("salt");
  expt.seedParticles(new float[] {-122.5,-122.48,-122.46,-122.44,-122.42,-122.4}, new float[] {47.6, 47.65, 47.7, 47.75, 47.8}, new float[] {0}, expt.run.firstTime(), 1);
  expt.saveInterval = 1; // to save _every_ particle position to the output file
  for (int i=0; i<expt.particles.length; i++) {
    expt.particles[i].dt = 900; // this is how you override the default timestep in the Particle class for a particular experiment
  }
  expt.calcToTime(expt.run.lastTime());  
  
}


// ---------------------------------------------


void ecohab5_returnMap() {
  float lunarHour = 44712/12;
  float internalTimestep = lunarHour / 5.;
  ReturnMapMaker expt = new ReturnMapMaker("/Volumes/vindaloo3/salish_2005_1_summer/", 5800, 6202, "/Volumes/vindaloo3/ecohab5_maps/ecohab5", 2);
  expt.surfaceTrapped = true;
  float[] cslevels = {0};
  int Nreps = 1;
  float startTime = expt.run.lastTime() - 14*24*lunarHour; // 14 tidal days, roughly Sep 3-17, 2005
  boolean includeIndices = true;
  for (int n=0; n<14*24/2; n++) {
    println("map #" + n + "---------------");
    expt.makeMap(startTime + n*2*lunarHour, 2*lunarHour, cslevels, Nreps, internalTimestep, ""+n, includeIndices);
    includeIndices = false; // for all except the first
  }   
}

void ecohab5_returnMap_daily() {
  float lunarHour = 44712/12;
  float internalTimestep = lunarHour / 5.;
  ReturnMapMaker expt = new ReturnMapMaker("/Volumes/vindaloo3/salish_2005_1_summer/", 5800, 6202, "/Volumes/vindaloo3/ecohab5_maps/ecohab5daily", 2);
  expt.surfaceTrapped = true;
  float[] cslevels = {0};
  int Nreps = 1;
  float startTime = expt.run.lastTime() - 14*24*lunarHour; // 14 tidal days, roughly Sep 3-17, 2005
  boolean includeIndices = true;
  for (int n=0; n<14; n++) {
    println("map #" + n + "---------------");
    expt.makeMap(startTime + n*24*lunarHour, 24*lunarHour, cslevels, Nreps, internalTimestep, ""+n, includeIndices);
    includeIndices = false; // for all except the first
  }   
}


// ------------------------------------------------------------------------------------------------------------------------


void returnmaps2005(Configuration config) {
  String runDir = config.getString("runDir");
  String outputDir = config.getString("outputDir");
  String outputPrefix = config.getString("outputPrefix");
  int fileStart = config.getInt("fileStart");
  int fileEnd = config.getInt("fileEnd");  
  float cs = config.getFloat("cs");
  float lunarHour = 44712/12;
  float internalTimestep = lunarHour * config.getFloat("internalTimestep_lunarHours");
  float mapTimestep = lunarHour * config.getFloat("mapTimestep_lunarHours");
  int numMaps = config.getInt("numMaps");
//  float internalTimestep = lunarHour / 5.;
//  float mapTimestep = lunarHour * 3;
  ReturnMapMaker expt = new ReturnMapMaker(runDir, fileStart, fileEnd, outputDir+outputPrefix, 2);
  expt.surfaceTrapped = (cs == 0);
  int Nreps = 1;
  float startTime = expt.run.firstTime();
  boolean includeIndices = true;
  for (int n=0; n<numMaps; n++) {
    println("map #" + n + "---------------");
    expt.makeMap(startTime + n*mapTimestep, mapTimestep, new float[] {cs}, Nreps, internalTimestep, ""+n, includeIndices);
    includeIndices = false; // for all except the first
  }   
}


// ------------------------------------------


void jdf2005(Configuration config) {
  String runDir = config.getString("runDir");
  String outputDir = config.getString("outputDir");
  int fileStart = config.getInt("fileStart");
  int fileEnd = config.getInt("fileEnd");  
  float asp = cos(48.55/180*PI);
  float[] x = new float[25];
  for (int i=-12; i<=12; i++) x[i+12] = -125.25 + (2./111.325/asp)*i; // 2 km spacing, for 50 km centered on 125.25
  float[] y = new float[25];
  for (int i=-12; i<=12; i++) y[i+12] = 48.55 + (2./111.325)*i; // 2 km spacing, for 50 km centered on 48.55
  ParticleExpt expt = new ParticleExpt();
  expt.linkToRun(runDir+"ocean_his_",fileStart,fileEnd);
  for (int yearday = 365*5/12; yearday < 365*11/12; yearday += 2) {
    expt.ncname = outputDir + "jdf_particles_2005_v1." + yearday + ".nc";
    expt.seedParticles(x, y, new float[] {0}, yearday*86400, 1);
    expt.saveInterval = 18;
    for (int i=0; i<expt.particles.length; i++) {
      expt.particles[i].dt = 1200;
    }
    expt.calcToTime((yearday+30)*86400);
  }
}


// -------------------------------------------


void riverYear(Configuration config) {
  String runDir = config.getString("runDir");
  String outputDir = config.getString("outputDir");
  int fileStart = config.getInt("fileStart");
  int fileEnd = config.getInt("fileEnd");  
  ParticleExpt expt = new ParticleExpt();
  expt.linkToRun(runDir+"ocean_his_",fileStart,fileEnd);
  // from ps_2006_riverFile_salish.mat (rivers.rpos)
  float[] x = {-122.2884, -122.0906, -122.2705, -122.5762, -122.4038, -122.2976, -122.6664, -122.8686, -122.3513, -123.1472, -122.9888, -122.6978, -122.9414, -123.0631, -122.3656, -122.4697, -122.4255};
  float[] y = {  48.3875,   47.9977,   48.2078,   46.2604,   47.2613,   47.5748,   47.0905,   47.0420,   48.3856,   47.3348,   47.6538,   49.1385,   47.6996,   47.5523,   47.6656,   48.7910,   48.5526};
  float[] t = new float[365];
  for (int n=0; n<t.length; n++) t[n] = 86400*n; // seconds since start of year; new particle every day
  expt.dt = 400;
  expt.saveInterval = 27; // 27 * 400 sec = every 3 h
  for (int r=0; r<x.length; r++) { // one particle expt/file for each river
    expt.ncname = outputDir + "river_particles_2006." + r + ".nc";
    expt.seedParticles(new float[] {x[r]}, new float[] {y[r]}, new float[] {0}, t, 4);
    expt.calcToTime(t[t.length-1]);
  }
}

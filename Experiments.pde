/* a collection of particular applications.

SURFACEBOX is the most general of these: if you want to release particles in a rectangular area at a single time,
you don't need to modify the code, just make a custom particular_config file with your settings. This would also
make a good, general template if you want to write your own experiment.

RIVERYEAR does three things differently from surfaceBox: first, the particles are released at an "arbitrary" list
of locations, rather than a regular grid. Second, they are started at multiple time points. Third, the output is
saved to a series of files in increments, rather than all in one. This is a complicated experiment to use as a
template if you're new to particulator.

QMHYEAR is similar to riverYear, but contains an example of how to read a list of lat-lon locations from a text file.

JDF2005 is an example of how to release particles from a series of start times to a series of end times (tracking particles
for a fixed duration, rather than fixing the end date). This is a bit awkward at the moment.

*/



void pugetSoundTest() {
  ParticleRelease release = new ParticleRelease();
  release.linkToRun("/Users/neil/pnwtox/psmid_jul8_18/ocean_his_", 1633, 1873);
  release.ncname = "/Users/neil/Desktop/particles_ps.nc";
  release.addTracer("salt");
  release.seedParticles(new float[] {-122.5,-122.48,-122.46,-122.44,-122.42,-122.4}, new float[] {47.6, 47.65, 47.7, 47.75, 47.8}, new float[] {0}, release.run.firstTime(), 1);
  release.saveInterval = 1; // to save _every_ particle position to the output file
  for (int i=0; i<release.particles.length; i++) {
    release.particles[i].dt = 900; // this is how you override the default timestep in the Particle class for a particular experiment
  }
  release.calcToTime(release.run.lastTime());  
}


// ---------------------------------------------


void ecohab5_returnMap() {
  float lunarHour = 44712/12;
  float internalTimestep = lunarHour / 5.;
  ReturnMapMaker release = new ReturnMapMaker("/Volumes/vindaloo3/salish_2005_1_summer/", 5800, 6202, "/Volumes/vindaloo3/ecohab5_maps/ecohab5", 2);
  release.surfaceTrap();
  float[] cslevels = {0};
  int Nreps = 1;
  float startTime = release.run.lastTime() - 14*24*lunarHour; // 14 tidal days, roughly Sep 3-17, 2005
  boolean includeIndices = true;
  for (int n=0; n<14*24/2; n++) {
    println("map #" + n + "---------------");
    release.makeMap(startTime + n*2*lunarHour, 2*lunarHour, cslevels, Nreps, internalTimestep, ""+n);
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
  ReturnMapMaker release = new ReturnMapMaker(runDir, fileStart, fileEnd, outputDir+outputPrefix, 2);
  release.trapToSigmaLevel(cs);
  int Nreps = 1;
  float startTime = release.run.firstTime();
  for (int n=0; n<numMaps; n++) {
    println("map #" + n + "---------------");
    release.makeMap(startTime + n*mapTimestep, mapTimestep, new float[] {cs}, Nreps, internalTimestep, ""+n);
  }   
}


// ------------------------------------------


void jdf2005(Configuration config) {
  String runDir = config.getString("runDir");
  String outputDir = config.getString("outputDir");
  int fileStart = config.getInt("fileStart");
  int fileEnd = config.getInt("fileEnd");  
  int Nreps = config.getInt("Nreps");
  float asp = cos(48.55/180*PI);
  float[] x = new float[25];
  for (int i=-12; i<=12; i++) x[i+12] = -125.25 + (2./111.325/asp)*i; // 2 km spacing, for 50 km centered on 125.25
  float[] y = new float[25];
  for (int i=-12; i<=12; i++) y[i+12] = 48.55 + (2./111.325)*i; // 2 km spacing, for 50 km centered on 48.55
  ParticleRelease release = new ParticleRelease();
  release.linkToRun(runDir+"ocean_his_",fileStart,fileEnd);
  release.dt = 1200;
  release.saveInterval = 18;
  // since we want the end time to vary with relase time, we can't just put a vector of start times into seedParticles():
  // instead, step through start times, treating each start time as a new release and saving it to its own file 
  for (int yearday = 365*5/12; yearday < 365*11/12; yearday += 2) {
    release.ncname = outputDir + "jdf_particles_2005_v1." + yearday + ".nc";
    release.seedParticles(x, y, 0, yearday*86400, Nreps);
    release.calcToTime((yearday+30)*86400);
  }
}


// -------------------------------------------


void riverYear(Configuration config) {
  String runDir = config.getString("runDir");
  String outputBasename = config.getString("outputBasename");
  int fileStart = config.getInt("fileStart");
  int fileEnd = config.getInt("fileEnd");
  int Nreps = config.getInt("Nreps");
  float releaseInterval_hours = config.getFloat("releaseInterval_hours");
  float releaseOffset_hours = config.getFloat("releaseOffset_hours");
  ParticleRelease release = new ParticleRelease();
  release.multithread = true;
  release.linkToRun(runDir+"ocean_his_",fileStart,fileEnd);
  // mouths of all rivers in salish grid except Columbia; updated may 2012
  float[] x = {-122.4629,-122.2084,-122.3611,-122.4124,-122.3465,-122.7045,-122.8992,-122.3766,-123.1271,-122.9317,-123.1832,-122.8925,-123.0375,-122.4094,-122.5590,-122.4654};
  float[] y = {  48.3552,  48.0172,  48.2033,  47.2615,  47.5863,  47.0994,  47.0509,  48.3026,  47.3434,  47.6456,  49.1154,  47.6906,  47.5514,  47.6733,  48.7804,  48.5627};
  // vector of start times
  float[] t = new float[floor((release.run.lastTime() - release.run.firstTime()) / (releaseInterval_hours*3600))];
  for (int n=0; n<t.length; n++) t[n] = release.run.firstTime() + releaseOffset_hours * 3600 + releaseInterval_hours * 3600 * n;
  release.dt = 400;
  release.saveInterval = 27; // 27 * 400 sec = every 3 h. Ignored if multithreading!
  release.seedParticles("lonLatList", x, y, 0, t, Nreps);
  for (int i=1; i<=365; i++) {
    for (int j=0; j<release.particles.length; j++) release.particles[j].step = 0;
    release.ncname = outputBasename + "_day" + i + ".nc";
    release.calcToTime(t[0] + 86400*i); // save every day in its own file
  }
}


void cystYear(Configuration config) {
  String runDir = config.getString("runDir");
  String outputBasename = config.getString("outputBasename");
  int fileStart = config.getInt("fileStart");
  int fileEnd = config.getInt("fileEnd");
  int Nreps = config.getInt("Nreps");
  float releaseInterval_hours = config.getFloat("releaseInterval_hours");
  float releaseOffset_hours = config.getFloat("releaseOffset_hours");
  ParticleRelease release = new ParticleRelease();
  release.multithread = true;
  release.linkToRun(runDir+"ocean_his_",fileStart,fileEnd);
  // 4 pts in QMH, 1 in Liberty Bay, 2 in Bellingham Bay, 3 in Sequim Bay
  float[] x = {-122.452, -122.442, -122.449, -122.471,    -122.641,    -122.612, -122.504,    -123.006, -123.032, -123.012};
  float[] y = {  47.400,   47.391,   47.382,   47.377,      47.721,      48.733,   48.681,      48.035,   48.057,   48.067};
  float[] t = new float[floor((release.run.lastTime() - release.run.firstTime()) / (releaseInterval_hours*3600))];
  for (int n=0; n<t.length; n++) t[n] = release.run.firstTime() + releaseOffset_hours * 3600 + releaseInterval_hours * 3600 * n;
  release.dt = 400;
  release.saveInterval = 27; // 27 * 400 sec = every 3 h
  release.seedParticles("lonLatList", x, y, 0, t, Nreps);
  for (int i=1; i<=365; i++) {
    for (int j=0; j<release.particles.length; j++) release.particles[j].step = 0;
    release.ncname = outputBasename + "_day" + i + ".nc";
    release.calcToTime(t[0] + 86400*i); // save every day in its own file
  }
}


void qmhYear(Configuration config) {
  String runDir = config.getString("runDir");
  String outputBasename = config.getString("outputBasename");
  int fileStart = config.getInt("fileStart");
  int fileEnd = config.getInt("fileEnd");
  int Nreps = config.getInt("Nreps");
  float releaseInterval_hours = config.getFloat("releaseInterval_hours");
  float releaseOffset_hours = config.getFloat("releaseOffset_hours");
  ParticleRelease release = new ParticleRelease();
  release.multithread = true;
  release.linkToRun(runDir+"ocean_his_",fileStart,fileEnd);
  String[] s = loadStrings("qmh.txt"); // this is how to read x and y from an external file...
  float[] x = new float[s.length];
  float[] y = new float[s.length];
  for (int i=0; i<s.length; i++) {
    float[] xy = float(splitTokens(s[i]));
    x[i] = xy[0];
    y[i] = xy[1];
  }
  float[] cs = {-5./6, -0.5, -1./6};
  float[] t = new float[floor((release.run.lastTime() - release.run.firstTime()) / (releaseInterval_hours*3600))];
  for (int n=0; n<t.length; n++) t[n] = release.run.firstTime() + releaseOffset_hours * 3600 + releaseInterval_hours * 3600 * n;
  release.dt = 400;
  release.saveInterval = 27; // 27 * 400 sec = every 3 h
  release.seedParticles("lonLatList", x, y, cs, t, Nreps);
  for (int i=1; i<=365; i++) {
    for (int j=0; j<release.particles.length; j++) release.particles[j].step = 0;
    release.ncname = outputBasename + "_day" + i + ".nc";
    release.calcToTime(t[0] + 86400*i); // save every day in its own file
  }
}


// ---------------------------------------------------------------------
void surfaceBox(Configuration config) {
  String runDir = config.getString("runDir");
  String outputDir = config.getString("outputDir");
  String outputBasename = config.getString("outputBasename");
  int fileStart = config.getInt("fileStart");
  int fileEnd = config.getInt("fileEnd");
  float startTime_yearday = config.getFloat("startTime_yearday");
  float duration_days = config.getFloat("duration_days");
  float latMin = config.getFloat("latMin");
  float latMax = config.getFloat("latMax");
  float lonMin = config.getFloat("lonMin");
  float lonMax = config.getFloat("lonMax");
  float spacing_km = config.getFloat("spacing_km");
  boolean doTS = config.getBoolean("doTS");
  boolean doBio = config.getBoolean("doBio");
  boolean do3D = config.getBoolean("do3D");
  boolean doSurface = config.getBoolean("doSurface"); 
  float asp = cos((latMin+latMax)/2/180*PI);
  float[] x = series(lonMin, lonMax, spacing_km/111.325/asp);
  float[] y = series(latMin, latMax, spacing_km/111.325);
  ParticleRelease release = new ParticleRelease();
  release.linkToRun(runDir+"ocean_his_",fileStart,fileEnd);
  release.dt = 1200;
  release.saveInterval = 18; // ignored if multithreading!
  if (doTS) {
    release.addTracer("salt");
    release.addTracer("temp");
  }
  if (doBio) {
    release.addTracer("NO3");
    release.addTracer("phytoplankton");
    release.addTracer("zooplankton");
  }
  if (do3D) {
    release.ncname = outputDir + outputBasename + ".3D.nc";
    release.seedParticles(x, y, 0, startTime_yearday*86400, 1);
    release.calcToTime((startTime_yearday+duration_days)*86400);
  }  
  if (doSurface) {
    release.ncname = outputDir + outputBasename + ".surface.nc";
    release.seedParticles(x, y, 0, startTime_yearday*86400, 1);
    release.surfaceTrap();
    release.calcToTime((startTime_yearday+duration_days)*86400);
  }
}

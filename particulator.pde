// particulator 3.0a
// neil banas, jun 2012

/* new in this version
  multithreaded!
  note that since the old autoSave method is apparently incompatible with multithreading, when multithreading is on,
  autoSave simply saves once per _frame of input data_, not particle timestep * saveInterval. The options for this
  need to be rethought.
*/

/* to do
  - save a simple flag that indicates which steps in the output file should be ignored (for example, particles that haven't started yet and thus are frozen at their start position)
  - flexible endTime vector that goes along with startTime
  - Particle.seedParticlesFromHotstart()
  - replacement for saveInterval that works with multithreading
  - Release.createNetcdf() for a subset of particles (pass the indices to createNetcdf(), or else 'all', which would also serve as a reminder that particles[] has to exist before this is called)
  - test higher-order stepping schemes (start with AB3)
  - divide seedParticles into cs and z versions
  - less hokey way to determine dz in dKs/dz
  - do a proper convergence test in the coastal domain; Puget Sound is so chaotically dispersive that it's a hard place to look at this
*/

boolean debug = true;
boolean useConfigFile = true;


void setup() {
  float tic = millis();
  if (useConfigFile) {
    Configuration config = new Configuration("particulator_config.txt"); // mian program always looks for a config file named this
    String exptname = config.getString("exptname");
    if (debug) config.echo();
    
    if (exptname.equals("surfaceBox")) surfaceBox(config);
    if (exptname.equals("riverYear")) riverYear(config);
    if (exptname.equals("riverYearHotstart")) riverYearHotstart(config);
    if (exptname.equals("cystYear")) cystYear(config);
    if (exptname.equals("qmhYear")) qmhYear(config);
    // etc etc etc
    
  } else {
    println("don't know what do to if not using a config file, although could be doing something...");
  }
  
  if (debug) println("elapsed time " + (millis()-tic)/1000. + " sec");
  exit();
}

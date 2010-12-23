class Configuration extends HashMap {
  
  Configuration(String filename) {
    String[] lines = loadStrings(filename);
    if (lines == null) println("configuration file " + filename + " not found.");
    for (int i=0; i<lines.length; i++) {
      String line = trim(lines[i]);
      int commentStart = line.indexOf('!');
      if (commentStart != 0) {
        if (commentStart > -1) line = line.substring(0,commentStart);
        String[] elements = split(line, '=');
        if (elements.length == 2) {
          put(trim(elements[0]), trim(elements[1]));
        }
      }
    }
  }
  
  void echo() {
    Iterator i = keySet().iterator();
    while(i.hasNext()) {
      String key = (String)i.next();
      println(key + " = " + getString(key));
    }
  }
  
  int getInt(String key) {
    return int((String)get(key));
  }
  
  float getFloat(String key) {
    return float((String)get(key));
  }
  
  boolean getBoolean(String key) {
    return getString(key).equals("true");
  }
  
  String getString(String key) {
    return (String)get(key);
  }
  
}

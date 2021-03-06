/**
 * TestedResultFile.java
 *
 * Parses and stores data from a ".tested" file generated by the search script.
 *
 * Original author:
 *   Mike Lam (lam@cs.umd.edu)
 *   Professor Jeffrey K. Hollingsworth, UMD
 *   February 2013
 */

// {{{ imports
import java.io.*;
import java.text.*;
import java.util.*;
import java.util.regex.*;
import javax.swing.*;
// }}}

public class TestedResultFile {

    public List<TestedResult> allResults;

    public TestedResultFile(File rfile) {
        allResults = new ArrayList<TestedResult>();
        parseResultFile(rfile);
    }

    public void parseResultFile(File rfile) {
        try {
            BufferedReader fin = new BufferedReader(new FileReader(rfile));
            String line = null;
            TestedResult latestResult = null;

            while ((line = fin.readLine()) != null) {
                line = line.trim();
                if (line.equals("--- !ruby/object:AppConfig")) {
                    latestResult = new TestedResult();
                    allResults.add(latestResult);
                } else if (line.startsWith("cuid:")) {
                    latestResult.cuid = Util.extractRegex(line, "cuid: (.*)$", 1);
                } else if (line.startsWith("default:")) {
                    latestResult.default_cfg = Util.extractRegex(line, "default: (.*)$", 1);
                } else if (line.startsWith("level:")) {
                    latestResult.level = Util.extractRegex(line, "level: (.*)$", 1);
                } else if (line.startsWith("result:")) {
                    latestResult.result = Util.extractRegex(line, "result: (.*)$", 1);
                } else if (line.startsWith("runtime:")) {
                    latestResult.runtime = Long.parseLong(Util.extractRegex(line, "runtime: (.*)$", 1));
                } else if (line.startsWith("error:")) {
                    latestResult.error = Double.parseDouble(Util.extractRegex(line, "error: (.*)$", 1));
                } else if (line.startsWith("label:")) {
                    latestResult.label = Util.extractRegex(line, "label: (.*)$", 1);
                } else if (line.matches("^\"(.*)\": (.)$")) {
                    String key = Util.extractRegex(line, "^\"(.*)\": (.)$", 1);
                    String value = Util.extractRegex(line, "^\"(.*)\": (.)$", 2);
                    latestResult.exceptions.put(key, value);
                }
            }

            fin.close();
        } catch (IOException e) {
            JOptionPane.showMessageDialog(null, "I/O error: " + e.getMessage());
        }
    }

    public List<TestedResult> getAllResults() {
        return allResults;
    }

    public List<TestedResult> getResult(String regex) {
        List<TestedResult> results = new ArrayList<TestedResult>();
        // TODO: search
        return results;
    }
}


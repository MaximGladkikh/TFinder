package ru.spbau.ablab.tagfinder.util;

import ru.spbau.ablab.tagfinder.util.io.FastScanner;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.HashMap;


public class ConfigReader {
	private static final String CONFIG_FILE_PATH = "stat.cfg";
	private static final HashMap<String, String> PROPERTY_MAP = new HashMap<String, String>();
	static {
		FastScanner scanner;
		try {
			scanner = new FastScanner(new File(CONFIG_FILE_PATH));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
			throw new RuntimeException("Can't read config file " + CONFIG_FILE_PATH);
		}
		for (String s; (s = scanner.nextLine()) != null;) {
			String[] ss = s.split("=");
			int index = ss[1].indexOf("//");
			if (index >= 0) {
				ss[1] = ss[1].substring(0, index);
			}
			PROPERTY_MAP.put(ss[0], ss[1]);
		}
		scanner.close();
	}
	
	public static int getIntProperty(String key) {
		return Integer.parseInt(getProperty(key));
	}
	
	public static boolean getBooleanProperty(String key) {
		return Boolean.parseBoolean(getProperty(key));
	}
	
	public static double getDoubleProperty(String key) {
		return Double.parseDouble(getProperty(key));
	}
	
	public static long getLongProperty(String key) {
		return Long.parseLong(getProperty(key));
	}

	public static String getProperty(String key) {
		String s = PROPERTY_MAP.get(key);
		if (s == null) {
			throw new RuntimeException("No property " + key + " in config file");
		}
		return s;
	}

    public static void setProperty(String key, String value) {
        PROPERTY_MAP.put(key, value);
    }
}

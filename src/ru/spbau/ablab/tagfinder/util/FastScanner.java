package ru.spbau.ablab.tagfinder.util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.StringReader;
import java.util.StringTokenizer;

public class FastScanner {
	private BufferedReader reader;
	private StringTokenizer tokenizer;

	public FastScanner(File file) throws FileNotFoundException {
		reader = new BufferedReader(new FileReader(file));
	}

	public FastScanner(String s) {
		reader = new BufferedReader(new StringReader(s));
	}

	public void skipToken(String token) {
		while (!nextToken().equals(token))
			;
	}

	public boolean skipLine(String line) {
		try {
			String s;
			while ((s = reader.readLine()) != null && !line.equals(s.trim()))
				;
			return s != null;
		} catch (IOException e) {
			e.printStackTrace();
			return false;
		}
	}
	
	public double[] getDoubleArray() {
		int n = nextInt();
		double[] a = new double[n];
		for (int i = 0; i < n; ++i) {
			a[i] = nextDouble();
		}
		return a;
	}

	public double getNextDoubleProperty(String key) {
		return Double.parseDouble(skipToToken(key).split("=")[1]);
	}

	public int getNextIntProperty(String key) {
		return Integer.parseInt(skipToToken(key).split("=")[1]);
	}

	private String skipToToken(String key) {
		String s;
		do {
			s = nextToken();
		} while (!s.startsWith(key));
		return s;
	}

	public int getNextIntValue(String key) {
		skipToken(key);
		return nextInt();
	}

	public double getNextDoubleValue(String key) {
		skipToken(key);
		return nextDouble();
	}

	public String nextToken() {
		while (tokenizer == null || !tokenizer.hasMoreTokens()) {
			String line;
			try {
				line = reader.readLine();
			} catch (IOException e) {
				e.printStackTrace();
				return null;
			}
			if (line == null) {
				return null;
			}
			tokenizer = new StringTokenizer(line);
		}
		return tokenizer.nextToken();
	}

	public String nextLine() {
		try {
			return reader.readLine();
		} catch (IOException e) {
			e.printStackTrace();
			return null;
		}
	}

	public boolean hasNextLine() {
		try {
			String s = reader.readLine();
			if (s != null) {
				tokenizer = new StringTokenizer(s);
			}
			return s != null;
		} catch (IOException e) {
			e.printStackTrace();
			return false;
		}
	}

	public int nextInt() {
		return Integer.parseInt(nextToken());
	}

	public double nextDouble() {
		return Double.parseDouble(nextToken());
	}

	public long nextLong() {
		return Long.parseLong(nextToken());
	}

	public void close() {
		try {
			reader.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}

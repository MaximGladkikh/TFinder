package ru.spbau.ablab.tagfinder;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Collections;

import ru.spbau.ablab.tagfinder.util.FastScanner;

public class TableReader {
	class Spectrum implements Comparable<Spectrum>{
		int id;
		int redDiff;
		int envelopes;
		public Spectrum(int id, int redDiff, int envelopes) {
			super();
			this.id = id;
			this.redDiff = redDiff;
			this.envelopes = envelopes;
		}
		@Override
		public int compareTo(Spectrum o) {
			int d = -(redDiff - o.redDiff);
			if (d == 0) {
				d = envelopes - o.envelopes;
				if (d == 0) {
					return id - o.id;
				}
				return d;
			}
			return d;
		}
	}
	public static void main(String[] args) throws FileNotFoundException {
		int[][] expData = new int[103][4];
		int[][] virtData = new int[103][4];
		FastScanner exp = new FastScanner(new File("table1.in"));
		FastScanner virt = new FastScanner(new File("table2.in"));
		for (int i =0 ; i < 103; ++i) {
			for (int j = 0; j < 4; ++j) {
				expData[i][j] = exp.nextInt();
				virtData[i][j] = virt.nextInt();
			}
		}
		exp.close();
		virt.close();
		ArrayList<Spectrum> list = new ArrayList<TableReader.Spectrum>();
		for (int i = 0; i < 103; ++i) {
			if (virtData[i][2] == 10 && expData[i][2] == 10) {
//				list.add(new Spectrum(virtData[i][0], virtData[i][3] - expData[i][3], virtData[i][1]));
			}
		}
		Collections.sort(list);
		for (Spectrum spectrum : list) {
			System.err.println(spectrum.id);
		}
	}
}

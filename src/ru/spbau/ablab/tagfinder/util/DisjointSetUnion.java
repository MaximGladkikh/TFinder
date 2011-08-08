package ru.spbau.ablab.tagfinder.util;

public class DisjointSetUnion {
	private int[] parent;
	
	public DisjointSetUnion(int n) {
		parent = new int[n];
		for (int i = 0; i < n; ++i) {
			parent[i] = i;
		}
	}
	
	public void unite(int v, int u) {
		v = findSet(v);
		u = findSet(u);
		parent[v] = u;
	}

	public int findSet(int v) {
		if (parent[v] == v) {
			return v;
		}
		return parent[v] = findSet(parent[v]);
	}
}

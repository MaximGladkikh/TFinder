package ru.spbau.ablab.tagfinder.path;

public interface Edge extends Comparable<Edge>{
	double getMass();
	Character getLetter();
}

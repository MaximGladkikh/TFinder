package ru.spbau.ablab.tagfinder;

import ru.spbau.ablab.tagfinder.util.ConfigReader;

public class StatRunner {
	public static void main(String[] args) {
		if (ConfigReader.getBoolProperty("ALIGN")) {
			new StatisticsGeneratorVirtual().run();
		} else {
			new StatisticsGeneratorExperimental().run();
		}
	}
}

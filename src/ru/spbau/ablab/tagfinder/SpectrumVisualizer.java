package ru.spbau.ablab.tagfinder;

import java.awt.Canvas;
import java.awt.Color;
import java.awt.Graphics;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.event.MouseEvent;
import java.awt.event.MouseMotionAdapter;
import java.io.File;
import java.io.FileNotFoundException;

import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.KeyStroke;

import ru.spbau.ablab.tagfinder.spectrum.Envelope;
import ru.spbau.ablab.tagfinder.spectrum.Spectrum;
import ru.spbau.ablab.tagfinder.util.FastScanner;

public class SpectrumVisualizer extends JFrame {
    private static final String envFile = StatisticsGenerator.SPECTRUM_FILE_SUFFIX.substring(1);
    private static final String FRAME_TITLE = "Spectrum Vizualizer";
    private static final long serialVersionUID = 1L;

    private final JFrame mainFrame = this;
    private final Painter painter = new Painter();

    public class Painter extends Canvas {
        private static final long serialVersionUID = 1L;
        private static final int BOX_SIZE = 5;

        private Spectrum spectrum;
        private Spectrum virtualSpectrum;
        private double maxInt;
        private double maxMass;
        private Envelope selectedPeak;
        private double coeff = 0.95;

        public Painter() {
            addMouseMotionListener(new MouseMotionAdapter() {
                @Override
                public void mouseMoved(MouseEvent e) {
                    if (spectrum == null) {
                        return;
                    }
                    int mx = e.getX();
                    int my = e.getY();
                    Envelope before = selectedPeak;
                    selectedPeak = null;
                    updateSelection(mx, my, spectrum);
                    if (virtualSpectrum != null) {
                        updateSelection(mx, my, virtualSpectrum);
                    }
                    if (before != selectedPeak) {
                        repaint();
                    }
                }

                private void updateSelection(int mx, int my, Spectrum spectrum) {
                    for (Envelope envelope : spectrum.envelopes) {
                        int xl = getPeakX(getXCoeff(), envelope);
                        int yl = getPeakY(getYCoeff(), envelope);
                        int xr = xl + BOX_SIZE;
                        int yr = yl + BOX_SIZE;
                        if (mx >= xl && mx <= xr && my >= yl && my <= yr) {
                            selectedPeak = envelope;
                            break;
                        }
                    }
                }
            });
        }

        @Override
        public void paint(Graphics g) {
            super.paint(g);
            if (spectrum == null) {
                return;
            }

            double cy = getYCoeff();
            double cx = getXCoeff();

            paintSpectrum(spectrum, g, cy, cx, Color.RED);
            if (virtualSpectrum != null) {
                g.setXORMode(Color.BLACK);
                paintSpectrum(virtualSpectrum, g, cy, cx, Color.MAGENTA);
                g.setPaintMode();
            }

            g.setColor(Color.BLACK);
            int fontHeight = g.getFontMetrics().getHeight();
            g.drawString("" + maxInt, 0, fontHeight);
            String massString = "" + maxMass;
            g.drawString(massString,
                    getWidth() - g.getFontMetrics().charsWidth(Double.toString(maxMass).toCharArray(), 0, massString.length()), getHeight());

            if (selectedPeak != null) {
                String message = selectedPeak.mass + " " + selectedPeak.intensity;
                int width = g.getFontMetrics().charsWidth(message.toCharArray(), 0, message.length());
                int height = g.getFontMetrics().getHeight();
                int xl = Math.min(getPeakX(cx, selectedPeak) + BOX_SIZE, getWidth() - width);
                int yl = getPeakY(cy, selectedPeak) - g.getFontMetrics().getHeight();
                g.setColor(Color.WHITE);
                g.fillRect(xl, yl, width, height);
                g.setColor(Color.BLACK);
                g.drawRect(xl, yl, width, height);
                g.drawString(message, xl, yl + height);
            }
        }

        private void paintSpectrum(Spectrum spectrum, Graphics g, double cy, double cx, Color color) throws AssertionError {
            g.setColor(color);
            for (Envelope envelope : spectrum.envelopes) {
                int y = getPeakY(cy, envelope);
                int x = getPeakX(cx, envelope);
                g.fillRect(x, y, BOX_SIZE, BOX_SIZE);
            }
        }

        private int getPeakX(double cx, Envelope envelope) {
            return (int) (cx * envelope.mass);
        }

        private int getPeakY(double cy, Envelope envelope) {
            return getHeight() - 1 - (int) (cy * envelope.intensity);
        }

        private double getXCoeff() {
            return 0.95 * getWidth() / maxMass;
        }

        private double getYCoeff() {
            return coeff * getHeight() / maxInt;
        }

        @Override
        public void update(Graphics g) {
            super.update(g);
            paint(g);
        }

        public void setSpectrum(Spectrum spectrum) {
            fillInfinity();
            init(spectrum);
        }

        private void fillInfinity() {
            maxInt = Double.NEGATIVE_INFINITY;
            maxMass = Double.NEGATIVE_INFINITY;
        }

        private void init(Spectrum spectrum) {
            this.spectrum = spectrum;
            relaxMaxValues(spectrum);
            selectedPeak = null;
            repaint();
        }

        private void relaxMaxValues(Spectrum spectrum) {
            for (Envelope envelope : spectrum.envelopes) {
                maxInt = Math.max(maxInt, envelope.intensity);
                maxMass = Math.max(maxMass, envelope.mass);
            }
        }

        public void setSpectrum(Spectrum spectrum, Spectrum virtualSpectrum) {
            this.virtualSpectrum = virtualSpectrum;
            fillInfinity();
            init(spectrum);
            relaxMaxValues(virtualSpectrum);
        }
    }

    public SpectrumVisualizer() {
        JMenuBar menuBar = new JMenuBar();
        JMenu menu = new JMenu("File");
        menu.setMnemonic('f');
        JMenuItem item = new JMenuItem("Open");
        item.setMnemonic('o');
        item.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_O, KeyEvent.CTRL_MASK));
        menu.add(item);
        menuBar.add(menu);
        setJMenuBar(menuBar);
        add(painter);
        item.setActionCommand("open");
        ActionListener listener = new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if (e.getActionCommand().equals("open")) {
                    JFileChooser fileChooser = new JFileChooser(new File(".").getAbsolutePath());
                    fileChooser.setFileFilter(new javax.swing.filechooser.FileFilter() {
                        @Override
                        public boolean accept(File f) {
                            return f.isDirectory() || f.getName().endsWith("." + envFile);
                        }

                        @Override
                        public String getDescription() {
                            return envFile;
                        }
                    });
                    if (fileChooser.showOpenDialog(mainFrame) == JFileChooser.APPROVE_OPTION) {
                        File file = fileChooser.getSelectedFile();
                        String path = file.getAbsolutePath();
                        String str = StatisticsGenerator.SPECTRUM_FILE_SUFFIX;
                        if (!path.contains(str)) {
                            return;
                        }
                        File virtual = new File(path.substring(0, path.indexOf(str)) + ".vir");
                        try {
                            FastScanner fastScanner = new FastScanner(virtual);
                            double[] masses = fastScanner.getDoubleArray();
                            double[] intensities = fastScanner.getDoubleArray();
                            Envelope[] virEnv = new Envelope[masses.length];
                            for (int i = 0; i < virEnv.length; ++i) {
                                virEnv[i] = new Envelope(masses[i], intensities[i], intensities[i]);
                            }
                            painter.setSpectrum(new Spectrum(0, file, 0), new Spectrum(0, virEnv, 0, 0));
                        } catch (FileNotFoundException e1) {
                            e1.printStackTrace();
                        }
                        setTitle(file.getName() + " - " + FRAME_TITLE);
                    }
                } else if (e.getActionCommand().equals("diff")) {
                    JFrame diffFrame = new JFrame("Difference");
                    int count = 0;
                    String[][] tableData = new String[painter.spectrum.envelopes.length][3];
                    for (Envelope envelope : painter.spectrum.envelopes) {
                        Envelope nearest = null;
                        for (Envelope envelope2 : painter.virtualSpectrum.envelopes) {
                            if (nearest == null || Math.abs(envelope.mass - envelope2.mass) < Math.abs(envelope.mass - nearest.mass)) {
                                nearest = envelope2;
                            }
                        }
                        tableData[count][0] = "" + envelope.mass;
                        assert nearest != null;
                        tableData[count][1] = "" + nearest.mass;
                        tableData[count][2] = "" + Math.abs(envelope.mass - nearest.mass);
                        ++count;
                    }
                    JTable table = new JTable(tableData, new String[]{"Virtual Peaks", "Closest Experimental Peaks", "Mass Difference"});
                    JScrollPane jScrollPane = new JScrollPane(table);
                    diffFrame.add(jScrollPane);
                    diffFrame.setSize((int) diffFrame.getPreferredSize().getWidth(), (int) diffFrame.getPreferredSize().getHeight() + 40);
                    diffFrame.setVisible(true);
                } else if (e.getActionCommand().equals("+")) {
                    painter.coeff *= 2;
                    painter.repaint();
                } else if (e.getActionCommand().equals("-")) {
                    painter.coeff /= 2;
                    painter.repaint();
                }
            }
        };
        item.addActionListener(listener);

        JMenu viewMenu = new JMenu("View");
        menuBar.add(viewMenu);
        viewMenu.setMnemonic('v');
        JMenuItem incItem = new JMenuItem("Increase scale");
        viewMenu.add(incItem);
        incItem.setMnemonic('i');
        incItem.setAccelerator(KeyStroke.getKeyStroke('='));
        incItem.setActionCommand("+");
        incItem.addActionListener(listener);
        JMenuItem decItem = new JMenuItem("Decrease scale");
        viewMenu.add(decItem);
        decItem.setMnemonic('e');
        decItem.setAccelerator(KeyStroke.getKeyStroke('-'));
        decItem.setActionCommand("-");
        decItem.addActionListener(listener);

        JMenu toolsMenu = new JMenu("Tools");
        menuBar.add(toolsMenu);
        toolsMenu.setMnemonic('t');
        JMenuItem diffItem = new JMenuItem("Show difference");
        toolsMenu.add(diffItem);
        diffItem.setMnemonic('d');
        diffItem.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_D, KeyEvent.CTRL_MASK));
        diffItem.setActionCommand("diff");
        diffItem.addActionListener(listener);
    }

    public static void main(String[] args) {
        SpectrumVisualizer mainFrame = new SpectrumVisualizer();
        mainFrame.setTitle(FRAME_TITLE);
        mainFrame.setSize(800, 600);
        mainFrame.setDefaultCloseOperation(EXIT_ON_CLOSE);
        mainFrame.setVisible(true);
    }
}

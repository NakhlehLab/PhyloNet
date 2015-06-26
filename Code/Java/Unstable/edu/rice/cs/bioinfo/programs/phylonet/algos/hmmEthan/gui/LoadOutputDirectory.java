package edu.rice.cs.bioinfo.programs.phylonet.algos.hmmEthan.gui;

import com.google.gson.*;
import edu.rice.cs.bioinfo.programs.phylonet.algos.hmmEthan.Configuration;
import edu.rice.cs.bioinfo.programs.phylonet.algos.hmmEthan.model.ProcessedOutputInformation;
import edu.rice.cs.bioinfo.programs.phylonet.algos.hmmEthan.model.RawOutputInformation;
import edu.rice.cs.bioinfo.programs.phylonet.algos.hmmEthan.model.FullOutputInformation;
import edu.rice.cs.bioinfo.programs.phylonet.algos.hmmEthan.util.JsonUtilities;

import javax.swing.*;
import java.nio.file.Path;

/**
 * This class shows a gui for loading a saved json file and creating a corresponding DialogWithTabs box.
 */
public class LoadOutputDirectory {

    private static Path selectDataFileToLoad() {
        JFileChooser chooser = new JFileChooser();
        chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
        int choice = chooser.showDialog(null, "Load");

        switch (choice) {
            case JFileChooser.CANCEL_OPTION:
                return null;

            case JFileChooser.ERROR_OPTION:
                throw new RuntimeException("ERROR in JFILECHOOSER");

            case JFileChooser.APPROVE_OPTION:
                return chooser.getSelectedFile().toPath();
        }
        throw new RuntimeException("JFileChooser choice not in correct range");
    }


    /**
     * The main method which provides a gui for loading saved data files.
     */
    public static void main(String[] args) throws ClassNotFoundException, UnsupportedLookAndFeelException, InstantiationException, IllegalAccessException {
        UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());

        Path pathToDirectory = selectDataFileToLoad();

        Gson g = new Gson();

        Path pathToRawData = pathToDirectory.resolve("rawOutput.json");
        RawOutputInformation info = JsonUtilities.readJson(pathToRawData,g,RawOutputInformation.class);

        Path pathToProcessedData = pathToDirectory.resolve("output.json");
        ProcessedOutputInformation process = JsonUtilities.readJson(pathToProcessedData,g,ProcessedOutputInformation.class);

        Path pathToConfig = pathToDirectory.resolveSibling("config.json");
        Configuration configuration = JsonUtilities.readJson(pathToConfig,g,Configuration.class);

        createTabs(info,process,configuration);
    }




    private static void createTabs(RawOutputInformation rawOutputInformation, ProcessedOutputInformation processedOutputInformation, Configuration config) {

        new PlotOutput(new FullOutputInformation(rawOutputInformation,processedOutputInformation),config).show();
    }
}

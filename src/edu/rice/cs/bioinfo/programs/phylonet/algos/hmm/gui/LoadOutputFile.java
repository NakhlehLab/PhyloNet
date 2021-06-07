package edu.rice.cs.bioinfo.programs.phylonet.algos.hmm.gui;

import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;

import javax.swing.*;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

/**
 * This class shows a gui for loading a saved json file and creating a corresponding DialogWithTabs box.
 */
public class LoadOutputFile {

    private static File selectDataFileToLoad() {
        JFileChooser chooser = new JFileChooser();
        int choice = chooser.showDialog(null, "Load");

        switch (choice) {
            case JFileChooser.CANCEL_OPTION:
                return null;

            case JFileChooser.ERROR_OPTION:
                throw new RuntimeException("ERROR in JFILECHOOSER");

            case JFileChooser.APPROVE_OPTION:
                return chooser.getSelectedFile();
        }
        throw new RuntimeException("Jfilechooser choice not in correct range");
    }


    /**
     * The main method which provides a gui for loading saved data files.
     */
    public static void main(String[] args) throws ClassNotFoundException, UnsupportedLookAndFeelException, InstantiationException, IllegalAccessException {
        UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());

        File jsonFile = selectDataFileToLoad();
        JsonElement json = loadFile(jsonFile);
        LoadOutputFile.createTabs(json);
    }


    /**
     * Creates a DialogWithTabs with the data from the provided json element.
     *
     * @param jsonElement A json element storing the data to load.
     */
    private static void createTabs(JsonElement jsonElement) {

        DialogWithTabs gui = new DialogWithTabs();

        JsonArray tabs = jsonElement.getAsJsonArray();

        for (JsonElement tab : tabs) {
            JsonObject obj = tab.getAsJsonObject();

            String name = obj.get("name").getAsString();
            JsonObject graph = obj.get("tab").getAsJsonObject();

            SavableTabData tabData = SavableTabData.load(graph);

            gui.addPanel(name, tabData);
        }
        gui.createAndShowGUI();
    }

    /**
     * Loads a parses a file into a JsonElement.
     *
     * @param file The file to read. The file will be closed after reading.
     * @return A JsonElement of the data in the file.
     */
    private static JsonElement loadFile(File file) {
        try (FileReader reader = new FileReader(file)) {
            return new JsonParser().parse(reader);
        } catch (IOException e) {
            e.printStackTrace();
            throw new RuntimeException(e);
        }
    }
}

package edu.rice.cs.bioinfo.programs.phylonet.algos.hmmEthan.gui;


import com.google.gson.JsonElement;
import com.google.gson.JsonObject;

import javax.swing.*;

/**
 * An interface for something that can both be showed in our DialogWithTabs and saved to a json structure.
 */
public abstract class SavableTabData {
    /**
     * deserialized the given JsonObject into a SavableTabData class.
     *
     * @param serialized The data to be deserialized
     * @return The SavebleTabData represented by the given json.
     */
    public static SavableTabData load(JsonObject serialized) {
        String type = serialized.get("type").getAsString();
        JsonObject data = serialized.get("data").getAsJsonObject();

        switch (type) {
            case "PlotAreaGraph":
                return new PlotAreaGraph(data);
            case "PlotBarGraph":
                return new PlotBarGraph(data);
            case "PlotLineGraph":
                return new PlotLineGraph(data);
            case "TextOutputTab":
                return new TextOutputTab(data);
        }

        throw new RuntimeException("Invalid type of " + type);
    }

    /**
     * Serializes the given object into json.
     *
     * @return The json representing the object.
     */
    public JsonElement save() {
        JsonObject obj = new JsonObject();
        obj.addProperty("type", getType());
        obj.add("data", getData());
        return obj;
    }

    /**
     * Shows the tab in a standalone frame.
     */
    public void show() {
        JFrame frame = new JFrame();
        frame.add(getComponent());
        frame.pack();
        frame.setVisible(true);
    }

    /**
     * Get a panel to display the current data.
     *
     * @return The gui panel.
     */
    public abstract JComponent getComponent();

    /**
     * Get which type of SavableTabData this is.
     *
     * @return A string representing the type.
     */
    public abstract String getType();

    /**
     * Get the associated data for a SavableTabObject.
     *
     * @return The associated data in json format.
     */
    public abstract JsonElement getData();
}

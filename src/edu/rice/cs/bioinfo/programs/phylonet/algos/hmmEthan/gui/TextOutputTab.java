package edu.rice.cs.bioinfo.programs.phylonet.algos.hmmEthan.gui;

import com.google.gson.JsonElement;
import com.google.gson.JsonObject;

import javax.swing.*;

public class TextOutputTab extends SavableTabData
{
    String data;
    public TextOutputTab(String inputdata){
        data = inputdata;
    }

    public TextOutputTab(JsonObject JObj){
        data = JObj.get("thedata").getAsString();
    }

    @Override
    public JComponent getComponent() {
        JTextArea textTabArea = new JTextArea(data);
        textTabArea.setEditable(false);
        JScrollPane scrollingTextTabArea = new JScrollPane(textTabArea);
        return scrollingTextTabArea;
    }

    @Override
    public String getType() {
        return "TextOutputTab";
    }

    @Override
    public JsonElement getData() {

        JsonObject obj = new JsonObject();
        obj.addProperty("thedata", data);
        return obj;
    }
}

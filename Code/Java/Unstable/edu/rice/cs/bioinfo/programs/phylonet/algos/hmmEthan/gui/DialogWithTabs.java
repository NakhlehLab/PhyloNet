package edu.rice.cs.bioinfo.programs.phylonet.algos.hmmEthan.gui;

import com.google.gson.JsonArray;
import com.google.gson.JsonObject;
import com.google.gson.JsonPrimitive;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;


public class DialogWithTabs extends JPanel {



    JsonArray array = new JsonArray();
    JTabbedPane tabbedPane;

    public DialogWithTabs() {
        super(new GridLayout(1, 1));

        tabbedPane = new JTabbedPane();

        //Add the tabbed pane to this panel.
        add(tabbedPane);

        //The following line enables to use scrolling tabs.
        tabbedPane.setTabLayoutPolicy(JTabbedPane.SCROLL_TAB_LAYOUT);
    }


    /**
     * Create the GUI and show it.  For thread safety,
     * this method should be invoked from
     * the event dispatch thread.
     */
    public void createAndShowGUI() {
        //Schedule a job for the event dispatch thread:
        //creating and showing this application's GUI.
        SwingUtilities.invokeLater(new Runnable()
        {
            @Override
            public void run()
            {
                try {
                    UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
                } catch (Exception e) {
                    e.printStackTrace();
                }
                //Create and set up the window.
                String TITLE1 = "PhyloNetHMM Output Information";

                JFrame frame = new JFrame(TITLE1);
                frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);

                JMenuBar bar = createMenuBar();
                frame.setJMenuBar(bar);

                //Add content to the window.
                frame.add(DialogWithTabs.this, BorderLayout.CENTER);

                //Display the window.
                frame.pack();
                frame.setVisible(true);

            }
        });
    }

    /**
     * Creates a menu bar with a save button.
     *
     * @return The JMenuBar created.
     */
    private JMenuBar createMenuBar() {
        JMenuBar bar = new JMenuBar();

        JMenu menu = new JMenu("SaveMenu");
        bar.add(menu);

        JMenuItem saveButton = new JMenuItem("SaveFile");
        saveButton.addActionListener(new ActionListener()
        {
            @Override
            public void actionPerformed(ActionEvent e)
            {
                saveData();
            }
        });

        menu.add(saveButton);
        return bar;
    }

    /**
     * Saves the data in the tabs by opening up a filechooser and then saving to the selected file.
     */
    private void saveData() {
        File saveFile = selectSaveFile();
        if (saveFile != null)
            writeDataToFile(saveFile);
    }

    /**
     * Selects a file that will be used to save the data.
     *
     * @return The selected file, or empty.
     */
    private File selectSaveFile() {
        JFileChooser chooser = new JFileChooser();
        int choice = chooser.showDialog(this, "Save");

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
     * Writes the current data to the given file.
     *
     * @param fileChoice The file to write to. Will by closed by this method.
     */
    public void writeDataToFile(File fileChoice) {
        try (FileWriter write = new FileWriter(fileChoice)) {
            write.write(array.toString());
        } catch (IOException e1) {
            e1.printStackTrace();
        }
    }

    /**
     * Adds a tab to the view as well as stores it for later saving.
     *
     * @param name The name for the tab.
     * @param tab  The tab to add.
     */
    public void addPanel(String name, SavableTabData tab) {
        JsonObject tabObject = new JsonObject();
        tabObject.add("name", new JsonPrimitive(name));
        tabObject.add("tab", tab.save());

        array.add(tabObject);
        tabbedPane.addTab(name, tab.getComponent());
    }


}
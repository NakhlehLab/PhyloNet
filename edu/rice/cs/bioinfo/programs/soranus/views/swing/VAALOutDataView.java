package edu.rice.cs.bioinfo.programs.soranus.views.swing;

import edu.rice.cs.bioinfo.programs.soranus.viewModels.VAALOutDataVM;

import javax.swing.*;
import java.awt.*;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 3/28/13
 * Time: 4:59 PM
 * To change this template use File | Settings | File Templates.
 */
public class VAALOutDataView extends JPanel {

    public VAALOutDataView(VAALOutDataVM vm)
    {
        this.setLayout(new BorderLayout());

        JTextPane textPane = new JTextPane();
        textPane.setText(vm.Content);
        textPane.setEditable(false);

        JScrollPane scrollPane = new JScrollPane(textPane);
        textPane.setCaretPosition(0);


        this.add(scrollPane);

    }
}
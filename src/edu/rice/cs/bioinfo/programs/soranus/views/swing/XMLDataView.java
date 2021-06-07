package edu.rice.cs.bioinfo.programs.soranus.views.swing;

import edu.rice.cs.bioinfo.programs.soranus.viewModels.XMLDataVM;

import javax.swing.*;
import java.awt.*;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 2/8/13
 * Time: 2:25 PM
 * To change this template use File | Settings | File Templates.
 */
public class XMLDataView extends JPanel {

    public XMLDataView(XMLDataVM xmlDataVM)
    {
        this.setLayout(new BorderLayout());

        JTextPane textPane = new JTextPane();
        textPane.setText(xmlDataVM.Content);
        textPane.setEditable(false);

        JScrollPane scrollPane = new JScrollPane(textPane);
        textPane.setCaretPosition(0);


        this.add(scrollPane);

    }
}

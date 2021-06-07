package edu.rice.cs.bioinfo.programs.soranus.views.swing;

import edu.rice.cs.bioinfo.programs.soranus.viewModels.RecombResultVM;

import javax.swing.*;
import java.awt.*;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 4/25/13
 * Time: 5:15 PM
 * To change this template use File | Settings | File Templates.
 */
public class RecombResultView extends JPanel
{
    private RecombResultVM _vm;

    public RecombResultView(RecombResultVM vm)
    {
        _vm = vm;

        JLabel htmlContent;
        if(vm != null)
        {

            htmlContent = new JLabel("<html><table>" +
                    "<tr><td>"                    + "</td><td>" + vm.Site1Label +     "</td><td>" + vm.Site2Label     + "</td></tr>" +
                    "<tr><td>"+ vm.Sequence1Label + "</td><td>" + vm.Sequence0Site0 + "</td><td>" + vm.Sequence0Site1 + "</td></tr>" +
                    "<tr><td>"+ vm.Sequence2Label + "</td><td>" + vm.Sequence1Site0 + "</td><td>" + vm.Sequence1Site1 + "</td></tr>" +
                    "<tr><td>"+ vm.Sequence3Label + "</td><td>" + vm.Sequence2Site0 + "</td><td>" + vm.Sequence2Site1 + "</td></tr>" +
                    "<tr><td>"+ vm.Sequence4Label + "</td><td>" + vm.Sequence3Site0 + "</td><td>" + vm.Sequence3Site1 + "</td></tr>" +
                                     "</table></html>");
        }
        else
        {
            htmlContent = new JLabel("No recombination detected.");
        }

        htmlContent.setBackground(Color.WHITE);
        this.add(htmlContent);
    }
}

package edu.rice.cs.bioinfo.programs.soranus.views.swing;

import edu.rice.cs.bioinfo.library.language.dot_2013_1.printing.DigraphDotPrinter;
import edu.rice.cs.bioinfo.programs.soranus.viewModels.TransMapVM;

import javax.swing.*;
import java.awt.*;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 4/2/13
 * Time: 4:31 PM
 * To change this template use File | Settings | File Templates.
 */
public class TransMapViewDot<N,E> extends TransMapViewBase<N,E>
{
    public TransMapViewDot(TransMapVM<N, E> transMapVM) {
        super(transMapVM);

        this.setLayout(new BorderLayout());

        String dotContent = generateDot(transMapVM);

        JTextPane textPane = new JTextPane();
        textPane.setText(dotContent);
        textPane.setEditable(false);

        JScrollPane scrollPane = new JScrollPane(textPane);
        textPane.setCaretPosition(0);


        this.add(scrollPane);
    }

    private String generateDot(final TransMapVM<N, E> transMapVM) {

        return new DigraphDotPrinter<N,E>()
        {
            @Override
            protected N getSource(E edge) {
                return transMapVM.getSource(edge);
            }

            @Override
            protected N getDestination(E edge) {
                return transMapVM.getDestination(edge);
            }

            @Override
            protected String getNodeLabel(N node)
            {
                return transMapVM.getNodeLabel(node);
            }


        }.toDot(transMapVM.getNodes(), transMapVM.getEdges());

    }
}

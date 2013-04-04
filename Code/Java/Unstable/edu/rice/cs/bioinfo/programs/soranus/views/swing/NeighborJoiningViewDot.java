package edu.rice.cs.bioinfo.programs.soranus.views.swing;

import edu.rice.cs.bioinfo.library.language.dot_2013_1.printing.DigraphDotPrinter;
import edu.rice.cs.bioinfo.programs.soranus.viewModels.NeighborJoiningVM;

import javax.swing.*;
import java.awt.*;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 4/1/13
 * Time: 12:52 PM
 * To change this template use File | Settings | File Templates.
 */
public class NeighborJoiningViewDot<N,E> extends JPanel {

    public NeighborJoiningViewDot(NeighborJoiningVM<N,E> vm)
    {
        this.setLayout(new BorderLayout());

        String dotContent = generateDot(vm);

        JTextPane textPane = new JTextPane();
        textPane.setText(dotContent);
        textPane.setEditable(false);

        JScrollPane scrollPane = new JScrollPane(textPane);
        textPane.setCaretPosition(0);


        this.add(scrollPane);

    }

    private String generateDot(final NeighborJoiningVM<N,E> vm) {

        return new DigraphDotPrinter<N,E>()
        {
            @Override
            protected N getSource(E edge) {
                return vm.getNodesOfEdge(edge).Item1;
            }

            @Override
            protected N getDestination(E edge) {
                return vm.getNodesOfEdge(edge).Item2;
            }

            @Override
            protected String getEdgeLabel(E edge) {
                return vm.getEdgeLabel(edge);
            }

            @Override
            protected String getNodeLabel(N node)
            {
                return vm.getNodeLabel(node);
            }


        }.toDot(vm.Nodes, vm.Edges);
    }
}

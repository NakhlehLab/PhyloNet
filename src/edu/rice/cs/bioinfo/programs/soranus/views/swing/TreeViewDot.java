package edu.rice.cs.bioinfo.programs.soranus.views.swing;

import edu.rice.cs.bioinfo.library.language.dot_2013_1.printing.GraphDotPrinter;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.soranus.viewModels.TreeVM;

import javax.swing.*;
import java.awt.*;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 4/5/13
 * Time: 2:09 PM
 * To change this template use File | Settings | File Templates.
 */
public class TreeViewDot<N,E> extends JPanel
{
    public TreeViewDot(TreeVM<N,E> vm)
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

    private String generateDot(final TreeVM<N,E> vm) {

        return new GraphDotPrinter<N,E>()
        {
            {
                this.setLegend(vm.Legend);
            }

            @Override
            protected N getEdgeRhs(E edge)
            {
                return vm.getNodesOfEdge(edge).Item2;
            }

            @Override
            protected N getEdgeLhs(E edge)
            {
                return vm.getNodesOfEdge(edge).Item1;
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


            @Override
            protected Tuple<N, N> getNodesOfEdge(E edge)
            {
                return vm.getNodesOfEdge(edge);  //To change body of implemented methods use File | Settings | File Templates.
            }

            @Override
            protected String getEdgeColor(E edge)
            {
                int edgeColor = vm.getEdgeColor(edge);

                switch (edgeColor)
                {
                    case 0:
                        return "black";
                    case 1:
                        return "blue";
                    case 2:
                        return "red";
                    default:
                        return null;

                }


            }
        }.toDot(vm.getEdges());
    }
}

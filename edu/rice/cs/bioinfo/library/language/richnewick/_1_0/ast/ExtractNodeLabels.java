package edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 10/19/11
 * Time: 3:08 PM
 * To change this template use File | Settings | File Templates.
 */
public class ExtractNodeLabels implements NetworkAlgo<Collection<String>, Object, RuntimeException> {

    public ExtractNodeLabels()
    {

    }

    public Collection<String> forNetworkEmpty(NetworkEmpty network, Object input) throws RuntimeException {
        return new ArrayList<String>();
    }

    public Collection<String> forNetworkNonEmpty(NetworkNonEmpty network, Object input) throws RuntimeException {

        Collection<String> labels = new LinkedList<String>();
        collectHelp(network.PrincipleInfo, network.PrincipleDescendants, labels);

        return labels;
    }

    private void collectHelp(NetworkInfo node, DescendantList children, final Collection<String> labels) {

        node.NodeLabel.execute(new NodeLabelAlgo<Object, Object, RuntimeException>() {
            public Object forNodeLabelNonEmpty(NodeLabelNonEmpty node, Object input) throws RuntimeException {
                labels.add(node.Label.Content);
                return null;  //To change body of implemented methods use File | Settings | File Templates.
            }

            public Object forNodeLabelEmpty(NodeLabelEmpty node, Object input) throws RuntimeException {
                return null;  //To change body of implemented methods use File | Settings | File Templates.
            }
        }, null);

        for(Subtree tree : children.Subtrees)
        {
            collectHelp(tree.NetworkInfo, tree.Descendants, labels);
        }

    }
}

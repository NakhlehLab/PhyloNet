package edu.rice.cs.bioinfo.library.language.richnewick._1_0.printing.ast;

import edu.rice.cs.bioinfo.library.language.richnewick._1_0.printing.HybridNodeType;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.printing.RichNewickPrinterCompact;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.*;
import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Func2;
import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;

import java.io.StringWriter;
import java.util.List;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 2/27/13
 * Time: 2:25 PM
 * To change this template use File | Settings | File Templates.
 */
public class AstRichNewickPrinterCompact extends RichNewickPrinterCompact<Subtree>
{
    private static final Func1<Subtree,String> _GET_LABEL = new Func1<Subtree, String>() {
        public String execute(Subtree node) {
            return node.NetworkInfo.NodeLabel.execute(new NodeLabelAlgo<String, Void, RuntimeException>() {
                public String forNodeLabelNonEmpty(NodeLabelNonEmpty node, Void input) throws RuntimeException {
                    return node.Label.Content;
                }

                public String forNodeLabelEmpty(NodeLabelEmpty node, Void input) throws RuntimeException {
                    return "";
                }
            }, null);
        }
    };

    private static final Func1<Subtree,Iterable<Subtree>> _GET_DESTINATION_NODES = new Func1<Subtree, Iterable<Subtree>>() {
        public Iterable<Subtree> execute(Subtree node)
        {
            return node.Descendants.Subtrees;
        }
    };

    private static final Func1<Subtree,String> _GET_HYBRID_INDEX = new Func1<Subtree, String>() {
        public String execute(Subtree node) {
            return node.NetworkInfo.HybridNodeQualifier.execute(new HybridNodeQualifierAlgo<String, Void, RuntimeException>() {

                public String forHybridNodeQualifierEmpty(HybridNodeQualifierEmpty qualifier, Void input) throws RuntimeException {
                    return null;
                }

                public String forHybridNodeQualifierNonEmpty(HybridNodeQualifierNonEmpty qualifier, Void input) throws RuntimeException {
                    return qualifier.HybridNodeIndex.Content;
                }

                public String forHybridNodeQualifierWithType(HybridNodeQualifierWithType qualifier, Void input) throws RuntimeException {
                    return qualifier.HybridNodeIndex.Content;
                }
            }, null);
        }
    };

    public AstRichNewickPrinterCompact()
    {
        this.setGetBranchLength(new Func2<Subtree, Subtree, String>() {
            public String execute(Subtree node1, Subtree node2) {

                List<Subtree> node1Desc = IterableHelp.toList(node1.Descendants.Subtrees);
                List<Subtree> node2Desc = IterableHelp.toList(node2.Descendants.Subtrees);

                Subtree blContainer =  node1Desc.contains(node2) ? node2 :
                        node2Desc.contains(node1) ? node1 : null;

                if(blContainer != null)
                {
                    return blContainer.NetworkInfo.BranchLength.execute(new BranchLengthAlgo<String, Void, RuntimeException>()
                    {
                        public String forBranchLengthEmpty(BranchLengthEmpty branchLength, Void input) throws RuntimeException {
                            return null;
                        }

                        public String forBranchLengthNonEmpty(BranchLengthNonEmpty branchLength, Void input) throws RuntimeException {
                            return branchLength.Length.Content;
                        }
                    }, null);
                }
                else
                {
                    return null;
                }
            }
        });

        this.setGetProbability(new Func2<Subtree, Subtree, String>() {
            public String execute(Subtree node1, Subtree node2) {

                List<Subtree> node1Desc = IterableHelp.toList(node1.Descendants.Subtrees);
                List<Subtree> node2Desc = IterableHelp.toList(node2.Descendants.Subtrees);

                Subtree probContainer =  node1Desc.contains(node2) ? node2 :
                        node2Desc.contains(node1) ? node1 : null;

                if(probContainer != null)
                {
                    return probContainer.NetworkInfo.Probability.execute(new ProbabilityAlgo<String, Object, RuntimeException>()
                    {
                        public String forProbabilityEmpty(ProbabilityEmpty prob, Object input) {
                            return null;
                        }

                        public String forProbabilityNonEmpty(ProbabilityNonEmpty prob, Object input) {
                            return prob.ProbabilityValue.Content;
                        }
                    }, null);
                }
                else
                {
                    return null;
                }
            }
        });

        this.setGetSupport(new Func2<Subtree, Subtree, String>() {
            public String execute(Subtree node1, Subtree node2) {

                List<Subtree> node1Desc = IterableHelp.toList(node1.Descendants.Subtrees);
                List<Subtree> node2Desc = IterableHelp.toList(node2.Descendants.Subtrees);

                Subtree supportContainer =  node1Desc.contains(node2) ? node2 :
                        node2Desc.contains(node1) ? node1 : null;

                if(supportContainer != null)
                {
                    return supportContainer.NetworkInfo.Support.execute(new SupportAlgo<String, Object, RuntimeException>()
                    {
                        public String forSupportNonEmpty(SupportNonEmpty support, Object input) {
                            return support.SupportValue.Content;
                        }

                        public String forSupportEmpty(SupportEmpty support, Object input) {
                            return null;
                        }
                    }, null);
                }
                else
                {
                    return null;
                }
            }
        });

    }

    private static final Func1<Subtree,HybridNodeType> _GET_HYBRID_TYPE = new Func1<Subtree, HybridNodeType>() {
        public HybridNodeType execute(Subtree node) {
            return node.NetworkInfo.HybridNodeQualifier.execute(new HybridNodeQualifierAlgo<HybridNodeType, Void, RuntimeException>() {

                public HybridNodeType forHybridNodeQualifierEmpty(HybridNodeQualifierEmpty qualifier, Void input) throws RuntimeException {
                    return null;
                }

                public HybridNodeType forHybridNodeQualifierNonEmpty(HybridNodeQualifierNonEmpty qualifier, Void input) throws RuntimeException {
                    return null;
                }

                public HybridNodeType forHybridNodeQualifierWithType(HybridNodeQualifierWithType qualifier, Void input) throws RuntimeException {

                    edu.rice.cs.bioinfo.library.language.richnewick._1_0.HybridNodeType type =
                            edu.rice.cs.bioinfo.library.language.richnewick._1_0.HybridNodeType.fromString(qualifier.HybridNodeType.Content);

                    if(type == edu.rice.cs.bioinfo.library.language.richnewick._1_0.HybridNodeType.Recombination)
                    {
                        return HybridNodeType.Recombination;
                    }
                    else if(type == edu.rice.cs.bioinfo.library.language.richnewick._1_0.HybridNodeType.Hybridization)
                    {
                        return HybridNodeType.Hybridization;
                    }
                    else if(type == edu.rice.cs.bioinfo.library.language.richnewick._1_0.HybridNodeType.LateralGeneTransfer)
                    {
                        return HybridNodeType.LateralGeneTransfer;
                    }

                    throw new IllegalStateException("Unknown hybrid type");

                }
            }, null);
        }
    };

    public void print(NetworkNonEmpty network, StringWriter writer)
    {
        boolean isRooted = network.RootageQualifier.execute(new IsRooted(), null);
        Subtree printRoot = new Subtree(network.PrincipleDescendants, network.PrincipleInfo);

        print(isRooted,printRoot, _GET_LABEL, _GET_DESTINATION_NODES, _GET_HYBRID_INDEX, _GET_HYBRID_TYPE, writer);


    }
}

package edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents;

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ast.*;
import org.w3c.dom.traversal.NodeIterator;

import java.io.FileInputStream;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/6/11
 * Time: 2:52 PM
 * To change this template use File | Settings | File Templates.
 */
public class BlockContentsFactoryFromAST {

    public static BlockContents make(Blocks blocks)
    {
        final Map<String,RichNewickAssignment> identToRichNewickAssignment = new HashMap<String,RichNewickAssignment>();
        final LinkedList<SyntaxCommand> commands = new LinkedList<SyntaxCommand>();


        for(Block block : blocks.Contents)
        {
            block.execute(new BlockAlgo<Object, Object, RuntimeException>() {
                public Object forTreesBlock(TreesBlockBody treesBlock, Object input) throws RuntimeException {
                   forRNewichAssignmentsBlockBody(treesBlock, true);
                   return null;
                }

                 public Object forNetworksBlock(NetworksBlockBody networksBlock, Object input) throws RuntimeException {

                   forRNewichAssignmentsBlockBody(networksBlock, false);
                   return null;
                }

                private <T extends RNewickAssigmentContainer> void forRNewichAssignmentsBlockBody(RNewichAssignmentsBlockBody<T> body, final boolean fromTreesBlock)
                {
                    for(RNewickAssigmentContainer assignment : body.getAssignments())
                    {
                        final edu.rice.cs.bioinfo.library.language.pyson._1_0.ast.RichNewickAssignment rAssignment = assignment.getAssignment();
                        final String assignmentIdentifier = rAssignment.LHSIdentifier.Content;
                        identToRichNewickAssignment.put(assignmentIdentifier, new RichNewickAssignment()
                        {
                            public String getIdentifier() {
                                return rAssignment.LHSIdentifier.Content;  //To change body of implemented methods use File | Settings | File Templates.
                            }

                            public boolean isDefinedByNetworksBlock() {
                                return !fromTreesBlock;  //To change body of implemented methods use File | Settings | File Templates.
                            }

                            public boolean isDefinedByTreesBlock() {
                                return fromTreesBlock;  //To change body of implemented methods use File | Settings | File Templates.
                            }

                            public String getRichNewickString() {
                                return rAssignment.RHSRichNewickString.String;
                            }

                            public int getRichNewickStringLine() {
                                return rAssignment.RHSRichNewickString.LineNumber;
                            }

                            public int getRichNewickStringColumn() {
                                return rAssignment.RHSRichNewickString.ColumnNumber;
                            }
                        });


                    }

                }

                public Object forPhylonetBlockBody(PhyloNetBlockBody phyloBlock, Object input) throws RuntimeException {

                    for(PhyloNetCommand command : phyloBlock.Commands)
                    {
                        PhyloNetCommandPart commandName = command.Parts.iterator().next();
                        final String nameFinal = commandName.execute(new PhyloNetCommandPartAlgo<String, Object, RuntimeException>() {
                            public String forIdentifier(PhyloNetCommandPartIdentifier ident, Object input) throws RuntimeException {
                               return ident.IdentContents;
                            }

                            public String forIdenList(PhyloNetCommandPartIdentList ident, Object input) throws RuntimeException {
                                return null;  //To change body of implemented methods use File | Settings | File Templates.
                            }

                            public String forQuote(PhyloNetCommandPartQuote quote, Object input) throws RuntimeException {
                                return null;  //To change body of implemented methods use File | Settings | File Templates.
                            }

                            public String forSetList(PhyloNetCommandPartSetList setList, Object input) throws RuntimeException {
                                return null;  //To change body of implemented methods use File | Settings | File Templates.
                            }

                            public String forIdentSet(PhyloNetCommandPartIdentSet identSet, Object input) throws RuntimeException {
                                return null;  //To change body of implemented methods use File | Settings | File Templates.
                            }

                            public String forTaxaMap(PhyloNetCommandPartTaxaMap taxaMap, Object input) throws RuntimeException {
                                return null;  //To change body of implemented methods use File | Settings | File Templates.
                            }
                        }, null);
                        final int lineFinal = commandName.Line;
                        final int colFinal = commandName.Col;


                        final LinkedList<Parameter> commandArguments = makeCommandArguments(command);


                        commands.add(new SyntaxCommand() {
                            public String getName() {
                                return nameFinal;
                            }

                             public int getLine()
                             {
                                 return lineFinal;
                             }

                             public int getColumn()
                             {
                                 return colFinal;
                             }

                            public Iterable<Parameter> getParameters() {
                                return commandArguments;
                            }
                        });
                    }

                    return null;  //To change body of implemented methods use File | Settings | File Templates.
                }
            }, null);
        }

        return new BlockContents() {

            public RichNewickAssignment getRichNewickAssigment(String identifier) {
                return identToRichNewickAssignment.get(identifier);
            }

            public Iterable<String> getRickNewickAssignmentIdentifiers()
            {
                return identToRichNewickAssignment.keySet();
            }

            public Iterable<SyntaxCommand> getCommands() {
                return commands;  //To change body of implemented methods use File | Settings | File Templates.
            }
        };
    }

    private static LinkedList<Parameter> makeCommandArguments(PhyloNetCommand command) {

        LinkedList<Parameter> commandArguments = new LinkedList<Parameter>();
        boolean first = true;

        for(PhyloNetCommandPart part : command.Parts)
        {
            if(first)
            {
                first = false;
                continue;
            }
            else
            {
                /*
                String partString = part.Content;
                if((partString.startsWith("\"") && partString.endsWith("\"")) ||
                        (partString.startsWith("\'") && partString.endsWith("\'")))
                {
                    partString = partString.substring(1, partString.length() - 1);
                }
                final String partStringFinal = partString;     */



                final int lineFinal = part.Line;
                final int colFinal = part.Col;

                Parameter p = part.execute(new PhyloNetCommandPartAlgo<Parameter, Object, RuntimeException>() {
                    public Parameter forIdentifier(PhyloNetCommandPartIdentifier ident, Object input) throws RuntimeException {

                        return new ParameterIdent(lineFinal,  colFinal, ident.IdentContents);

                    }

                    public Parameter forQuote(PhyloNetCommandPartQuote quote, Object input) throws RuntimeException {

                        return new ParameterQuote(lineFinal,  colFinal,quote.TotalText);
                    }

                    public Parameter forSetList(PhyloNetCommandPartSetList setList, Object input) throws RuntimeException {
                        return new ParameterTaxonSetList(lineFinal, colFinal, setList.PartContents);
                    }

                    public Parameter forIdentSet(PhyloNetCommandPartIdentSet set, Object input) throws RuntimeException {
                        return new ParameterIdentSet(lineFinal, colFinal, set.PartContents);
                    }

                    public Parameter forTaxaMap(PhyloNetCommandPartTaxaMap taxaMap, Object input) throws RuntimeException {
                        return new ParameterTaxaMap(lineFinal, colFinal, taxaMap.Map);
                    }

                      public Parameter forIdenList(PhyloNetCommandPartIdentList ident, Object input) throws RuntimeException {

                          LinkedList<String> elements = new LinkedList<String>();

                          for(Identifier identElement : ident.List.Elements)
                          {
                                elements.add(identElement.Content);
                          }

                          return new ParameterIdentList(lineFinal, colFinal, elements);
                    }
                }, null);


                commandArguments.add(p);
            }
        }

        return commandArguments;
    }
}

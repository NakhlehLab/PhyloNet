package edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.keyedstringsandcommands;

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ast.*;

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
public class SACFactoryFromAST {

    public static StringsAndCommands make(Blocks blocks)
    {
        final Map<String,String> identToRichNewickString = new HashMap<String,String>();
        final LinkedList<SyntaxCommand> commands = new LinkedList<SyntaxCommand>();


        for(Block block : blocks.Contents)
        {
            block.execute(new BlockAlgo<Object, Object, RuntimeException>() {
                public Object forTreesBlock(TreesBlockBody treesBlock, Object input) throws RuntimeException {

                    for(TreeAssignment assignment : treesBlock.getAssignments())
                    {
                        identToRichNewickString.put(assignment.Assignment.LHSIdentifier.Content,
                                                    assignment.Assignment.RHSRichNewickString.String);
                    }

                    return null;  //To change body of implemented methods use File | Settings | File Templates.
                }

                public Object forPhylonetBlockBody(PhyloNetBlockBody phyloBlock, Object input) throws RuntimeException {

                    for(PhyloNetCommand command : phyloBlock.Commands)
                    {
                        boolean first = true;
                        String name = null;
                        int line = -1, col = -1;


                        final LinkedList<Parameter> commandArguments = new LinkedList<Parameter>();

                        for(PhyloNetCommandPart part : command.Parts)
                        {
                            if(first)
                            {
                                name = part.Content;
                                line = part.Line;
                                col = part.Col;
                                first = false;
                            }
                            else
                            {
                                String partString = part.Content;
                                if((partString.startsWith("\"") && partString.endsWith("\"")) ||
                                   (partString.startsWith("\'") && partString.endsWith("\'")))
                                {
                                    partString = partString.substring(1, partString.length() - 1);
                                }
                                final String partStringFinal = partString;

                                final int lineFinal = part.Line;
                                final int colFinal = part.Col;
                                commandArguments.add(new Parameter() {
                                    public int getLine() {
                                        return lineFinal;  //To change body of implemented methods use File | Settings | File Templates.
                                    }

                                    public int getColumn() {
                                        return colFinal;  //To change body of implemented methods use File | Settings | File Templates.
                                    }

                                    public String getValue() {
                                        return partStringFinal;  //To change body of implemented methods use File | Settings | File Templates.
                                    }
                                });
                            }
                        }

                        final String nameFinal = name;
                        final int lineFinal = line;
                        final int colFinal = col;
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

        return new StringsAndCommands() {
            public String getRickNewickString(String identifier) {
                return identToRichNewickString.get(identifier);
            }

            public Iterable<String> getRickNewickStringIdentifiers()
            {
                return identToRichNewickString.keySet();
            }

            public Iterable<SyntaxCommand> getCommands() {
                return commands;  //To change body of implemented methods use File | Settings | File Templates.
            }
        };
    }
}

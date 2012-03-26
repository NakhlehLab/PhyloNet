/*
 * Copyright (c) 2012 Rice University.
 *
 * This file is part of PhyloNet.
 *
 * PhyloNet is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PhyloNet is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PhyloNet.  If not, see <http://www.gnu.org/licenses/>.
 */

package edu.rice.cs.bioinfo.programs.phylonet;

import edu.rice.cs.bioinfo.library.language.parsing.CoordinateParseError;
import edu.rice.cs.bioinfo.library.language.parsing.CoordinateParseErrorsException;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ast.Blocks;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.*;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.RichNewickReadResult;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast.Network;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast.Networks;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.csa.CSAError;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.parsers.antlr.ast.RichNewickParser;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.parsers.antlr.ast.RichNewickReaderAST_ANTLR;
import edu.rice.cs.bioinfo.library.programming.*;
import edu.rice.cs.bioinfo.programs.phylonet.commands.Command;
import edu.rice.cs.bioinfo.programs.phylonet.commands.CommandFactory;
import edu.rice.cs.bioinfo.programs.phylonet.commands.NexusOut;

import java.io.*;
import java.security.KeyStore;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/6/11
 * Time: 1:22 PM
 * To change this template use File | Settings | File Templates.
 */
public class Program {



    public static void main(String[] args) throws FileNotFoundException, IOException
    {
        System.out.print("");

        if(args.length !=1 && args.length != 2) // we only expect one or two parameters, the input nexus file and a random seed
        {
            showUsage();
            return;
        }

        File nexusFile = new File(args[0]);

        Random rand;
        if(args.length == 2)
        {
            try
            {
                long seed = Long.parseLong(args[1]);
                rand = new Random(seed);
            }
            catch(NumberFormatException e)
            {
                System.err.println("Unknown random seed: " +args[1]);
                showUsage();
                return;
            }
        }
        else
        {
            rand = new Random();
        }


        if(!nexusFile.isFile()) // assert the file the user gave us actually exists
        {
            showFileDoesNotExist(nexusFile);
        }
        else
        {
            run(new FileInputStream(nexusFile), System.err, System.out, rand);
        }

        System.out.print("\n");

    }

    static void run(InputStream nexusStream, final PrintStream errorStream,
                                             final PrintStream displaySteam, Random rand) throws IOException
    {
             Proc<String> display = new Proc<String>()
             {
                public void execute(String s) {

                     displaySteam.print(s);
                }
             };

             final Container<Boolean> allowCommandExecution = new Container<Boolean>(true);
             final Proc3<String, Integer, Integer> errorDetected = new Proc3<String, Integer, Integer>()
             {
                public void execute(String message, Integer line, Integer col) {

                    int oneBasedColumn = col + 1;
                    errorStream.print(String.format("\n\nError at [%s,%s]: %s", line, oneBasedColumn, message));
                    allowCommandExecution.setContents(false);
                }
             };

        /*
         * Parse input nexus file to AST. Report and terminate on syntax errors.
         */
        Blocks blocks;
        try
        {
            blocks = parseToBlocks(nexusStream); // Exception thrown on syntax errors.
        }
        catch(CoordinateParseErrorsException e)  // syntax errors detected
        {
            for(CoordinateParseError error : e.Errors)
            {
                errorDetected.execute(error.getMessage(), error.getLineNumber(), error.getColumnNumber() );
            }
            return;
        }

        // no syntax errors at this point

        /*
         * Perform Context Sensitive Analysis
         */
        boolean containsDuplicates = ContextSensitiveAnalyser.checkforDuplicateAssignmentIdentifiers(blocks, errorDetected);

        // if we already have errors, don't continue.  Else we will get confusing cascades.
        if(containsDuplicates)
        {
            return;
        }

        BlockContents blockContents = BlockContentsFactoryFromAST.make(blocks); // covert AST to our IR


        Map<String,NetworkNonEmpty> sourceIdentToNetwork = makeNetworks(blockContents, errorDetected); // make a Network representation of each defined Rich Newick string
                                                                                        // map is keyed by source code identifier of the string

        // if we already have errors, don't continue.  Else we will get confusing cascades.
        if(!allowCommandExecution.getContents())
        {
            return;
        }

        ContextSensitiveAnalyser.checkforHybridNodesInTrees(sourceIdentToNetwork, blockContents, errorDetected);
        ContextSensitiveAnalyser.analyseNetworks(sourceIdentToNetwork, blockContents, errorDetected); // check the Networks for any context errors

        LinkedList<Command> commands = new LinkedList<Command>();

        for(SyntaxCommand sCommand : blockContents.getCommands())
        {
            try
            {
                commands.add(CommandFactory.make(sCommand, sourceIdentToNetwork, errorDetected, rand));
            }
            catch(IllegalArgumentException e)
            {
                errorDetected.execute(e.getMessage(), sCommand.getLine(), sCommand.getColumn());
            }
        }

        for(Command command : commands)
        {
            command.checkParams();
        }

        /*
         * Execute commands if no errors detected
         */

        boolean first = true;
         BufferedWriter nexusOut = null;
        if(allowCommandExecution.getContents())
        {
           int commandNumber = 1;
           for(final Command command : commands)
           {
               final int commandNumberCapture = commandNumber;
               final BufferedWriter nexusOutCapture = nexusOut;
               if(commandNumber != 1)
               {
                  display.execute("\n");

               }
               try
               {
                    command.addSTTreeGeneratedListener(new Proc<String>()
                    {
                        private int _treeNumber = 1;

                        public void execute(String newickTree) {

                            if(nexusOutCapture != null)
                            {
                                try
                                {
                                    if(commandNumberCapture == 2)
                                    {
                                       nexusOutCapture.write("\n");
                                    }
                                    nexusOutCapture.write("\n" + commandNumberCapture + "_" + command.getDefiningSyntaxCommand().getName() + "_" + _treeNumber + " = " + newickTree);
                                    _treeNumber++;
                                }
                                catch(IOException e)
                                {
                                   SyntaxCommand motivatingSyntaxCommand = command.getDefiningSyntaxCommand();
                                   errorDetected.execute("Error writing to nexus out file. (" + e.getMessage() + ")", motivatingSyntaxCommand.getLine(), motivatingSyntaxCommand.getColumn());
                                }
                            }
                        }
                    });
                    showCommand(command.getDefiningSyntaxCommand(), display);
                    command.executeCommand(display);

                   if(command instanceof NexusOut)
                   {
                      if(nexusOut != null)
                      {
                           nexusOut.write("\n\nEND;");
                           nexusOut.flush();
                           nexusOut.close();
                      }

                      nexusOut = new BufferedWriter(new FileWriter((((NexusOut)command).getNexusOutFile())));
                      nexusOut.write("#NEXUS");
                      nexusOut.write("\n\nBEGIN TREES;");
                   }
               }
               catch(IOException e)
               {
                   SyntaxCommand motivatingSyntaxCommand = command.getDefiningSyntaxCommand();
                   errorDetected.execute(String.format("Error executing command '%s' (%s).", motivatingSyntaxCommand.getName(), e.getMessage()),
                                          motivatingSyntaxCommand.getLine(), motivatingSyntaxCommand.getColumn());
               }


               commandNumber++;
           }
        }

        if(nexusOut != null)
        {
           nexusOut.write("\n\nEND;");
           nexusOut.flush();
           nexusOut.close();
        }
    }

    private static void showCommand(SyntaxCommand definingSyntaxCommand, Proc<String> displayResult) {

        StringBuffer accum = new StringBuffer("\n" + definingSyntaxCommand.getName());



        for(edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.Parameter p : definingSyntaxCommand.getParameters())
        {
              String paramValue = p.execute(new ParameterAlgo<String, Object, RuntimeException>() {

                  public String forIdentifier(ParameterIdent parameterIdent, Object o) throws RuntimeException {
                      return parameterIdent.Content;
                  }

                  public String forIdentList(ParameterIdentList parameterIdentList, Object o) throws RuntimeException {
                     StringBuilder b = new StringBuilder();
                     b.append("(");

                    Iterator<String> elements = parameterIdentList.Elements.iterator();

                    if(elements.hasNext())
                    {
                        b.append(elements.next());
                    }

                    while(elements.hasNext())
                    {
                        b.append(", " + elements.next());
                    }

                     b.append(")");

                      return b.toString();
                  }

                  public String forQuote(ParameterQuote parameterQuote, Object o) throws RuntimeException {
                      return "\"" + parameterQuote.UnquotedText + "\"";
                  }

                  public String forTaxonSetList(ParameterTaxonSetList parameterTaxonSetList, Object o) throws RuntimeException {
                      return parameterTaxonSetList.OriginalSource;
                  }

                  public String forIdentSet(ParameterIdentSet parameterIdentSet, Object o) throws RuntimeException {
                      return parameterIdentSet.OriginalSource;
                  }

                  public String forTaxaMap(ParameterTaxaMap parameterTaxaMap, Object o) throws RuntimeException {

                      StringBuilder b = new StringBuilder();
                      b.append("<");

                      boolean firstKey = true;
                      for(Map.Entry<String,List<String>> entry : parameterTaxaMap._mappings)
                      {
                          if(!firstKey)
                          {
                              b.append("; ");

                          }

                           b.append(entry.getKey() + ":");
                           firstKey = false;


                          boolean firstValueEntry = true;
                          for(String mapEntry : entry.getValue())
                          {
                             if(!firstValueEntry)
                             {
                               b.append(",");
                             }

                             b.append(mapEntry);
                             firstValueEntry = false;
                          }
                      }

                      b.append(">");

                      return b.toString();

                  }
              }, null);
              accum.append(" " + paramValue);
        }

        displayResult.execute(accum.toString());
    }


    private static Map<String, NetworkNonEmpty> makeNetworks(BlockContents blockContents, Proc3<String, Integer, Integer> errorDetected) throws IOException {

        HashMap<String,NetworkNonEmpty> tbr = new HashMap<String, NetworkNonEmpty>();
        for(String richNewickSourceIdent : blockContents.getRickNewickAssignmentIdentifiers())
        {
            RichNewickAssignment assignment = blockContents.getRichNewickAssigment(richNewickSourceIdent);


            RichNewickReadResult<Networks> readResult = null;
            try
            {
               RichNewickReaderAST_ANTLR reader = new RichNewickReaderAST_ANTLR();
               readResult = reader.read(new ByteArrayInputStream(assignment.getRichNewickString().getBytes()));
             /*  Iterator<NetworkNonEmpty> oneNetwork = RichNewickParser.parse(
                        new ByteArrayInputStream(assignment.getRichNewickString().getBytes())).Networks.iterator();  */
            }
            catch(CoordinateParseErrorsException e)
            {
                for(CoordinateParseError error : e.Errors)
                {
                    errorDetected.execute(error.getMessage(), assignment.getRichNewickStringLine() + error.getLineNumber() -1, assignment.getRichNewickStringColumn() + error.getColumnNumber());
                }
                return new HashMap<String, NetworkNonEmpty>();
            }
            catch(RuntimeException e)
            {
                errorDetected.execute("Invalid Rich Newick string.", assignment.getRichNewickStringLine(), assignment.getRichNewickStringColumn());
                return new HashMap<String, NetworkNonEmpty>();
            }

            boolean sawError = false;
            for(CSAError error : readResult.getContextErrors())
            {
                errorDetected.execute(error.Message, assignment.getRichNewickStringLine() + error.LineNumber -1, assignment.getRichNewickStringColumn() + error.ColumnNumber);
                sawError = true;
            }

            if(sawError)
            {
                return new HashMap<String, NetworkNonEmpty>();
            }

            Iterator<NetworkNonEmpty> networks = readResult.getNetworks().Networks.iterator();
            if(!networks.hasNext())
            {
                errorDetected.execute(
                        String.format("Rich Newick string '%s' does not define a network.", richNewickSourceIdent),
                        assignment.getRichNewickStringLine(), assignment.getRichNewickStringColumn());
            }

            NetworkNonEmpty network = networks.next();

            if(networks.hasNext())
            {
                errorDetected.execute(
                        String.format("Rich Newick string '%s' defines more than one network.", richNewickSourceIdent),
                        assignment.getRichNewickStringLine(), assignment.getRichNewickStringColumn());
                continue;
            }

            tbr.put(richNewickSourceIdent, network);

        }

        return tbr;
    }

    private static Blocks parseToBlocks(InputStream nexusStream) throws IOException, CoordinateParseErrorsException {

        return edu.rice.cs.bioinfo.library.language.pyson._1_0.parsers.antlr.ast.Parser.parse(nexusStream);

    }

    private static void showFileDoesNotExist(File inFile) {
        System.err.println(String.format("\nNo such file '%s'", inFile.getName())); 
    }

    private static void showUsage() {

        System.out.println("\nUsage: java -jar phylonet.jar nexus_file [random seed integer]");

    }
}
